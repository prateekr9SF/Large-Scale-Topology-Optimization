// filter_sens_pthread.c
// Sensitivity filtering (no volume weighting) using pthreads + OpenMP atomics.
// Math: df/dx_i = sum_j [ H_{j i} / sum_k H_{j k} ] * (df/dx_tilde_j)
// Triplet files (text) read from CWD: drow.dat, dcol.dat, dval.dat
// Indices in files are assumed 1-based; converted to 0-based internally.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>

#ifdef _OPENMP
  #include <omp.h>   // for #pragma omp atomic (portable atomics on doubles)
#endif

#ifndef BLOCK_SIZE
#define BLOCK_SIZE (1<<20)   // ~1M triplets per block; tune for your I/O/cache
#endif

typedef struct {
    int thread_id, num_threads;
    int *drow, *dcol;          // buffered 1-based indices
    double *dval;              // buffered weights
    int block_read, ne;
    double q;                  // exponent on weights (use 1.0 typically)
    // PASS 1:
    double *row_sum;           // size ne: donor sums, row_sum[j] = sum_k H_{j,k}^q
    // PASS 2:
    const double *SensIn;      // size ne: df/dx_tilde (input)
    double *SensOut;           // size ne: df/dx        (output, accumulated)
} ThreadArgs;

// ---------- helpers ----------
static inline int inb(int x,int n){ return (unsigned)x < (unsigned)n; }

// Count nnz by scanning triplet text files in lockstep.
static long long count_nnz_from_files(void){
    FILE *fr=fopen("drow.dat","r");
    FILE *fc=fopen("dcol.dat","r");
    FILE *fv=fopen("dval.dat","r");
    if(!fr||!fc||!fv){
        fprintf(stderr,"ERROR: cannot open drow.dat/dcol.dat/dval.dat: %s\n", strerror(errno));
        if(fr) fclose(fr); if(fc) fclose(fc); if(fv) fclose(fv);
        return -1;
    }
    long long cnt=0; int r,c; double v;
    for(;;){
        int okr=fscanf(fr,"%d",&r);
        int okc=fscanf(fc,"%d",&c);
        int okv=fscanf(fv,"%lf",&v);
        if(okr==1 && okc==1 && okv==1) ++cnt; else break;
    }
    fclose(fr); fclose(fc); fclose(fv);
    return cnt;
}

// ---------- workers ----------
static void *worker_build_row_sums(void *args_ptr)
{
    ThreadArgs *a = (ThreadArgs*)args_ptr;
    int s = (a->block_read * a->thread_id) / a->num_threads;
    int e = (a->block_read * (a->thread_id + 1)) / a->num_threads;
    const int use_pow = (a->q != 1.0);

    for (int t = s; t < e; ++t) {
        int j = a->drow[t] - 1;    // donor row in H
        int i = a->dcol[t] - 1;    // receiver column (unused except bounds)
        if (!inb(j,a->ne) || !inb(i,a->ne)) continue;
        double w = use_pow ? pow(a->dval[t], a->q) : a->dval[t];
        #pragma omp atomic
        a->row_sum[j] += w;
    }
    return NULL;
}

static void *worker_accumulate_sens(void *args_ptr)
{
    ThreadArgs *a = (ThreadArgs*)args_ptr;
    int s = (a->block_read * a->thread_id) / a->num_threads;
    int e = (a->block_read * (a->thread_id + 1)) / a->num_threads;
    const int use_pow = (a->q != 1.0);

    for (int t = s; t < e; ++t) {
        int j = a->drow[t] - 1;    // donor j (filtered index)
        int i = a->dcol[t] - 1;    // design index i
        if (!inb(j,a->ne) || !inb(i,a->ne)) continue;

        double denom = a->row_sum[j];
        if (denom <= 0.0) continue;

        double w = use_pow ? pow(a->dval[t], a->q) : a->dval[t];
        double contrib = (w / denom) * a->SensIn[j];

        // Accumulate into column bin (dc/dx)[i]
        #pragma omp atomic
        a->SensOut[i] += contrib;
    }
    return NULL;
}

// ---------- public API ----------
void filterSensitivity_buffered_mt(const double *SensIn,  // df/dx_tilde[j]
                                   double *SensOut,       // df/dx[i] (output)
                                   int ne,
                                   long long nnz_total   // <=0 → auto-count
                                   )              // weight exponent
{

    double q = 1;
    if (ne <= 0) { fprintf(stderr,"ERROR: ne<=0\n"); exit(EXIT_FAILURE); }

    // threads from OMP_NUM_THREADS or default 4
    int num_threads = 4;
    const char *env = getenv("OMP_NUM_THREADS");
    if (env && *env) {
        int tmp = atoi(env);
        if (tmp > 0) num_threads = tmp;
    }
    printf("filterSensitivity_buffered_mt: using %d thread(s)\n", num_threads);

    if (nnz_total <= 0) {
        nnz_total = count_nnz_from_files();
        if (nnz_total <= 0) {
            fprintf(stderr,"ERROR: could not determine nnz_total from files.\n");
            exit(EXIT_FAILURE);
        }
    }

    // open triplets
    FILE *frow = fopen("drow.dat", "r");
    FILE *fcol = fopen("dcol.dat", "r");
    FILE *fval = fopen("dval.dat", "r");
    if (!frow || !fcol || !fval) {
        fprintf(stderr,"ERROR: opening triplet files: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    // allocate buffers/accumulators
    int *drow_block = (int*)malloc(BLOCK_SIZE * sizeof(int));
    int *dcol_block = (int*)malloc(BLOCK_SIZE * sizeof(int));
    double *dval_block = (double*)malloc(BLOCK_SIZE * sizeof(double));
    if (!drow_block || !dcol_block || !dval_block) {
        fprintf(stderr,"alloc BLOCK buffers failed\n"); exit(EXIT_FAILURE);
    }

    double *row_sum = (double*)calloc((size_t)ne, sizeof(double));
    if (!row_sum) { fprintf(stderr,"alloc row_sum failed\n"); exit(EXIT_FAILURE); }

    pthread_t *threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    ThreadArgs *targs  = (ThreadArgs*)malloc(num_threads * sizeof(ThreadArgs));
    if (!threads || !targs) { fprintf(stderr,"alloc thread structs failed\n"); exit(EXIT_FAILURE); }

    // ---------------- PASS 1: donor row sums ----------------
    {   
        printf("First pass for evaluating donor row sums...");
        long long total_read = 0;
        while (total_read < nnz_total) {
            int n = (int)((nnz_total - total_read) < BLOCK_SIZE ? (nnz_total - total_read) : BLOCK_SIZE);

            for (int i = 0; i < n; ++i) {
                if (fscanf(frow,"%d",&drow_block[i])!=1 ||
                    fscanf(fcol,"%d",&dcol_block[i])!=1 ||
                    fscanf(fval,"%lf",&dval_block[i])!=1) {
                    fprintf(stderr,"triplet read error (pass1) at %lld\n", total_read+i);
                    exit(EXIT_FAILURE);
                }
            }

            for (int t = 0; t < num_threads; ++t) {
                targs[t] = (ThreadArgs){
                    .thread_id=t, .num_threads=num_threads,
                    .drow=drow_block, .dcol=dcol_block, .dval=dval_block,
                    .block_read=n, .ne=ne, .q=q,
                    .row_sum=row_sum,
                    .SensIn=NULL, .SensOut=NULL
                };
                if (pthread_create(&threads[t], NULL, worker_build_row_sums, &targs[t]) != 0) {
                    perror("pthread_create pass1"); exit(EXIT_FAILURE);
                }
            }
            for (int t = 0; t < num_threads; ++t) pthread_join(threads[t], NULL);

            total_read += n;
        }

        printf("done! \n");
    }

    fclose(frow); fclose(fcol); fclose(fval);

    // ---------------- PASS 2: accumulate sensitivities ----------------
    for (int i = 0; i < ne; ++i) SensOut[i] = 0.0;

    frow = fopen("drow.dat", "r");
    fcol = fopen("dcol.dat", "r");
    fval = fopen("dval.dat", "r");
    if (!frow || !fcol || !fval) {
        fprintf(stderr,"ERROR: reopening triplet files (pass2): %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    {
        long long total_read = 0;
        while (total_read < nnz_total) {
            int n = (int)((nnz_total - total_read) < BLOCK_SIZE ? (nnz_total - total_read) : BLOCK_SIZE);

            for (int i = 0; i < n; ++i) {
                if (fscanf(frow,"%d",&drow_block[i])!=1 ||
                    fscanf(fcol,"%d",&dcol_block[i])!=1 ||
                    fscanf(fval,"%lf",&dval_block[i])!=1) {
                    fprintf(stderr,"triplet read error (pass2) at %lld\n", total_read+i);
                    exit(EXIT_FAILURE);
                }
            }

            for (int t = 0; t < num_threads; ++t) {
                targs[t] = (ThreadArgs){
                    .thread_id=t, .num_threads=num_threads,
                    .drow=drow_block, .dcol=dcol_block, .dval=dval_block,
                    .block_read=n, .ne=ne, .q=q,
                    .row_sum=row_sum,
                    .SensIn=SensIn, .SensOut=SensOut
                };
                if (pthread_create(&threads[t], NULL, worker_accumulate_sens, &targs[t]) != 0) {
                    perror("pthread_create pass2"); exit(EXIT_FAILURE);
                }
            }
            for (int t = 0; t < num_threads; ++t) pthread_join(threads[t], NULL);

            total_read += n;
        }
    }

    fclose(frow); fclose(fcol); fclose(fval);

    // cleanup
    free(drow_block); free(dcol_block); free(dval_block);
    free(row_sum); free(threads); free(targs);
}

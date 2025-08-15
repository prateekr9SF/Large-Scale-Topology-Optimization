// -----------------------------------------------------------------------------
// filter_sens_pthread_onepass.c
//
// One-pass sensitivity filtering (no volume weighting) using pthreads + atomics
// with precomputed donor row sums read from dsum.dat.
//
// Math (adjoint):
//   df/dx_i = Σ_j [ H_{j i}^q / (Σ_k H_{j k}^q) ] * (df/dx̃_j)
//
// Files (text, CWD, 1-based):
//   drow.dat : donor row j per line
//   dcol.dat : receiver col i per line
//   dval.dat : weight H_{j i} per line
//   dsum.dat : donor row sums Σ_k H_{j k}^q per line (size ne)
//
// Build:
//   gcc -O3 -march=native -ffast-math -pthread -fopenmp \
//       filter_sens_pthread_onepass.c -o filter_sens_onepass
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>

#ifdef _OPENMP
  #include <omp.h>   // for #pragma omp atomic
#endif

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 100000000    // tune for your I/O/cache
#endif

#ifndef SHARDS
#define SHARDS 8                // power-of-two (4..16 typical) to reduce atomics
#endif

#ifdef _OPENMP
  #define ATOMIC_ADD(var, val) do { _Pragma("omp atomic") (var) += (val); } while(0)
#else
  #define ATOMIC_ADD(var, val) __atomic_fetch_add(&(var), (val), __ATOMIC_RELAXED)
#endif

typedef struct {
    int thread_id, num_threads;
    int *drow, *dcol;          // buffered 1-based indices for current block
    double *w;                 // buffered weights (already ^q) for current block
    int block_read, ne;
    const double *row_sum;     // size ne: donor row sums (from dsum.dat)
    const double *SensIn;      // size ne: df/dx̃
    double **out_shard;        // SHARDS × ne accumulation buffers
} ThreadArgs;

static inline int inb(int x,int n){ return (unsigned)x < (unsigned)n; }

// Optional helper: count nnz by scanning files once (with big stdio buffers)
static long long count_nnz_from_files(void){
    FILE *fr=fopen("drow.dat","r");
    FILE *fc=fopen("dcol.dat","r");
    FILE *fv=fopen("dval.dat","r");
    if(!fr||!fc||!fv){
        fprintf(stderr,"ERROR: open drow/dcol/dval: %s\n", strerror(errno));
        if(fr) fclose(fr); if(fc) fclose(fc); if(fv) fclose(fv);
        return -1;
    }
    setvbuf(fr,NULL,_IOFBF,8<<20);
    setvbuf(fc,NULL,_IOFBF,8<<20);
    setvbuf(fv,NULL,_IOFBF,8<<20);
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

// Worker: single pass accumulate into sharded outputs using preloaded row_sum
static void *worker_accumulate_onepass(void *args_ptr)
{
    ThreadArgs *a = (ThreadArgs*)args_ptr;
    int s = (a->block_read * a->thread_id) / a->num_threads;
    int e = (a->block_read * (a->thread_id + 1)) / a->num_threads;

    for (int t = s; t < e; ++t) {
        int j = a->drow[t] - 1;             // donor j (0-based)
        int i = a->dcol[t] - 1;             // receiver i (0-based)
        if (!inb(j,a->ne) || !inb(i,a->ne)) continue;

        double denom = a->row_sum[j];       // Σ_k H_{j k}^q (precomputed)
        if (denom <= 0.0) continue;

        double contrib = (a->w[t] / denom) * a->SensIn[j];

        // shard by i’s low bits to reduce atomic contention
        int shard = i & (SHARDS - 1);
        ATOMIC_ADD(a->out_shard[shard][i], contrib);
    }
    return NULL;
}

// Public API
//  SensIn   : df/dx̃ (size ne)
//  SensOut  : df/dx  (size ne) ← single-pass result
//  ne       : number of elements
//  nnz_total: total nnz in triplet files (<=0 → auto-count)
//  q        : exponent used for weights; MUST MATCH the q used to make dsum.dat
void filterSensitivity_buffered_mts(const double *SensIn,
                                  double *SensOut,
                                  int ne,
                                  long long nnz_total)
{

    double q = 1.0;

    if (ne <= 0) { fprintf(stderr,"ERROR: ne<=0\n"); exit(EXIT_FAILURE); }

    // threads from OMP_NUM_THREADS or default 4
    int num_threads = 4;
    const char *env = getenv("OMP_NUM_THREADS");
    if (env && *env) {
        int tmp = atoi(env);
        if (tmp > 0) num_threads = tmp;
    }
    printf("filterSensitivity_onepass_mt: using %d thread(s)\n", num_threads);

    // Determine nnz if needed
    if (nnz_total <= 0) {
        nnz_total = count_nnz_from_files();
        if (nnz_total <= 0) {
            fprintf(stderr,"ERROR: could not determine nnz_total from files.\n");
            exit(EXIT_FAILURE);
        }
    }

    // Load row sums (dsum.dat) into memory
    FILE *frs = fopen("dsum.dat","r");
    if (!frs) {
        fprintf(stderr,"ERROR: open dsum.dat: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    double *row_sum = (double*)malloc((size_t)ne * sizeof(double));
    if (!row_sum) { fprintf(stderr,"alloc row_sum failed\n"); exit(EXIT_FAILURE); }
    for (int j=0; j<ne; ++j) {
        if (fscanf(frs, "%lf", &row_sum[j]) != 1) {
            fprintf(stderr,"ERROR: reading dsum.dat at line %d\n", j+1);
            exit(EXIT_FAILURE);
        }
    }
    fclose(frs);

    // Open triplets for a single streaming pass
    FILE *frow = fopen("drow.dat", "r");
    FILE *fcol = fopen("dcol.dat", "r");
    FILE *fval = fopen("dval.dat", "r");
    if (!frow || !fcol || !fval) {
        fprintf(stderr,"ERROR: opening triplet files: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    setvbuf(frow,NULL,_IOFBF,8<<20);
    setvbuf(fcol,NULL,_IOFBF,8<<20);
    setvbuf(fval,NULL,_IOFBF,8<<20);

    // Allocate block buffers & thread structs
    int    *drow_block = (int*)   malloc((size_t)BLOCK_SIZE * sizeof(int));
    int    *dcol_block = (int*)   malloc((size_t)BLOCK_SIZE * sizeof(int));
    double *dval_block = (double*)malloc((size_t)BLOCK_SIZE * sizeof(double));
    double *w_block    = (double*)malloc((size_t)BLOCK_SIZE * sizeof(double));
    if (!drow_block || !dcol_block || !dval_block || !w_block) {
        fprintf(stderr,"alloc block buffers failed\n"); exit(EXIT_FAILURE);
    }

    // Sharded outputs
    double *out_shard[SHARDS];
    for (int s=0; s<SHARDS; ++s) {
        out_shard[s] = (double*)calloc((size_t)ne, sizeof(double));
        if (!out_shard[s]) { fprintf(stderr,"alloc out_shard failed\n"); exit(EXIT_FAILURE); }
    }

    pthread_t  *threads = (pthread_t*) malloc((size_t)num_threads * sizeof(pthread_t));
    ThreadArgs *targs   = (ThreadArgs*)malloc((size_t)num_threads * sizeof(ThreadArgs));
    if (!threads || !targs) { fprintf(stderr,"alloc thread structs failed\n"); exit(EXIT_FAILURE); }

    const int use_pow = (q != 1.0);

    // Single streaming pass
    long long total_read = 0;
    while (total_read < nnz_total) {
        long long rem = nnz_total - total_read;
        int n = (int)(rem < BLOCK_SIZE ? rem : (long long)BLOCK_SIZE);

        for (int i=0; i<n; ++i) {
            if (fscanf(frow,"%d",&drow_block[i])!=1 ||
                fscanf(fcol,"%d",&dcol_block[i])!=1 ||
                fscanf(fval,"%lf",&dval_block[i])!=1) {
                fprintf(stderr,"triplet read error (one-pass) at %lld\n", total_read+i);
                exit(EXIT_FAILURE);
            }
        }

        // Precompute w = (dval)^q
        if (use_pow) {
            #pragma omp parallel for schedule(static)
            for (int i=0; i<n; ++i) w_block[i] = pow(dval_block[i], q);
        } else {
            memcpy(w_block, dval_block, (size_t)n*sizeof(double));
        }

        // Launch workers on this block
        for (int t=0; t<num_threads; ++t) {
            targs[t] = (ThreadArgs){
                .thread_id=t, .num_threads=num_threads,
                .drow=drow_block, .dcol=dcol_block, .w=w_block,
                .block_read=n, .ne=ne,
                .row_sum=row_sum, .SensIn=SensIn, .out_shard=out_shard
            };
            if (pthread_create(&threads[t], NULL, worker_accumulate_onepass, &targs[t]) != 0) {
                perror("pthread_create"); exit(EXIT_FAILURE);
            }
        }
        for (int t=0; t<num_threads; ++t) pthread_join(threads[t], NULL);

        total_read += n;
    }

    fclose(frow); fclose(fcol); fclose(fval);

    // Reduce shards → SensOut
    for (int i=0; i<ne; ++i) {
        double sum = 0.0;
        for (int s=0; s<SHARDS; ++s) sum += out_shard[s][i];
        SensOut[i] = sum;
    }

    // Cleanup
    for (int s=0; s<SHARDS; ++s) free(out_shard[s]);
    free(drow_block); free(dcol_block); free(dval_block); free(w_block);
    free(row_sum); free(threads); free(targs);
}

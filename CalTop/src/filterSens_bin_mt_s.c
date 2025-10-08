// -----------------------------------------------------------------------------
// filter_sens_pthread_onepass_bin.c
//
// One-pass sensitivity filtering (no volume weighting) using pthreads + atomics
// with precomputed donor row sums read from dsum.bin.
//
// Math (adjoint):
//   df/dx_i = Σ_j [ H_{j i}^q / (Σ_k H_{j k}^q) ] * (df/dx̃_j)
//
// Files (binary, CWD, 1-based):
//   drow.bin : donor row j (int32)
//   dcol.bin : receiver col i (int32)
//   dval.bin : weight H_{j i} (double)
//   dsum.bin : donor row sums Σ_k H_{j k}^q (double, size ne)
//
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <inttypes.h>

#include <sys/stat.h>
#include <stdint.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 100000000    // tune for I/O/cache
#endif

#ifndef SHARDS
#define SHARDS 8                // power-of-two (4..16 typical) to reduce atomics
#endif

#ifdef _OPENMP
  #define ATOMIC_ADD(var, val) do { _Pragma("omp atomic") (var) += (val); } while(0)
#else
  #define ATOMIC_ADD(var, val) __atomic_fetch_add(&(var), (val), __ATOMIC_RELAXED)
#endif

static long long file_elems(const char *path, size_t item_size) 
{
    struct stat st;
    if (stat(path, &st) != 0) { perror(path); exit(EXIT_FAILURE); }
    if (st.st_size % item_size) {
        fprintf(stderr, "Size mismatch in %s (not a multiple of item size)\n", path);
        exit(EXIT_FAILURE);
    }
    return (long long)(st.st_size / (long long)item_size);
}



typedef struct {
    int thread_id, num_threads;
    int64_t *drow, *dcol;          // buffered 1-based indices
    double *w;                 // buffered weights (already ^q)
    int block_read, ne;
    const double *row_sum;     // size ne: donor row sums (from dsum.bin)
    const double *SensIn;      // size ne: df/dx̃
    double **out_shard;        // SHARDS × ne accumulation buffers
} ThreadArgs;

static inline int inb(int x,int n){ return (unsigned)x < (unsigned)n; }

static void *worker_accumulate_onepass(void *args_ptr)
{
    ThreadArgs *a = (ThreadArgs*)args_ptr;
    int s = (a->block_read * a->thread_id) / a->num_threads;
    int e = (a->block_read * (a->thread_id + 1)) / a->num_threads;

    for (int t = s; t < e; ++t) 
    {

        int64_t j64 = a->drow[t] - 1;
        int64_t i64 = a->dcol[t] - 1;
        if (j64 < 0 || j64 >= (int64_t)a->ne || i64 < 0 || i64 >= (int64_t)a->ne) continue;


        int j = (int)j64;   // safe if ne fits in int
        int i = (int)i64;

        //int j = a->drow[t] - 1;   // donor j (0-based)
        //int i = a->dcol[t] - 1;   // receiver i
        if (!inb(j,a->ne) || !inb(i,a->ne)) continue;

        double denom = a->row_sum[j];
        if (denom <= 0.0) continue;

        double contrib = (a->w[t] / denom) * a->SensIn[j];

        int shard = i & (SHARDS - 1);
        ATOMIC_ADD(a->out_shard[shard][i], contrib);
    }
    return NULL;
}

// Public API
void filterSensitivity_bin_buffered_mts(const double *SensIn,
                                  double *SensOut,
                                  int ne,
                                  long long nnz_total)
{
    if (ne <= 0) { fprintf(stderr,"ERROR: ne<=0\n"); exit(EXIT_FAILURE); }


    int num_threads = 4;
    const char *env = getenv("OMP_NUM_THREADS");

    if (env && *env) { int tmp = atoi(env); if (tmp > 0) num_threads = tmp; }

    printf("using %d thread(s) ", num_threads);

    // Load row sums from binary
    FILE *frs = fopen("dsum.bin","rb");
    if (!frs) { perror("open dsum.bin"); exit(EXIT_FAILURE); }

    double *row_sum = (double*)malloc((size_t)ne*sizeof(double));

    if (!row_sum) { fprintf(stderr,"alloc row_sum failed\n"); exit(EXIT_FAILURE); }

    size_t got = fread(row_sum,sizeof(double),(size_t)ne,frs);
    if (got != (size_t)ne) { fprintf(stderr,"ERROR: dsum.bin short read (%zu/%d)\n",got,ne); exit(EXIT_FAILURE); }
    fclose(frs);

    // Open triplet binaries
    FILE *frow = fopen("drow.bin","rb");
    FILE *fcol = fopen("dcol.bin","rb");
    FILE *fval = fopen("dval.bin","rb");
    if (!frow||!fcol||!fval){ perror("open drow/dcol/dval.bin"); exit(EXIT_FAILURE); }
    setvbuf(frow,NULL,_IOFBF,8<<20);
    setvbuf(fcol,NULL,_IOFBF,8<<20);
    setvbuf(fval,NULL,_IOFBF,8<<20);

    // Buffers
    int64_t    *drow_block=(int64_t*)   malloc((size_t)BLOCK_SIZE*sizeof(int64_t));
    int64_t    *dcol_block=(int64_t*)   malloc((size_t)BLOCK_SIZE*sizeof(int64_t));
    double *dval_block=(double*)malloc((size_t)BLOCK_SIZE*sizeof(double));


    if(!drow_block||!dcol_block||!dval_block){fprintf(stderr,"alloc blocks failed\n");exit(EXIT_FAILURE);}
    double *w_block=(double*)malloc((size_t)BLOCK_SIZE*sizeof(double));

    // Sharded outputs
    double *out_shard[SHARDS];
    for(int s=0;s<SHARDS;++s){
        out_shard[s]=(double*)calloc((size_t)ne,sizeof(double));
        if(!out_shard[s]){fprintf(stderr,"alloc shard failed\n");exit(EXIT_FAILURE);}
    }

    pthread_t  *threads=(pthread_t*)malloc((size_t)num_threads*sizeof(pthread_t));
    ThreadArgs *targs  =(ThreadArgs*)malloc((size_t)num_threads*sizeof(ThreadArgs));

    //long long total_read=0;

    int64_t total_read = 0;

    while(1){
        size_t n1=fread(drow_block,sizeof(int64_t),   BLOCK_SIZE,frow);
        size_t n2=fread(dcol_block,sizeof(int64_t),   BLOCK_SIZE,fcol);
        size_t n3=fread(dval_block,sizeof(double),BLOCK_SIZE,fval);
        if(n1==0||n2==0||n3==0) break;
        if(n1!=n2||n1!=n3){fprintf(stderr,"ERROR: triplet block mismatch\n");exit(EXIT_FAILURE);}
        int n=(int)n1;

        memcpy(w_block,dval_block,(size_t)n*sizeof(double));

        for(int t=0;t<num_threads;++t){
            targs[t]=(ThreadArgs){t,num_threads,drow_block,dcol_block,w_block,n,ne,row_sum,SensIn,out_shard};
            if(pthread_create(&threads[t],NULL,worker_accumulate_onepass,&targs[t])!=0){perror("pthread_create");exit(EXIT_FAILURE);}
        }
        for(int t=0;t<num_threads;++t) pthread_join(threads[t],NULL);

        total_read+=n;
    }

    fclose(frow); fclose(fcol); fclose(fval);

    // Reduce shards
    for(int i=0;i<ne;++i){
        double sum=0.0;
        for(int s=0;s<SHARDS;++s) sum+=out_shard[s][i];
        SensOut[i]=sum;
    }

    for(int s=0;s<SHARDS;++s) free(out_shard[s]);
    free(drow_block); free(dcol_block); free(dval_block); free(w_block);
    free(row_sum); free(threads); free(targs);

    //printf("Processed %lld triplets (binary)\n", total_read);

    //printf("Processed %" PRId64 " triplets (binary)\n", total_read); // new
}
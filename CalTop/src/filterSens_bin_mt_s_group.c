// -----------------------------------------------------------------------------
// filter_sens_pthread_onepass_bin_multi.c
//
// One-pass sensitivity filtering (no volume weighting) using pthreads + atomics,
// with precomputed donor row sums read from dsum.bin.
//
// Math (adjoint):
//   df/dx_i = Σ_j [ H_{j i}^q / (Σ_k H_{j k}^q) ] * (df/dx̃_j)
//
// Files (binary, CWD, 1-based):
//   drow.bin : donor row j (int64)
//   dcol.bin : receiver col i (int64)
//   dval.bin : weight H_{j i}^q (double)   // already raised to q
//   dsum.bin : donor row sums Σ_k H_{j k}^q (double, size ne)
//
// Public API (multi & 3-vector wrapper):
//   void filterSensitivity_bin_buffered_mts_multi(const double* const* in,
//                                                 double* const* out,
//                                                 int nsens, int ne,
//                                                 long long nnz_total);
//   void filterSensitivity_bin_buffered_mts3(const double* dCGx,
//                                            const double* dCGy,
//                                            const double* dCGz,
//                                            double* dCGxF,
//                                            double* dCGyF,
//                                            double* dCGzF,
//                                            int ne, long long nnz_total);
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <pthread.h>

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

typedef struct {
    int thread_id, num_threads;
    int block_read, ne;

    // buffered triplets (1-based indices in files)
    const int64_t *drow;     // donor j
    const int64_t *dcol;     // receiver i
    const double  *w;        // H_{ji}^q

    // normalization and inputs
    const double *row_sum;         // size ne, Σ_k H_{jk}^q
    int nsens;                     // number of sensitivities
    const double * const *SensIn;  // [nsens][ne]

    // outputs: nsens * SHARDS shard arrays, each of size ne
    double **out_shard;
} ThreadArgs;

static inline int inb(int x, int n) { return (unsigned)x < (unsigned)n; }

static void *worker_accumulate_onepass(void *args_ptr)
{
    ThreadArgs *a = (ThreadArgs*)args_ptr;

    const int s = (a->block_read * a->thread_id) / a->num_threads;
    const int e = (a->block_read * (a->thread_id + 1)) / a->num_threads;

    for (int t = s; t < e; ++t) {

        int64_t j64 = a->drow[t] - 1; // to 0-based
        int64_t i64 = a->dcol[t] - 1; // to 0-based
        if (j64 < 0 || j64 >= (int64_t)a->ne || i64 < 0 || i64 >= (int64_t)a->ne) continue;

        const int j = (int)j64;
        const int i = (int)i64;
        if (!inb(j, a->ne) || !inb(i, a->ne)) continue;

        const double denom = a->row_sum[j];
        if (denom <= 0.0) continue;

        const double base = a->w[t] / denom;    // reuse for all sensitivities
        const int shard = i & (SHARDS - 1);

        // accumulate into each sensitivity's sharded buffer
        for (int sidx = 0; sidx < a->nsens; ++sidx) {
            const double contrib = base * a->SensIn[sidx][j];
            ATOMIC_ADD(a->out_shard[sidx*SHARDS + shard][i], contrib);
        }
    }
    return NULL;
}

static void die(const char *msg) { perror(msg); exit(EXIT_FAILURE); }

void filterSensitivity_bin_buffered_mts_multi(const double * const *SensInArr,
                                              double * const *SensOutArr,
                                              int nsens,
                                              int ne,
                                              long long nnz_total)
{
    if (ne <= 0) { fprintf(stderr,"ERROR: ne<=0\n"); exit(EXIT_FAILURE); }
    if (nsens <= 0) { fprintf(stderr,"ERROR: nsens<=0\n"); exit(EXIT_FAILURE); }

    // threads
    int num_threads = 4;
    const char *env = getenv("OMP_NUM_THREADS");
    if (env && *env) {
        int tmp = atoi(env);
        if (tmp > 0) num_threads = tmp;
    }
    printf("filterSensitivity_onepass_bin (multi): using %d thread(s)\n", num_threads);

    // load donor row sums
    FILE *frs = fopen("dsum.bin","rb");
    if (!frs) die("open dsum.bin");
    double *row_sum = (double*)malloc((size_t)ne * sizeof(double));
    if (!row_sum) { fprintf(stderr,"alloc row_sum failed\n"); exit(EXIT_FAILURE); }
    size_t got = fread(row_sum, sizeof(double), (size_t)ne, frs);
    if (got != (size_t)ne) {
        fprintf(stderr, "ERROR: dsum.bin short read (%zu/%d)\n", got, ne);
        exit(EXIT_FAILURE);
    }
    fclose(frs);

    // open triplets
    FILE *frow = fopen("drow.bin","rb");
    FILE *fcol = fopen("dcol.bin","rb");
    FILE *fval = fopen("dval.bin","rb");
    if (!frow || !fcol || !fval) die("open drow/dcol/dval.bin");

    setvbuf(frow, NULL, _IOFBF, 8u<<20);
    setvbuf(fcol, NULL, _IOFBF, 8u<<20);
    setvbuf(fval, NULL, _IOFBF, 8u<<20);

    // buffers
    int64_t *drow_block = (int64_t*)malloc((size_t)BLOCK_SIZE * sizeof(int64_t));
    int64_t *dcol_block = (int64_t*)malloc((size_t)BLOCK_SIZE * sizeof(int64_t));
    double  *dval_block = (double*)  malloc((size_t)BLOCK_SIZE * sizeof(double));
    double  *w_block    = (double*)  malloc((size_t)BLOCK_SIZE * sizeof(double));
    if (!drow_block || !dcol_block || !dval_block || !w_block) {
        fprintf(stderr, "alloc blocks failed\n"); exit(EXIT_FAILURE);
    }

    // sharded outputs for each sensitivity
    double **out_shard = (double**)malloc((size_t)(nsens * SHARDS) * sizeof(double*));
    if (!out_shard) { fprintf(stderr, "alloc out_shard ptrs failed\n"); exit(EXIT_FAILURE); }
    for (int s = 0; s < nsens * SHARDS; ++s) {
        out_shard[s] = (double*)calloc((size_t)ne, sizeof(double));
        if (!out_shard[s]) { fprintf(stderr, "alloc shard failed\n"); exit(EXIT_FAILURE); }
    }

    // threads
    pthread_t  *threads = (pthread_t*)malloc((size_t)num_threads * sizeof(pthread_t));
    ThreadArgs *targs   = (ThreadArgs*) malloc((size_t)num_threads * sizeof(ThreadArgs));
    if (!threads || !targs) { fprintf(stderr,"alloc thread objects failed\n"); exit(EXIT_FAILURE); }

    long long total_read = 0;

    while (1) {
        const size_t n1 = fread(drow_block, sizeof(int64_t), BLOCK_SIZE, frow);
        const size_t n2 = fread(dcol_block, sizeof(int64_t), BLOCK_SIZE, fcol);
        const size_t n3 = fread(dval_block, sizeof(double),  BLOCK_SIZE, fval);
        if (n1 == 0 || n2 == 0 || n3 == 0) break;
        if (n1 != n2 || n1 != n3) {
            fprintf(stderr, "ERROR: triplet block mismatch\n"); exit(EXIT_FAILURE);
        }

        const int n = (int)n1;
        memcpy(w_block, dval_block, (size_t)n * sizeof(double));

        for (int t = 0; t < num_threads; ++t) {
            targs[t] = (ThreadArgs){
                .thread_id = t, .num_threads = num_threads,
                .block_read = n, .ne = ne,
                .drow = drow_block, .dcol = dcol_block, .w = w_block,
                .row_sum = row_sum,
                .nsens = nsens,
                .SensIn = SensInArr,
                .out_shard = out_shard
            };
            if (pthread_create(&threads[t], NULL, worker_accumulate_onepass, &targs[t]) != 0) {
                die("pthread_create");
            }
        }
        for (int t = 0; t < num_threads; ++t) pthread_join(threads[t], NULL);

        total_read += n;
    }

    fclose(frow); fclose(fcol); fclose(fval);

    // reduce shards per sensitivity
    for (int sidx = 0; sidx < nsens; ++sidx) {
        double *out = SensOutArr[sidx];
        for (int i = 0; i < ne; ++i) {
            double sum = 0.0;
            for (int sh = 0; sh < SHARDS; ++sh) sum += out_shard[sidx*SHARDS + sh][i];
            out[i] = sum;
        }
    }

    // cleanup
    for (int s = 0; s < nsens * SHARDS; ++s) free(out_shard[s]);
    free(out_shard);
    free(drow_block); free(dcol_block); free(dval_block); free(w_block);
    free(row_sum); free(threads); free(targs);

    printf("Processed %lld triplets (binary)\n", total_read);
}

// Convenience wrapper for exactly 3 sensitivities (CGx, CGy, CGz)
void filterSensitivity_bin_buffered_mts3(const double *dCGx,
                                         const double *dCGy,
                                         const double *dCGz,
                                         double *dCGxFiltered,
                                         double *dCGyFiltered,
                                         double *dCGzFiltered,
                                         int ne, long long nnz_total)
{
    const double *in[3]  = { dCGx, dCGy, dCGz };
    double *out[3]       = { dCGxFiltered, dCGyFiltered, dCGzFiltered };
    filterSensitivity_bin_buffered_mts_multi(in, out, 3, ne, nnz_total);
}

/*
Usage:

// Before: three separate passes
// filterSensitivity_bin_buffered_mts(dCGx, dCGxF, ne, nnz);
// filterSensitivity_bin_buffered_mts(dCGy, dCGyF, ne, nnz);
// filterSensitivity_bin_buffered_mts(dCGz, dCGzF, ne, nnz);

// Now: one pass
filterSensitivity_bin_buffered_mts3(
    dCGx, dCGy, dCGz,
    dCGxFiltered, dCGyFiltered, dCGzFiltered,
    ne, filternnz
);
*/

// pthread_filtering.c  (minimal changes for binary I/O + dsum.bin)
#include "filterDensity_mt.h"
#include <pthread.h>
#include <math.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>   

pthread_mutex_t *row_locks;

#include <sys/stat.h>
#include <stdint.h> 

#ifdef PROFILING_ON
#include <TAU.h>
#endif

/*
static long long file_elems(const char *path, size_t item_size) {
    struct stat st;
    if (stat(path, &st) != 0) { perror(path); exit(EXIT_FAILURE); }
    if (st.st_size % item_size) {
        fprintf(stderr, "Size mismatch in %s (not a multiple of item size)\n", path);
        exit(EXIT_FAILURE);
    }
    return (long long)(st.st_size / (long long)item_size);
}

*/

static int64_t file_elems64(const char *path, size_t item_size) 
{
    struct stat st;
    if (stat(path, &st) != 0)
    { 
        perror(path); exit(EXIT_FAILURE); 
    }
    if ((size_t)st.st_size % item_size) 
    {
        fprintf(stderr, "Size mismatch in %s (not a multiple of item size)\n", path);
        exit(EXIT_FAILURE);
    }
    return (int64_t)((size_t)st.st_size / item_size);
}


void *thread_filter_bin_worker_mutex(void *args_ptr)
{
    ThreadArgs *args = (ThreadArgs *)args_ptr;

    int start = (args->block_read * args->thread_id) / args->num_threads;
    int end   = (args->block_read * (args->thread_id + 1)) / args->num_threads;

    for (int i = start; i < end; ++i)
    {
        int64_t row64 = args->drow[i] - 1;
        int64_t col64 = args->dcol[i] - 1;

        if (row64 < 0 || row64 >= (int64_t)args->ne || col64 < 0 || col64 >= (int64_t)args->ne) 
        {
            continue;
        }

        int row = (int)row64;   // safe if args->ne fits in int
        int col = (int)col64;

        //double w = pow(args->dval[i], args->q);
        double w = args->dval[i];
        double contrib = w * args->Vector[col];

        pthread_mutex_lock(&row_locks[row]);
        args->VectorFiltered[row] += contrib;
        // NOTE: we leave this line intact, but we'll normalize with dsum.bin later
        if (args->weight_sum) args->weight_sum[row] += w;
        pthread_mutex_unlock(&row_locks[row]);
    }

    return NULL;
}

void filterDensity_buffered_bin_mt(double *Vector, double *VectorFiltered,
                               int *filternnzElems,
                               int *ne_ptr, int *fnnzassumed_ptr,
                               double *q_ptr, int filternnz_total)
{
    
    #ifdef PROFILING_ON
        // Timer for entire density filter function
        TAU_PROFILE_TIMER(t_filter_total,"Density Filter: Total","",TAU_USER);

        // Timer for Disk I/O
        TAU_PROFILE_TIMER(t_filter_io,"Density Filter: Binary I/O","",TAU_USER);

        // Timer for parallel filter application
        TAU_PROFILE_TIMER(t_filter_apply,"Density Filter: Parallel Apply","",TAU_USER);
    #endif

    // Profile the entire density filter 
    #ifdef PROFILING_ON
        TAU_PROFILE_START(t_filter_total);
    #endif
    
    const int ne = *ne_ptr;
    const double q = *q_ptr;

    int64_t nrow = file_elems64("drow.bin", sizeof(int64_t));
    int64_t ncol = file_elems64("dcol.bin", sizeof(int64_t));
    int64_t nval = file_elems64("dval.bin", sizeof(double));

    //filternnz_total = (int)nrow;   // use this instead of a passed-in value

    int64_t filternnz_total_64 = nrow;             // total triplets (symmetric already)

    printf("Filternnz total: %" PRId64 "\n", filternnz_total_64);

    if (!(nrow == ncol && ncol == nval)) 
    {
        fprintf(stderr, "Triplet file length mismatch: drow=%" PRId64 " dcol=%" PRId64 " dval=%" PRId64 "\n",nrow, ncol, nval);
        exit(EXIT_FAILURE);
    }

    const char *env_threads = getenv("OMP_NUM_THREADS");
    int num_threads = (env_threads ? atoi(env_threads) : 4);
    if (num_threads <= 0) num_threads = 4;

    printf("Using %d thread(s) to filter vector\n", num_threads);

    // Time the binary file I/O
    #ifdef PROFILING_ON
        TAU_PROFILE_START(t_filter_io);
    #endif

    // --- Open BINARY triplet files ---
    FILE *frow = fopen("drow.bin", "rb");
    FILE *fcol = fopen("dcol.bin", "rb");
    FILE *fval = fopen("dval.bin", "rb");

    #ifdef PROFILING_ON
        TAU_PROFILE_STOP(t_filter_io);
    #endif

    if (!frow || !fcol || !fval) 
    {
        perror("Error opening drow.bin/dcol.bin/dval.bin");
        exit(EXIT_FAILURE);
    }

    // Large stdio buffers help when files are big
    setvbuf(frow, NULL, _IOFBF, 8<<20);
    setvbuf(fcol, NULL, _IOFBF, 8<<20);
    setvbuf(fval, NULL, _IOFBF, 8<<20);

    // Allocate memory for dsum here to avoid in I/O profiling
    double *dsum = (double*)malloc((size_t)ne * sizeof(double));
    if (!dsum) 
    {
        fprintf(stderr, "Failed to allocate dsum\n");
        exit(EXIT_FAILURE);
    }

    // Time dsum I/O
    #ifdef PROFILING_ON
        TAU_PROFILE_START(t_filter_io);
    #endif

    // --- Read dsum.bin once (row-wise denominators, length = ne) ---
    FILE *fdsum = fopen("dsum.bin", "rb");
    if (!fdsum) 
    {
        perror("Error opening dsum.bin");
        exit(EXIT_FAILURE);
    }

    size_t got = fread(dsum, sizeof(double), (size_t)ne, fdsum);
    fclose(fdsum);

    #ifdef PROFILING_ON
        TAU_PROFILE_STOP(t_filter_io);
    #endif

    if (got != (size_t)ne) 
    {
        fprintf(stderr, "Short read from dsum.bin (got %zu of %d)\n", got, ne);
        free(dsum);
        exit(EXIT_FAILURE);
    }

    // --- Allocate block buffers (triplets) ---
    int64_t    *drow_block = (int64_t*)   malloc((size_t)BLOCK_SIZE * sizeof(int64_t));
    int64_t    *dcol_block = (int64_t*)   malloc((size_t)BLOCK_SIZE * sizeof(int64_t));
    double *dval_block = (double*)malloc((size_t)BLOCK_SIZE * sizeof(double));
    
    if (!drow_block || !dcol_block || !dval_block) 
    {
        fprintf(stderr, "Block buffer allocation failed.\n");
        free(dsum);
        exit(EXIT_FAILURE);
    }

    // Output + (legacy) weight_sum buffer (left in place; normalization uses dsum)
    for (int i = 0; i < ne; ++i) VectorFiltered[i] = 0.0;
    
    double *weight_sum = (double*)calloc((size_t)ne, sizeof(double)); // optional; kept to avoid touching worker
    
    if (!weight_sum)
    {
        fprintf(stderr, "Allocation failed for weight_sum\n");
        free(dsum); free(drow_block); free(dcol_block); free(dval_block);
        exit(EXIT_FAILURE);
    }

    // Threads & locks
    pthread_t *threads = (pthread_t*)malloc((size_t)num_threads * sizeof(pthread_t));
    ThreadArgs *thread_args = (ThreadArgs*)malloc((size_t)num_threads * sizeof(ThreadArgs));

    if (!threads || !thread_args) 
    {
        fprintf(stderr, "Thread arrays allocation failed\n");
        free(dsum); free(drow_block); free(dcol_block); free(dval_block); free(weight_sum);
        exit(EXIT_FAILURE);
    }

    row_locks = (pthread_mutex_t*)malloc((size_t)ne * sizeof(pthread_mutex_t));
    for (int i = 0; i < ne; ++i) pthread_mutex_init(&row_locks[i], NULL);

    // --- Stream blocks in binary ---
    int64_t total_read = 0;
    while (total_read < filternnz_total_64) 
    {
        int64_t remaining = filternnz_total_64 - total_read;
        int block_read = (remaining > BLOCK_SIZE) ? BLOCK_SIZE : (int)remaining;

        // Time block-read of binary files
        #ifdef PROFILING_ON
            TAU_PROFILE_START(t_filter_io);
        #endif

        //if (block_read > BLOCK_SIZE) block_read = BLOCK_SIZE;
        size_t r1 = fread(drow_block, sizeof(int64_t),    (size_t)block_read, frow);
        size_t r2 = fread(dcol_block, sizeof(int64_t),    (size_t)block_read, fcol);
        size_t r3 = fread(dval_block, sizeof(double), (size_t)block_read, fval);

        #ifdef PROFILING_ON
            TAU_PROFILE_STOP(t_filter_io);
        #endif

        if (r1 != (size_t)block_read || r2 != (size_t)block_read || r3 != (size_t)block_read) 
        {
            fprintf(stderr, "Binary read error at triplet offset %" PRId64 " (got %zu/%zu/%zu)\n", total_read, r1, r2, r3);
            fclose(frow);
            fclose(fcol); 
            fclose(fval);
            free(dsum); free(drow_block); free(dcol_block); free(dval_block);
            free(weight_sum);
            
            for (int i = 0; i < ne; ++i) pthread_mutex_destroy(&row_locks[i]);
            free(row_locks); free(threads); free(thread_args);
            exit(EXIT_FAILURE);
        }

        // Profile parallel filter application
        #ifdef PROFILING_ON
            TAU_PROFILE_START(t_filter_apply);
        #endif

        // Launch threads on this block
        for (int t = 0; t < num_threads; ++t) 
        {
            thread_args[t] = (ThreadArgs)
            {
                .thread_id       = t,
                .num_threads     = num_threads,
                .drow            = drow_block,
                .dcol            = dcol_block,
                .dval            = dval_block,
                .Vector          = Vector,
                .VectorFiltered  = VectorFiltered,
                .weight_sum      = weight_sum,   // kept; final norm uses dsum
                .q               = q,
                .ne              = ne,
                .block_read      = block_read
            };
            
            if (pthread_create(&threads[t], NULL, thread_filter_bin_worker_mutex, &thread_args[t]) != 0) 
            {
                perror("Failed to create thread");
                exit(EXIT_FAILURE);
            }
        }
        for (int t = 0; t < num_threads; ++t) 
        {
            if (pthread_join(threads[t], NULL) != 0) 
            {
                perror("Failed to join thread");
                exit(EXIT_FAILURE);
            }
        }

        #ifdef PROFILING_ON
            TAU_PROFILE_STOP(t_filter_apply);
        #endif

        total_read += block_read;
    }

    #ifdef PROFILING_ON
        TAU_PROFILE_START(t_filter_io);
    #endif

    fclose(frow);
    fclose(fcol);
    fclose(fval);

    #ifdef PROFILING_ON
        TAU_PROFILE_STOP(t_filter_io);
    #endif

    // --- Normalize using precomputed dsum.bin ---
    for (int i = 0; i < ne; ++i) 
    {
        VectorFiltered[i] = (dsum[i] > 0.0) ? (VectorFiltered[i] / dsum[i]) : 0.0;
    }

    // Cleanup
    free(dsum);
    free(drow_block); free(dcol_block); free(dval_block);
    free(weight_sum);
    for (int i = 0; i < ne; ++i) pthread_mutex_destroy(&row_locks[i]);
    free(row_locks);
    free(threads);
    free(thread_args);

    #ifdef PROFILING_ON
    TAU_PROFILE_STOP(t_filter_total);
    #endif
}
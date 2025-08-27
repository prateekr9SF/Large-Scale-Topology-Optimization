// pthread_filtering.c
#include "filterDensity_mt.h"
#include <pthread.h>
#include <math.h>
#include <stdatomic.h>


pthread_mutex_t *row_locks;


void *thread_filter_worker_mutex(void *args_ptr) 
{

    /**
    * @brief Thread worker that applies one block of filter triplets to
    *        accumulate weighted contributions into filtered element densities.
    *
 * Each thread is assigned a contiguous block of triplets [start, end)
 * within the current buffer. For each triplet i:
 *
 *    row = drow[i] - 1;    // receiver element (0-based)
 *    col = dcol[i] - 1;    // donor element (0-based)
 *    w   = dval[i];        // weight, precomputed as (rmin - dist)
 *
 * The contribution from donor "col" is accumulated as:
 *
 *    VectorFiltered[row] += w * Vector[col];
 *    weight_sum[row]     += w;
 *
 * Updates to per-row accumulators are protected by a mutex
 * (row_locks[row]) to avoid race conditions across threads.
 *
 * @param args_ptr Pointer to a ThreadArgs struct with the following fields:
 *        - int thread_id        : Thread index (0 .. num_threads-1).
 *        - int num_threads      : Total number of threads.
 *        - int block_read       : Number of triplets in the current block.
 *        - int *drow, *dcol     : Arrays of row and col indices (1-based).
 *        - double *dval         : Array of triplet weights.
 *        - int ne               : Number of elements (vector length).
 *        - const double *Vector : Input vector of element densities.
 *        - double *VectorFiltered : Output vector (accumulated values).
 *        - double *weight_sum   : Row-wise sum of weights (for normalization).
 *
 * @return Always returns NULL (pthread signature).
 *
 * @note
 * - Indices in drow/dcol are 1-based; they are converted to 0-based here.
 * - Rows with no incident triplets will retain weight_sum[row] == 0,
 *   so normalization must guard against divide-by-zero.
 * - Duplicate triplets are allowed: they scale both numerator and
 *   denominator equally, so normalized results are unaffected.
 */

    ThreadArgs *args = (ThreadArgs *)args_ptr;

    int start = (args->block_read * args->thread_id) / args->num_threads;
    int end   = (args->block_read * (args->thread_id + 1)) / args->num_threads;

    printf("Thread %d started: processing from index %d to %d\n", args->thread_id, start, end);

    for (int i = start; i < end; ++i) 
    {
        int row = args->drow[i] - 1;
        int col = args->dcol[i] - 1;

        if (row < 0 || row >= args->ne || col < 0 || col >= args->ne) 
        {
            fprintf(stderr, "Invalid row/col index: row=%d, col=%d\n", row, col);
            continue;
        }

        //double w = pow(args->dval[i], args->q);
        double w = args-> dval[i];
        double contrib = w * args->Vector[col];

        pthread_mutex_lock(&row_locks[row]);
        args->VectorFiltered[row] += contrib;
        args->weight_sum[row] += w;
        pthread_mutex_unlock(&row_locks[row]);
    }

    return NULL;
}


void filterDensity_buffered_dat_mt(double *Vector, double *VectorFiltered,
                              int *filternnzElems,
                              int *ne_ptr, int *fnnzassumed_ptr,
                              double *q_ptr, int filternnz_total)
{

    /**
 * @brief Apply the density filter to a vector using buffered, multi-threaded processing.
 *
 * This routine reads the sparse filter matrix triplets (row, col, value) from
 * disk files written previously (`drow.dat`, `dcol.dat`, `dval.dat`), partitions
 * them into blocks of size BLOCK_SIZE, and processes each block in parallel using
 * pthreads. Each triplet contributes:
 *
 *    VectorFiltered[row] += w * Vector[col];
 *    weight_sum[row]     += w;
 *
 * where w = dval[i] (the precomputed linear weight rmin − dist).
 *
 * After all blocks are processed, the result is normalized row-wise:
 *
 *    VectorFiltered[i] = (weight_sum[i] > 0.0) ?
 *                        VectorFiltered[i] / weight_sum[i] : 0.0;
 *
 * @param Vector            [in]  Original per-element values (length = ne).
 * @param VectorFiltered    [out] Filtered values (length = ne, zeroed here then filled).
 * @param filternnzElems    [in]  Array of neighbor counts per element (not used here).
 * @param ne_ptr            [in]  Pointer to number of elements.
 * @param fnnzassumed_ptr   [in]  Pointer to assumed maximum nnz per row (not used here).
 * @param q_ptr             [in]  Pointer to exponent q (ignored in current implementation).
 * @param filternnz_total   [in]  Total number of nonzero triplets to read and process.
 *
 * @details
 * - Triplet files (`drow.dat`, `dcol.dat`, `dval.dat`) must contain exactly
 *   filternnz_total entries each, with 1-based indices.
 * - Duplicate triplets are allowed; they scale numerator and denominator equally,
 *   so normalized results are unaffected.
 * - Rows with no incident weights yield VectorFiltered[i] = 0.0.
 *
 * @sideeffects
 * - Allocates temporary buffers and per-row mutexes.
 * - Opens and reads drow.dat, dcol.dat, dval.dat from disk.
 * - Prints progress messages to stdout.
 *
 * @note
 * - Thread safety across rows is enforced with one pthread_mutex_t per row.
 * - If BLOCK_SIZE is large, memory usage for buffers may be significant.
 *
 * @return void
 */


    int ne = *ne_ptr;
    int fnnzassumed = *fnnzassumed_ptr;
    double q = *q_ptr;

    printf("Number of elements:%d \n", ne );
    size_t bytes = (size_t)ne * sizeof(double);
    printf("Attempting to allocate %.2f MB for weight_sum\n", bytes / (1024.0 * 1024.0));

    const char *env_threads = getenv("OMP_NUM_THREADS");

    int num_threads = (env_threads!= NULL) ? atoi(env_threads): 4;

    if (num_threads <= 0)
    {
        fprintf(stderr,"Invalid OMP_NUM_THREADS setting; falling back to 4 threads. \n");
        num_threads = 4;
    }

    printf("Using %d thread(s) to filter vector \n", num_threads);

    printf("Opening filter matrix files...");
    FILE *frow = fopen("drow.dat", "r");
    FILE *fcol = fopen("dcol.dat", "r");
    FILE *fval = fopen("dval.dat", "r");
    printf("Done!\n");

    if (!frow || !fcol || !fval) {
        perror("Error opening filter input files");
        exit(EXIT_FAILURE);
    }

    int *drow_block = malloc(BLOCK_SIZE * sizeof(int));
    int *dcol_block = malloc(BLOCK_SIZE * sizeof(int));
    double *dval_block = malloc(BLOCK_SIZE * sizeof(double));
    

    if (!drow_block || !dcol_block || !dval_block) 
    {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    printf("Allocating memory for weights...");
    double *weight_sum = calloc(ne, sizeof(double));
    if (!weight_sum) 
    {
        fprintf(stderr, "Memory allocation for weights failed.\n");
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");

    // Zero-initialize VectorFiltered (caller owns allocation)
    for (int i = 0; i < ne; ++i) VectorFiltered[i] = 0.0;

    

    // Allocate arrays for pthread management
    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    ThreadArgs *thread_args = malloc(num_threads * sizeof(ThreadArgs));

    row_locks = malloc(ne * sizeof(pthread_mutex_t));
    for (int i = 0; i < ne; ++i) pthread_mutex_init(&row_locks[i], NULL);


    int total_read = 0;
    while (total_read < filternnz_total) 
    {
        int block_read = (filternnz_total - total_read < BLOCK_SIZE) ? (filternnz_total - total_read) : BLOCK_SIZE;

        for (int i = 0; i < block_read; ++i) {
            if (fscanf(frow, "%d", &drow_block[i]) != 1 ||
                fscanf(fcol, "%d", &dcol_block[i]) != 1 ||
                fscanf(fval, "%lf", &dval_block[i]) != 1) {
                fprintf(stderr, "Error reading triplet %d from disk\n", i + total_read);
                exit(EXIT_FAILURE);
            }
        }
        printf("Launching threads to process the block...\n");
        for (int t = 0; t < num_threads; ++t) 
        {
            thread_args[t] = (ThreadArgs){
                .thread_id = t,
                .num_threads = num_threads,
                .drow = drow_block,
                .dcol = dcol_block,
                .dval = dval_block,
                .Vector = Vector,
                .VectorFiltered = VectorFiltered,
                .weight_sum = weight_sum,
                .ne = ne,
                .block_read = block_read
            };

            /*
            if (pthread_create(&threads[t], NULL, thread_filter_worker_atomic, &thread_args[t]) != 0) {
                perror("Failed to create thread");
                exit(EXIT_FAILURE);
            }
            */

            
            if (pthread_create(&threads[t], NULL, thread_filter_worker_mutex, &thread_args[t]) != 0) 
            {
                perror("Failed to create thread");
                exit(EXIT_FAILURE);
            }

            

        }

        printf("Waiting on all threads to finish...\n");

        for (int t = 0; t < num_threads; ++t) 
        {
            if (pthread_join(threads[t], NULL) != 0) 
            {
                perror("Failed to join thread");
                exit(EXIT_FAILURE);
            }
        }

        total_read += block_read;
    }

    for (int i = 0; i < ne; ++i) {
        VectorFiltered[i] = (weight_sum[i] > 0.0) ? VectorFiltered[i] / weight_sum[i] : 0.0;
    }

    fclose(frow); fclose(fcol); fclose(fval);
    free(drow_block); free(dcol_block); free(dval_block); free(weight_sum);

    for (int i = 0; i < ne; ++i) pthread_mutex_destroy(&row_locks[i]);
    free(row_locks);
}


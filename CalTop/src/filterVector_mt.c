// pthread_filtering.c
#include "filterVector_mt.h"
#include <pthread.h>
#include <math.h>
#include <stdatomic.h>

void *thread_filter_worker_atomic(void *args_ptr) 
{
    ThreadArgs *args = (ThreadArgs *)args_ptr;

    int start = (args->block_read * args->thread_id) / args->num_threads;
    int end   = (args->block_read * (args->thread_id + 1)) / args->num_threads;

    


    for (int i = start; i < end; ++i) {
        int row = args->drow[i] - 1;
        int col = args->dcol[i] - 1;
        double w = pow(args->dval[i], args->q);
        double contrib = w * args->Vector[col];

        // Atomic operations using OpenMP atomic for portability
        #pragma omp atomic
        args->VectorFiltered[row] += contrib;

        #pragma omp atomic
        args->weight_sum[row] += w;
    }

    return NULL;
}



pthread_mutex_t *row_locks;


void *thread_filter_worker_mutex(void *args_ptr) 
{
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

        double w = pow(args->dval[i], args->q);
        double contrib = w * args->Vector[col];

        pthread_mutex_lock(&row_locks[row]);
        args->VectorFiltered[row] += contrib;
        args->weight_sum[row] += w;
        pthread_mutex_unlock(&row_locks[row]);
    }

    return NULL;
}

// tanh-projection thread worker

void *thread_filter_worker_projected(void *args_ptr) 
{
    ThreadArgs *args = (ThreadArgs *)args_ptr;

    int start = (args->block_read * args->thread_id) / args->num_threads;
    int end   = (args->block_read * (args->thread_id + 1)) / args->num_threads;

    for (int i = start; i < end; ++i) 
    {
        int row = args->drow[i] - 1;
        int col = args->dcol[i] - 1;

        if (row < 0 || row >= args->ne || col < 0 || col >= args->ne) continue;

        double w = pow(args->dval[i], args->q);
        double contrib = w * args->Vector[col];

        pthread_mutex_lock(&row_locks[row]);
        args->VectorFiltered[row] += contrib;
        args->weight_sum[row] += w;
        pthread_mutex_unlock(&row_locks[row]);
    }

    return NULL;
}



void filterVector_buffered_mt(double *Vector, double *VectorFiltered,
                              int *filternnzElems,
                              int *ne_ptr, int *fnnzassumed_ptr,
                              double *q_ptr, int filternnz_total)
{

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
                .q = q,
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


void filterVector_projected_mt(double *rho, double *rho_projected,
                               int *ne_ptr, int *fnnzassumed_ptr,
                               double *q_ptr, int filternnz_total,
                               double beta, double eta)
{
    int ne = *ne_ptr;
    int fnnzassumed = *fnnzassumed_ptr;
    double q = *q_ptr;

    const char *env_threads = getenv("OMP_NUM_THREADS");
    int num_threads = (env_threads != NULL) ? atoi(env_threads) : 4;
    if (num_threads <= 0) num_threads = 4;

    printf("Using %d thread(s) for projected density filtering...\n", num_threads);

    FILE *frow = fopen("drow.dat", "r");
    FILE *fcol = fopen("dcol.dat", "r");
    FILE *fval = fopen("dval.dat", "r");
    if (!frow || !fcol || !fval) {
        perror("Error opening filter input files");
        exit(EXIT_FAILURE);
    }

    int *drow_block = malloc(BLOCK_SIZE * sizeof(int));
    int *dcol_block = malloc(BLOCK_SIZE * sizeof(int));
    double *dval_block = malloc(BLOCK_SIZE * sizeof(double));

    double *weight_sum = calloc(ne, sizeof(double));
    double *filtered = calloc(ne, sizeof(double));
    if (!drow_block || !dcol_block || !dval_block || !weight_sum || !filtered) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < ne; ++i) rho_projected[i] = 0.0;

    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    ThreadArgs *thread_args = malloc(num_threads * sizeof(ThreadArgs));
    row_locks = malloc(ne * sizeof(pthread_mutex_t));
    for (int i = 0; i < ne; ++i) pthread_mutex_init(&row_locks[i], NULL);

    int total_read = 0;
    while (total_read < filternnz_total) {
        int block_read = (filternnz_total - total_read < BLOCK_SIZE) ? (filternnz_total - total_read) : BLOCK_SIZE;

        for (int i = 0; i < block_read; ++i) {
            if (fscanf(frow, "%d", &drow_block[i]) != 1 ||
                fscanf(fcol, "%d", &dcol_block[i]) != 1 ||
                fscanf(fval, "%lf", &dval_block[i]) != 1) {
                fprintf(stderr, "Error reading triplet at line %d\n", total_read + i);
                exit(EXIT_FAILURE);
            }
        }

        for (int t = 0; t < num_threads; ++t) {
            thread_args[t] = (ThreadArgs){
                .thread_id = t,
                .num_threads = num_threads,
                .drow = drow_block,
                .dcol = dcol_block,
                .dval = dval_block,
                .Vector = rho,
                .VectorFiltered = filtered,
                .weight_sum = weight_sum,
                .q = q,
                .ne = ne,
                .block_read = block_read
            };

            if (pthread_create(&threads[t], NULL, thread_filter_worker_projected, &thread_args[t]) != 0) {
                perror("Thread creation failed");
                exit(EXIT_FAILURE);
            }
        }

        for (int t = 0; t < num_threads; ++t) {
            pthread_join(threads[t], NULL);
        }

        total_read += block_read;
    }

    // Apply tanh projection
    double tanh_beta_eta = tanh(beta * eta);
    double tanh_beta_1_eta = tanh(beta * (1.0 - eta));
    double denom = tanh_beta_eta + tanh_beta_1_eta;

    for (int i = 0; i < ne; ++i) {
        double rho_bar = (weight_sum[i] > 0.0) ? (filtered[i] / weight_sum[i]) : 0.0;
        rho_projected[i] = (tanh_beta_eta + tanh(beta * (rho_bar - eta))) / denom;
    }

    fclose(frow); fclose(fcol); fclose(fval);
    free(drow_block); free(dcol_block); free(dval_block);
    free(weight_sum); free(filtered);
    for (int i = 0; i < ne; ++i) pthread_mutex_destroy(&row_locks[i]);
    free(row_locks); free(threads); free(thread_args);
}

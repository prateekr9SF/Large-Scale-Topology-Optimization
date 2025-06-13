// pthread_filtering.c
#include "filterVector_mt.h"
#include <pthread.h>
#include <math.h>
#include <stdatomic.h>

void *thread_filter_worker_atomic(void *args_ptr) {
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

    printf("Opening files...\n");
    FILE *frow = fopen("drow.dat", "r");
    FILE *fcol = fopen("dcol.dat", "r");
    FILE *fval = fopen("dval.dat", "r");
    printf("done opening files \n");
    if (!frow || !fcol || !fval) {
        perror("Error opening filter input files");
        exit(EXIT_FAILURE);
    }

    int *drow_block = malloc(BLOCK_SIZE * sizeof(int));
    int *dcol_block = malloc(BLOCK_SIZE * sizeof(int));
    double *dval_block = malloc(BLOCK_SIZE * sizeof(double));
    

    if (!drow_block || !dcol_block || !dval_block) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }
    printf("Allocating memeory for weights...\n");
    double *weight_sum = calloc(ne, sizeof(double));
    if (!weight_sum) 
    {
        fprintf(stderr, "Memory allocation for weights failed.\n");
        exit(EXIT_FAILURE);
    }
    printf("Done!");

    // Zero-initialize VectorFiltered (caller owns allocation)
    for (int i = 0; i < ne; ++i) VectorFiltered[i] = 0.0;

    int num_threads = 8;
    pthread_t threads[num_threads];
    ThreadArgs thread_args[num_threads];

    int total_read = 0;
    while (total_read < filternnz_total) {
        int block_read = (filternnz_total - total_read < BLOCK_SIZE) ? (filternnz_total - total_read) : BLOCK_SIZE;

        for (int i = 0; i < block_read; ++i) {
            if (fscanf(frow, "%d", &drow_block[i]) != 1 ||
                fscanf(fcol, "%d", &dcol_block[i]) != 1 ||
                fscanf(fval, "%lf", &dval_block[i]) != 1) {
                fprintf(stderr, "Error reading triplet %d from disk\n", i + total_read);
                exit(EXIT_FAILURE);
            }
        }

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
            if (pthread_create(&threads[t], NULL, thread_filter_worker_atomic, &thread_args[t]) != 0) {
                perror("Failed to create thread");
                exit(EXIT_FAILURE);
            }
        }

        for (int t = 0; t < num_threads; ++t) {
            if (pthread_join(threads[t], NULL) != 0) {
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
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "CalculiX.h"

#define MAX_NNZ_PER_THREAD 1000  // maximum buffer space for all threads combined

// Struct to pass arguments to each pthread worker
typedef struct {
    int thread_id;
    int ne0;                  // total number of elements
    int ne_start, ne_end;     // element range for this thread
    double *elCentroid;       // element centroids
    double rmin_local;        // filter radius
    int *fnnzassumed;         // max nnz allowed per row (safety)
    int *filternnzElems;      // output: number of neighbors per element

    int *drow;                // thread-local buffer for row indices
    int *dcol;                // thread-local buffer for col indices
    double *dval;             // thread-local buffer for weights
    int *nnz_count;           // nnz count per element for this thread
} ThreadArgs;

// Worker function: computes neighbors and writes filter triplets into thread-local buffers
void *filter_thread_worker(void *args_ptr) {
    ThreadArgs *args = (ThreadArgs *)args_ptr;
    int count = 0;  // offset within thread-local buffers

    for (int i = args->ne_start; i < args->ne_end; ++i) {
        int local_nnz = 0;

        // Get centroid of element i
        double xi = args->elCentroid[3 * i + 0];
        double yi = args->elCentroid[3 * i + 1];
        double zi = args->elCentroid[3 * i + 2];

        // Loop over all other elements to compute neighbors within rmin
        for (int j = 0; j < args->ne0; ++j) {
            if (i == j) continue; // skip self

            double xj = args->elCentroid[3 * j + 0];
            double yj = args->elCentroid[3 * j + 1];
            double zj = args->elCentroid[3 * j + 2];

            double dx = xi - xj;
            double dy = yi - yj;
            double dz = zi - zj;
            double dist = sqrt(dx * dx + dy * dy + dz * dz);

            if (dist <= args->rmin_local) {
                double w = args->rmin_local - dist;

                // Write (i+1, j+1, w)
                int idx = count + local_nnz * 2;
                args->drow[idx]     = i + 1;
                args->dcol[idx]     = j + 1;
                args->dval[idx]     = w;

                // Write symmetric entry (j+1, i+1, w)
                args->drow[idx + 1] = j + 1;
                args->dcol[idx + 1] = i + 1;
                args->dval[idx + 1] = w;

                local_nnz++;
            }
        }

        // Save neighbor count and row-wise nnz count
        args->filternnzElems[i] = local_nnz;
        args->nnz_count[i] = local_nnz * 2;

        count += local_nnz * 2;

        // Safety check
        if (local_nnz > *args->fnnzassumed) {
            printf("WARNING: Element %d has %d neighbors. Increase fnnzassumed.\n", i, local_nnz);
            exit(EXIT_FAILURE);
        }
    }

    return NULL;
}

// Main function to build the symmetric filter matrix using pthreads
void densityfilterFast_mt(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
                         ITG *ne, double *ttime, double *timepar,
                         ITG *mortar, double *rmin, ITG *filternnz,
                         ITG *filternnzElems, ITG itertop, ITG *fnnzassumed)
{
    // Determine number of threads
    int num_threads = 1;
    char *env = getenv("OMP_NUM_THREADS");
    if (env) num_threads = atoi(env);
    if (num_threads <= 0) num_threads = 1;

    printf("Using %d threads to build filter matrix.\n", num_threads);
    ITG ne0 = *ne;
    double time = timepar[1];

    // Compute element centroids
    double *elCentroid = NULL;
    NNEW(elCentroid, double, 3 * ne0);
    mafillsmmain_filter(co, nk, *konp, *ipkonp, *lakonp, ne, ttime, &time, mortar, &ne0, elCentroid);

    // Allocate threading resources
    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    ThreadArgs *args = malloc(num_threads * sizeof(ThreadArgs));

    int elems_per_thread = (ne0 + num_threads - 1) / num_threads;

    printf("Elements per thread: %d \n", elems_per_thread);

    printf("Initializing global shared output buffers for all threads...");
    // Global shared output buffers for all threads
    int *global_drow = malloc(MAX_NNZ_PER_THREAD * sizeof(int));
    int *global_dcol = malloc(MAX_NNZ_PER_THREAD * sizeof(int));
    double *global_dval = malloc(MAX_NNZ_PER_THREAD * sizeof(double));
    int *global_nnz_count = calloc(ne0, sizeof(int));

    printf("done!\n");

    int offset = 0;
    for (int t = 0; t < num_threads; ++t) 
    {
        int start = t * elems_per_thread;
        int end = (t + 1) * elems_per_thread;
        if (end > ne0) end = ne0;

        // Slice global buffers for this thread
        args[t] = (ThreadArgs){
            .thread_id = t,
            .ne0 = ne0,
            .ne_start = start,
            .ne_end = end,
            .elCentroid = elCentroid,
            .rmin_local = *rmin,
            .fnnzassumed = fnnzassumed,
            .filternnzElems = filternnzElems,
            .drow = &global_drow[offset],
            .dcol = &global_dcol[offset],
            .dval = &global_dval[offset],
            .nnz_count = &global_nnz_count[start]
        };

        printf("Creating pthread...");
        pthread_create(&threads[t], NULL, filter_thread_worker, &args[t]);
        printf("done!\n");
    }

    
    printf("Joint threads...");
    for (int t = 0; t < num_threads; ++t)
        pthread_join(threads[t], NULL);
    printf("Done!\n");
    // Write matrix to files
    FILE *frow = fopen("drow.dat", "w");
    FILE *fcol = fopen("dcol.dat", "w");
    FILE *fval = fopen("dval.dat", "w");
    FILE *fdnnz = fopen("dnnz.dat", "w");

    int total_nnz = 0;
    for (int i = 0; i < ne0; ++i) {
        fprintf(fdnnz, "%d\n", filternnzElems[i]);
        total_nnz += global_nnz_count[i];
    }

    for (int i = 0; i < total_nnz; ++i) {
        fprintf(frow, "%d\n", global_drow[i]);
        fprintf(fcol, "%d\n", global_dcol[i]);
        fprintf(fval, "%.6f\n", global_dval[i]);
    }

    // Output summary
    *filternnz = total_nnz;

    fclose(frow); fclose(fcol); fclose(fval); fclose(fdnnz);
    SFREE(elCentroid);
    free(global_drow); free(global_dcol); free(global_dval); free(global_nnz_count);
    free(threads); free(args);
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include "CalculiX.h"

#define PROGRESS_WIDTH 40

// Global mutex for synchronized output
pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    int thread_id;
    int ne0;
    int ne_start, ne_end;
    double *elCentroid;
    double rmin_local;
    int *fnnzassumed;
    int *filternnzElems;
} ThreadArgs;

void *filter_thread_streamed(void *args_ptr) 
{
    ThreadArgs *args = (ThreadArgs *)args_ptr;
    char fname_row[64], fname_col[64], fname_val[64], fname_dnnz[64];

    // Local thread files
    sprintf(fname_row, "drow_%d.dat", args->thread_id);
    sprintf(fname_col, "dcol_%d.dat", args->thread_id);
    sprintf(fname_val, "dval_%d.dat", args->thread_id);
    sprintf(fname_dnnz, "dnnz_%d.dat", args->thread_id);

    FILE *frow = fopen(fname_row, "w");
    FILE *fcol = fopen(fname_col, "w");
    FILE *fval = fopen(fname_val, "w");
    FILE *fdnnz = fopen(fname_dnnz, "w");

    if (!frow || !fcol || !fval || !fdnnz) 
    {
        fprintf(stderr, "Thread %d: Failed to open output files.\n", args->thread_id);
        exit(EXIT_FAILURE);
    }

    int total = args->ne_end - args->ne_start;


    for (int i = args->ne_start; i < args->ne_end; ++i) 
    {
        int count = 0;
        double xi = args->elCentroid[3 * i + 0];
        double yi = args->elCentroid[3 * i + 1];
        double zi = args->elCentroid[3 * i + 2];

        for (int j = 0; j < args->ne0; ++j) 
        {
            if (i == j) continue;
            double xj = args->elCentroid[3 * j + 0];
            double yj = args->elCentroid[3 * j + 1];
            double zj = args->elCentroid[3 * j + 2];
            double dx = xi - xj, dy = yi - yj, dz = zi - zj;
            double dist = sqrt(dx * dx + dy * dy + dz * dz);

            if (dist <= args->rmin_local) 
            {
                double w = args->rmin_local - dist;
                fprintf(frow, "%d\n%d\n", i + 1, j + 1);
                fprintf(fcol, "%d\n%d\n", j + 1, i + 1);
                fprintf(fval, "%.6f\n%.6f\n", w, w);
                count++;
            }
        }


        fprintf(fdnnz, "%d\n", count);
        if (count > *args->fnnzassumed) 
        {
            printf("WARNING: Element %d has %d neighbors. Increase fnnzassumed.\n", i, count);
            exit(EXIT_FAILURE);
        }
        args->filternnzElems[i] = count;

        if ((i - args->ne_start) % 100 == 0 || i == args->ne_end - 1) 
        {
            int done = i - args->ne_start + 1;
            double percent = (double)done / total;
            int barwidth = (int)(percent * PROGRESS_WIDTH);

            char buffer[128];
            int pos = 0;
            pos += snprintf(buffer + pos, sizeof(buffer) - pos, "Thread %2d [", args->thread_id);
            for (int k = 0; k < PROGRESS_WIDTH; ++k)
                pos += snprintf(buffer + pos, sizeof(buffer) - pos, "%c", k < barwidth ? '#' : ' ');
            snprintf(buffer + pos, sizeof(buffer) - pos, "] %3.0f%%", percent * 100);

            pthread_mutex_lock(&print_mutex);
            fprintf(stderr, "\033[%d;1H", args->thread_id + 1);  // Move to line for this thread
            fprintf(stderr, "%s", buffer);                      // No newline
            fflush(stderr);
            pthread_mutex_unlock(&print_mutex);
        }
    }
    

    // Final newline after the progress bar
    pthread_mutex_lock(&print_mutex);
    fprintf(stderr, "Thread %2d completed.\n", args->thread_id);
    pthread_mutex_unlock(&print_mutex);
 

    fclose(frow); fclose(fcol); fclose(fval); fclose(fdnnz);
    return NULL;
}

void densityfilterFast_mt(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
                                ITG *ne, double *ttime, double *timepar,
                                ITG *mortar, double *rmin, ITG *filternnz,
                                ITG *filternnzElems, ITG itertop, ITG *fnnzassumed) 
                                
{
    // Determine number of threads from environment variable
    int num_threads = 1;
    char *env = getenv("OMP_NUM_THREADS");
    if (env) num_threads = atoi(env);
    if (num_threads <= 0) num_threads = 1;

    printf("\nUsing %d threads to build the filter matrix. \n", num_threads);

    
    ITG ne0 = *ne;
    double time = timepar[1];

    double *elCentroid = NULL;
    NNEW(elCentroid, double, 3 * ne0);

    // Compute element centroids
    mafillsmmain_filter(co, nk, *konp, *ipkonp, *lakonp, ne, ttime, &time, mortar, &ne0, elCentroid);

    // Allcate handles for pthread handles and thread arguments
    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    ThreadArgs *args = malloc(num_threads * sizeof(ThreadArgs));

    // Divide the total elements across threads
    int elems_per_thread = (ne0 + num_threads - 1) / num_threads;

    printf("Number of elements per thread %d \n", elems_per_thread);
    printf("Creating threads and excuting thread-local streaming...");
    // Clear screen and move cursor to top
    printf("\033[2J\033[H");
    // Launch a thread per chunk of elements
    for (int t = 0; t < num_threads; ++t) 
    {
        int start = t * elems_per_thread;
        int end = (t + 1) * elems_per_thread;
        if (end > ne0) end = ne0;

        // Set up argument struct for thread
        args[t] = (ThreadArgs)
        {
            .thread_id = t,
            .ne0 = ne0,
            .ne_start = start,
            .ne_end = end,
            .elCentroid = elCentroid,
            .rmin_local = *rmin,
            .fnnzassumed = fnnzassumed,
            .filternnzElems = filternnzElems
        };

        // Launch thread for streaming filter entries
        pthread_create(&threads[t], NULL, filter_thread_streamed, &args[t]);
    }
    printf("done. \n");

    printf("Waiting on all threads to finish...");
    // Wait for all threads to complete
    for (int t = 0; t < num_threads; ++t)
        pthread_join(threads[t], NULL);
    
    printf("done \n");

    // After threads finish, aggregate number of neighbours per element into dnnz.dat
    FILE *fdnnz = fopen("dnnz.dat", "w");
    for (int t = 0; t < num_threads; ++t) 
    {
        char fname[64];
        sprintf(fname, "dnnz_%d.dat", t);
        FILE *in = fopen(fname, "r");
        if (!in) continue;
        char line[64];
        while (fgets(line, sizeof(line), in)) fputs(line, fdnnz);
        fclose(in);
        remove(fname);
    }

    fclose(fdnnz);

    *filternnz = 0;
    for (int i = 0; i < ne0; ++i)
        *filternnz += 2 * filternnzElems[i]; // Each connection mirrored

    SFREE(elCentroid);
    free(threads); free(args);

    printf("\nThread-local filter triplet files written to disk.\n");
    printf("Merging filter triplet files from all threads...\n");
    int status = 0;

    status = system("cat drow_*.dat > drow.dat");
    if (status != 0) fprintf(stderr, "Failed to merge drow files\n");

    status = system("cat dcol_*.dat > dcol.dat");
    if (status != 0) fprintf(stderr, "Failed to merge dcol files\n");

    status = system("cat dval_*.dat > dval.dat");
    if (status != 0) fprintf(stderr, "Failed to merge dval files\n");

    // Remove intermediate thread-local files
    system("rm -f drow_*.dat dcol_*.dat dval_*.dat");

    printf("Filter matrix written: %d total nonzeros (symmetric)\n", *filternnz);
}

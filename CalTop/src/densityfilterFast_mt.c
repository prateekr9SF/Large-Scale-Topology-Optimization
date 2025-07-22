#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include "CalculiX.h"

typedef struct {
    int thread_id;
    int ne0;
    int ne_start, ne_end;
    double *elCentroid;
    double rmin_local;
    int *fnnzassumed;
    int *filternnzElems;
} ThreadArgs;

void *filter_thread_streamed(void *args_ptr) {
    ThreadArgs *args = (ThreadArgs *)args_ptr;
    char fname_row[64], fname_col[64], fname_val[64], fname_dnnz[64];

    sprintf(fname_row, "drow_%d.dat", args->thread_id);
    sprintf(fname_col, "dcol_%d.dat", args->thread_id);
    sprintf(fname_val, "dval_%d.dat", args->thread_id);
    sprintf(fname_dnnz, "dnnz_%d.dat", args->thread_id);

    FILE *frow = fopen(fname_row, "w");
    FILE *fcol = fopen(fname_col, "w");
    FILE *fval = fopen(fname_val, "w");
    FILE *fdnnz = fopen(fname_dnnz, "w");

    if (!frow || !fcol || !fval || !fdnnz) {
        fprintf(stderr, "Thread %d: Failed to open output files.\n", args->thread_id);
        exit(EXIT_FAILURE);
    }

    for (int i = args->ne_start; i < args->ne_end; ++i) {
        int count = 0;
        double xi = args->elCentroid[3 * i + 0];
        double yi = args->elCentroid[3 * i + 1];
        double zi = args->elCentroid[3 * i + 2];

        for (int j = 0; j < args->ne0; ++j) {
            if (i == j) continue;
            double xj = args->elCentroid[3 * j + 0];
            double yj = args->elCentroid[3 * j + 1];
            double zj = args->elCentroid[3 * j + 2];
            double dx = xi - xj, dy = yi - yj, dz = zi - zj;
            double dist = sqrt(dx * dx + dy * dy + dz * dz);

            if (dist <= args->rmin_local) {
                double w = args->rmin_local - dist;
                fprintf(frow, "%d\n%d\n", i + 1, j + 1);
                fprintf(fcol, "%d\n%d\n", j + 1, i + 1);
                fprintf(fval, "%.6f\n%.6f\n", w, w);
                count++;
            }
        }

        fprintf(fdnnz, "%d\n", count);
        if (count > *args->fnnzassumed) {
            printf("WARNING: Element %d has %d neighbors. Increase fnnzassumed.\n", i, count);
            exit(EXIT_FAILURE);
        }
        args->filternnzElems[i] = count;
    }

    fclose(frow); fclose(fcol); fclose(fval); fclose(fdnnz);
    return NULL;
}

void densityfilterFast_mt(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
                                ITG *ne, double *ttime, double *timepar,
                                ITG *mortar, double *rmin, ITG *filternnz,
                                ITG *filternnzElems, ITG itertop, ITG *fnnzassumed) {
    int num_threads = 1;
    char *env = getenv("OMP_NUM_THREADS");
    if (env) num_threads = atoi(env);
    if (num_threads <= 0) num_threads = 1;

    ITG ne0 = *ne;
    double time = timepar[1];

    double *elCentroid = NULL;
    NNEW(elCentroid, double, 3 * ne0);
    mafillsmmain_filter(co, nk, *konp, *ipkonp, *lakonp, ne, ttime, &time, mortar, &ne0, elCentroid);

    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    ThreadArgs *args = malloc(num_threads * sizeof(ThreadArgs));

    int elems_per_thread = (ne0 + num_threads - 1) / num_threads;

    for (int t = 0; t < num_threads; ++t) {
        int start = t * elems_per_thread;
        int end = (t + 1) * elems_per_thread;
        if (end > ne0) end = ne0;

        args[t] = (ThreadArgs){
            .thread_id = t,
            .ne0 = ne0,
            .ne_start = start,
            .ne_end = end,
            .elCentroid = elCentroid,
            .rmin_local = *rmin,
            .fnnzassumed = fnnzassumed,
            .filternnzElems = filternnzElems
        };

        pthread_create(&threads[t], NULL, filter_thread_streamed, &args[t]);
    }

    for (int t = 0; t < num_threads; ++t)
        pthread_join(threads[t], NULL);

    FILE *fdnnz = fopen("dnnz.dat", "w");
    for (int t = 0; t < num_threads; ++t) {
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
        *filternnz += 2 * filternnzElems[i];

    SFREE(elCentroid);
    free(threads); free(args);
}

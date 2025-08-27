#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include "CalculiX.h"

#include <stdint.h> 

#define PROGRESS_WIDTH 40

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

// --- Small helper: copy a whole file (binary) into an already open out fp ---
static void copy_file_into(FILE *out, const char *path, size_t buf_bytes) {
    FILE *in = fopen(path, "rb");
    if (!in) return;
    if (buf_bytes < 1) buf_bytes = 1<<20; // 1 MB default
    unsigned char *buf = (unsigned char*)malloc(buf_bytes);
    if (!buf) { fclose(in); return; }
    size_t nread;
    while ((nread = fread(buf, 1, buf_bytes, in)) > 0) {
        fwrite(buf, 1, nread, out);
    }
    free(buf);
    fclose(in);
}


/**
 * Builds a symmetric density filter matrix (in triplet format) in parallel.
 * Each thread processes a block of elements, computes filter weights
 * based on distance to neighbors within a radius `rmin_local`, and writes
 * output to thread-local files: drow_*.dat, dcol_*.dat, dval_*.dat, dnnz_*.dat.
 */
void *filter_thread_bin(void *args_ptr)
{
    ThreadArgs *args = (ThreadArgs *)args_ptr;

    // Prepare thread-local filenames
    char fname_row[64], fname_col[64], fname_val[64], fname_dnnz[64], fname_dsum[64];

    sprintf(fname_row,  "drow_%d.bin", args->thread_id);
    sprintf(fname_col,  "dcol_%d.bin", args->thread_id);
    sprintf(fname_val,  "dval_%d.bin", args->thread_id);
    sprintf(fname_dnnz, "dnnz_%d.bin", args->thread_id);
    sprintf(fname_dsum, "dsum_%d.bin", args->thread_id);

    FILE *frow = fopen(fname_row, "wb");
    FILE *fcol = fopen(fname_col, "wb");
    FILE *fval = fopen(fname_val, "wb");
    FILE *fdnnz = fopen(fname_dnnz, "wb");


    if (!frow || !fcol || !fval || !fdnnz) 
    {
        pthread_mutex_lock(&print_mutex);
        fprintf(stderr, "Thread %d: Failed to open output files.\n", args->thread_id);
        pthread_mutex_unlock(&print_mutex);
        exit(EXIT_FAILURE);
    }

    // Thread-local row-sum accumulator (full length ne0)
    double *row_sum_local = (double*)calloc((size_t)args->ne0, sizeof(double));

    if (!row_sum_local) 
    {
        pthread_mutex_lock(&print_mutex);
        fprintf(stderr, "Thread %d: Failed to allocate row_sum_local.\n", args->thread_id);
        pthread_mutex_unlock(&print_mutex);
        exit(EXIT_FAILURE);
    }

    int total = args->ne_end - args->ne_start;

    // Loop over this thread's assigned element block
    for (int i = args->ne_start; i < args->ne_end; ++i) 
    {
        int count = 0;

        // Get centroid of element i
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

            // Compute distance
            double dist = sqrt(dx * dx + dy * dy + dz * dz);
            
            // If within filter radius, perform evals
            if (dist <= args->rmin_local) 
            {
                double w = args->rmin_local - dist;

                // Write symmetric pair in binary (1-based ints)(64-bin indices)                
                int64_t row1 = (int64_t)i + 1, col1 = (int64_t)j + 1;
                int64_t row2 = (int64_t)j + 1, col2 = (int64_t)i + 1;
                fwrite(&row1, sizeof(int64_t), 1, frow);
                fwrite(&row2, sizeof(int64_t), 1, frow);
                fwrite(&col1, sizeof(int64_t), 1, fcol);
                fwrite(&col2, sizeof(int64_t), 1, fcol);

                fwrite(&w, sizeof(double), 1, fval);
                fwrite(&w, sizeof(double), 1, fval);

               // Accumulate row-wise sums for symmetric H:
                // row i gains w for (i,j), row j gains w for (j,i)
                row_sum_local[i] += w;  // row i gets w for (i,j)
                row_sum_local[j] += w;  // row j gets w for (j,i)

                count++;
            }
        }

        // Write number of neighbors for element i
        //fprintf(fdnnz, "%d\n", count);

        // Per-row neighbour count (for i-only -- excludes mirroring )
        fwrite (&count, sizeof(int), 1, fdnnz);

        // Guard: check if fnnzassumed is too small
        if (count > *args->fnnzassumed) 
        {
            pthread_mutex_lock(&print_mutex);
            printf("WARNING: Element %d has %d neighbors. Increase fnnzassumed.\n", i, count);
            pthread_mutex_unlock(&print_mutex);
            exit(EXIT_FAILURE);
        }

        // Store neighbor count for postprocessing
        args->filternnzElems[i] = count;

        // Progress bar update every 100 elements
        if ((i - args->ne_start) % 100 == 0 || i == args->ne_end - 1) 
        {
            int done = i - args->ne_start + 1;
            double percent = (double)done / total;
            int barwidth = (int)(percent * PROGRESS_WIDTH);

            char buffer[128];
            int pos = 0;
            pos += snprintf(buffer + pos, sizeof(buffer), "Thread %2d [", args->thread_id);
            for (int k = 0; k < PROGRESS_WIDTH; ++k)
                pos += snprintf(buffer + pos, sizeof(buffer) - pos, "%c", k < barwidth ? '#' : ' ');
            snprintf(buffer + pos, sizeof(buffer) - pos, "] %3.0f%%", percent * 100);

            pthread_mutex_lock(&print_mutex);
            fprintf(stderr, "\033[%d;1H\033[2K%s", args->thread_id + 1, buffer);
            fflush(stderr);
            pthread_mutex_unlock(&print_mutex);
        }
    }

    
    pthread_mutex_lock(&print_mutex);
    fprintf(stderr, "\033[%d;1H\033[2K\033[32mThread %2d [########################################] 100%% - completed\033[0m", args->thread_id + 1, args->thread_id);
    fflush(stderr);
    pthread_mutex_unlock(&print_mutex);

    // Write thread-local rown sums (length = ne0) to binary file
    FILE *fdsum = fopen(fname_dsum, "wb");
    if (!fdsum) 
    {
        pthread_mutex_lock(&print_mutex);
        fprintf(stderr, "Thread %d: Failed to open %s for write.\n", args->thread_id, fname_dsum);
        pthread_mutex_unlock(&print_mutex);
        exit(EXIT_FAILURE);
    }


    fwrite(row_sum_local, sizeof(double), (size_t)args->ne0, fdsum);
    fclose(fdsum);


    // Close this thread's files
    fclose(frow); 
    fclose(fcol); 
    fclose(fval); 
    fclose(fdnnz);

    free(row_sum_local);

    return NULL;
}



/* Top-level driver: launches threads that stream/write triplets and per-row nnz,
   then merges thread-local files into drow.dat, dcol.dat, dval.dat, dnnz.dat,
   and reduces row sums into dsum.dat. */
void densityfilterFast_bin_mt(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
                          ITG *ne, double *ttime, double *timepar,
                          ITG *mortar, double *rmin, ITG *filternnz,
                          ITG *filternnzElems, ITG itertop, ITG *fnnzassumed) 
{
    int num_threads = 1;
    char *env = getenv("OMP_NUM_THREADS");
    if (env) num_threads = atoi(env);
    if (num_threads <= 0) num_threads = 1;

    printf("\033[2J\033[H"); // Clear screen
    printf("Using %d threads to build the filter matrix.\n", num_threads);

    ITG ne0 = *ne;
    double time = timepar[1];

    double *elCentroid = NULL;
    NNEW(elCentroid, double, 3 * ne0);

    mafillsmmain_filter(co, nk, *konp, *ipkonp, *lakonp, ne, ttime, &time, mortar, &ne0, elCentroid);

    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    ThreadArgs *args = malloc(num_threads * sizeof(ThreadArgs));

    int elems_per_thread = (ne0 + num_threads - 1) / num_threads;

    for (int t = 0; t < num_threads; ++t) printf("\n");
    printf("\033[%d;1H", num_threads + 2);
    fflush(stdout);

    printf("Number of elements per thread: %d\n", elems_per_thread);
    printf("Creating threads and streaming thread-local binary shards...\n");

    for (int t = 0; t < num_threads; ++t) 
    {
        int start = t * elems_per_thread;
        int end = (t + 1) * elems_per_thread;
        if (end > ne0) end = ne0;

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

        pthread_create(&threads[t], NULL, filter_thread_bin, &args[t]);
    }

    for (int t = 0; t < num_threads; ++t)
        pthread_join(threads[t], NULL);

    // Move below progress bar region before printing
    printf("\033[%d;1H", num_threads + 5);
    fflush(stdout);
    printf("All threads completed.\n");

    // Compute filternnz (symmetric): 2 * sum of per-row neighbor counts
    long long sum_dnnz = 0;
    for (int i = 0; i < ne0; ++i) sum_dnnz += filternnzElems[i];
    *filternnz = (ITG)(2LL * sum_dnnz);

    //SFREE(elCentroid);
    //free(threads); free(args);

    printf("\n\033[36m==================== Filter Matrix Summary ====================\033[0m\n");
    //printf("Thread-local filter triplet files written to disk.\n");
    //printf("Merging filter triplet files from all threads...\n");

    // Begin merging binary files
    printf("Merging thread-local binaries...");
    FILE *Fdrow = fopen("drow.bin","wb");
    FILE *Fdcol = fopen("dcol.bin","wb");
    FILE *Fdval = fopen("dval.bin","wb");
    FILE *Fdnnz = fopen("dnnz.bin","wb");

    if (!Fdrow||!Fdcol||!Fdval||!Fdnnz) 
    {
        fprintf(stderr,"Failed to open final binary files.\n");
        exit(EXIT_FAILURE);
    }


    // Merge drow/dcol/dval
    for (int t = 0; t < num_threads; ++t) {
        char fpath[64];
        sprintf(fpath, "drow_%d.bin", t); copy_file_into(Fdrow, fpath, 8<<20); remove(fpath);
        sprintf(fpath, "dcol_%d.bin", t); copy_file_into(Fdcol, fpath, 8<<20); remove(fpath);
        sprintf(fpath, "dval_%d.bin", t); copy_file_into(Fdval, fpath, 8<<20); remove(fpath);
    }
    fclose(Fdrow); fclose(Fdcol); fclose(Fdval);

        // --- Preview: first 5 entries of dval.bin (after merge) ---
    {
        FILE *fin = fopen("dval.bin", "rb");
        if (fin) {
            double buf[5];
            size_t got = fread(buf, sizeof(double), 5, fin);
            fclose(fin);
            int nshow = (got < 5) ? (int)got : 5;
            printf("First %d entries of dval.bin: ", nshow);
            for (int i = 0; i < nshow; ++i) {
                printf("%.10g%s", buf[i], (i == nshow - 1) ? "\n" : ", ");
            }
        } else {
            fprintf(stderr, "Warning: could not reopen dval.bin for preview\n");
        }
    }

    // Merge dnnz (already contiguous per thread)
    for (int t = 0; t < num_threads; ++t) {
        char fpath[64];
        sprintf(fpath, "dnnz_%d.bin", t); copy_file_into(Fdnnz, fpath, 1<<20); remove(fpath);
    }
    fclose(Fdnnz);

    // Reduce dsum_T.bin → dsum.bin (sum across threads)
    printf("Reducing row-wise sums into dsum.bin...\n");
    double *dsum = (double*)calloc((size_t)ne0, sizeof(double));
    if (!dsum) { fprintf(stderr, "Failed to allocate dsum accumulator.\n"); exit(EXIT_FAILURE); }

    double *tmpbuf = (double*)malloc((size_t)ne0 * sizeof(double));
    if (!tmpbuf) { fprintf(stderr,"Failed to allocate tmpbuf.\n"); free(dsum); exit(EXIT_FAILURE); }

    for (int t = 0; t < num_threads; ++t) {
        char fpath[64];
        sprintf(fpath, "dsum_%d.bin", t);
        FILE *in = fopen(fpath, "rb");
        if (!in) { fprintf(stderr,"Warning: missing %s; treating as zeros\n", fpath); continue; }
        size_t got = fread(tmpbuf, sizeof(double), (size_t)ne0, in);
        if (got != (size_t)ne0) {
            fprintf(stderr,"Error: short read in %s (%zu of %d)\n", fpath, got, (int)ne0);
            fclose(in); free(tmpbuf); free(dsum); exit(EXIT_FAILURE);
        }
        fclose(in);
        remove(fpath);
        for (int r = 0; r < ne0; ++r) dsum[r] += tmpbuf[r];
    }
    free(tmpbuf);

    // --- Preview: first 5 entries of dsum (after merge/reduction) ---
    {
    int nshow = (ne0 < 5) ? ne0 : 5;
    printf("First %d entries of dsum: ", nshow);
    for (int i = 0; i < nshow; ++i) {
        printf("%.10g%s", dsum[i], (i == nshow - 1) ? "\n" : ", ");
    }
    }

    FILE *Fdsum = fopen("dsum.bin","wb");
    if (!Fdsum) { fprintf(stderr,"Failed to open dsum.bin for write.\n"); free(dsum); exit(EXIT_FAILURE); }
    fwrite(dsum, sizeof(double), (size_t)ne0, Fdsum);
    fclose(Fdsum);
    free(dsum);

    // Cleanup
    SFREE(elCentroid);
    free(threads);
    free(args);

    printf("Done!\n");

    printf("Final binaries: drow.bin, dcol.bin, dval.bin, dnnz.bin, dsum.bin\n");
    printf("Filter matrix nnz (symmetric): %d\n", *filternnz);
    printf("===============================================================\n");
}
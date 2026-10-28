/* CalFilt: disk-backed CSR density filter assembly.
 * Zero-based CSR columns and row pointers. Requires CalculiX.h and pthreads.
 * Retains original O(ne^2) centroid search for controlled validation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include "CalculiX.h"
#ifdef PROFILING_ON
#include <TAU.h>
#endif

#define PROGRESS_WIDTH 40
#define PROGRESS_INTERVAL 100

static pthread_mutex_t progress_mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    int id, ne, begin, end, excludePassive;
    const int *passive;
    const double *centroid;
    double radius;
    ITG maxNeighbors;
    ITG *counts;
    double *sums;
    int failed;
} CSRWorker;

static int copy_shard(FILE *out, const char *name) 
{
    unsigned char buffer[1 << 20];
    FILE *in = fopen(name, "rb");
    if (!in) { perror(name); return -1; }
    size_t n;

    while ((n = fread(buffer, 1, sizeof buffer, in)) != 0) 
    {
        if (fwrite(buffer, 1, n, out) != n) { perror("write merged CSR"); fclose(in); return -1; }
    }

    if (ferror(in)) { perror("read CSR shard"); fclose(in); return -1; }
    if (fclose(in) != 0) return -1;

    /* Shards are removed only after final CSR output files close successfully. */
    return 0;
}

static void print_csr_progress(int thread_id, int done, int total)
{
    if (total <= 0) return;

    double percent = 100.0 * (double)done / (double)total;
    int filled = (int)(PROGRESS_WIDTH * percent / 100.0);

    pthread_mutex_lock(&progress_mutex);

    //fprintf(stderr, "\rThread %2d [", thread_id);

    fprintf(stderr, "\033[%d;1H\033[2KThread %2d [",
        thread_id + 1, thread_id);

    for (int k = 0; k < PROGRESS_WIDTH; ++k)
        fputc(k < filled ? '#' : ' ', stderr);

    fprintf(stderr, "] %6.2f%% (%d/%d)",
            percent, done, total);

    if (done == total)
        fprintf(stderr, " - completed\n");

    fflush(stderr);

    pthread_mutex_unlock(&progress_mutex);
}

static void *csr_worker(void *p) 
{
    CSRWorker *a = (CSRWorker *)p;

#ifdef PROFILING_ON
    TAU_REGISTER_THREAD();
    TAU_PROFILE_TIMER(timer_filter_thread, "densityFilter_eval_CSR", "", TAU_USER);
    TAU_PROFILE_START(timer_filter_thread);
#endif

    char colname[80], valname[80];
    snprintf(colname, sizeof colname, "csr_col_%d.tmp", a->id);
    snprintf(valname, sizeof valname, "csr_val_%d.tmp", a->id);
    FILE *fc = fopen(colname, "wb"), *fv = fopen(valname, "wb");
    if (!fc || !fv) { perror("open CSR shard"); a->failed = 1; goto done; }
    for (int i = a->begin; i < a->end; ++i) {
        int64_t count = 0;
        double sum = 0.0;
        if (a->excludePassive && a->passive[i]) {
            int64_t j = i;
            double w = 1.0;
            if (fwrite(&j, sizeof j, 1, fc) != 1 || fwrite(&w, sizeof w, 1, fv) != 1) { a->failed = 1; break; }
            count = 1;
            sum = 1.0;
        }
        else
        {
            const double xi = a->centroid[3*(size_t)i];
            const double yi = a->centroid[3*(size_t)i+1];
            const double zi = a->centroid[3*(size_t)i+2];
            
            for (int j = 0; j < a->ne; ++j) 
            {
                if (a->excludePassive && a->passive[j]) continue;
                double dx = xi - a->centroid[3*(size_t)j];
                double dy = yi - a->centroid[3*(size_t)j+1];
                double dz = zi - a->centroid[3*(size_t)j+2];
                double dist2 = dx*dx + dy*dy + dz*dz;
                if (dist2 > a->radius*a->radius) continue;

                double w = a->radius - sqrt(dist2);
                /* A zero-weight boundary neighbor is not a nonzero. */
                if (w <= 0.0) continue;
                int64_t col = j;
                if (fwrite(&col, sizeof col, 1, fc) != 1 || fwrite(&w, sizeof w, 1, fv) != 1) { a->failed = 1; break; }
                sum += w;
                ++count;
            }
            if (a->failed) break;
        }

        if (count > a->maxNeighbors || count > INT_MAX) 
        {
            fprintf(stderr, "CSR row %d has %" PRId64 " neighbors; fnnzassumed=%lld\n", i+1, count, (long long)a->maxNeighbors);
            a->failed = 1; break;
        }
        a->counts[i] = (ITG)count;
        a->sums[i] = sum;

        int completed = i - a->begin + 1;
        int total = a->end - a->begin;

        if (completed % PROGRESS_INTERVAL == 0 ||completed == total)
        {
            print_csr_progress(a->id, completed, total);
        }

    }
done:
    if (fc && fclose(fc) != 0) a->failed = 1;
    if (fv && fclose(fv) != 0) a->failed = 1;
#ifdef PROFILING_ON
    TAU_PROFILE_STOP(timer_filter_thread);
#endif
    return NULL;
}

void densityfilterFast_bin_mt_csr(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
                              ITG *ne, double *ttime, double *timepar, ITG *mortar,
                              double *rmin, ITG *filternnz, ITG *filternnzElems,
                              ITG itertop, ITG *fnnzassumed,
                              int *passiveIDs, int numPassive, int excludePassive) {
    (void)itertop;
    if (*ne <= 0 || *rmin <= 0 || *fnnzassumed <= 0) {
        fprintf(stderr, "Invalid CSR assembly parameters\n"); exit(EXIT_FAILURE);
    }
    const int n = (int)*ne;
    if ((ITG)n != *ne) { fprintf(stderr, "Element count exceeds int range\n"); exit(EXIT_FAILURE); }
    int nt = 1;
    const char *env = getenv("OMP_NUM_THREADS");
    if (env) nt = atoi(env);
    if (nt < 1) nt = 1;
    if (nt > n) nt = n;
    int *passive = calloc((size_t)n, sizeof *passive);
    double *centroid = NULL;
    double *sums = calloc((size_t)n, sizeof *sums);
    int64_t *rowptr = malloc(((size_t)n + 1) * sizeof *rowptr);
    pthread_t *threads = calloc((size_t)nt, sizeof *threads);
    CSRWorker *args = calloc((size_t)nt, sizeof *args);
    if (!passive || !sums || !rowptr || !threads || !args) { perror("CSR allocation"); exit(EXIT_FAILURE); }
    for (int p = 0; p < numPassive; ++p) {
        int id = passiveIDs[p] - 1;
        if (id >= 0 && id < n) passive[id] = 1;
    }

    ITG ne0 = *ne;
    double time = timepar[1];
    NNEW(centroid, double, 3 * (size_t)n);
    mafillsmmain_filter(co, nk, *konp, *ipkonp, *lakonp, ne, ttime, &time, mortar, &ne0, centroid);

    fprintf(stderr, "\033[2J\033[H");

    for (int t = 0; t < nt; ++t)
        fprintf(stderr, "\n");

    fflush(stderr);

    int started = 0;
    for (int t = 0; t < nt; ++t) {
        args[t] = (CSRWorker){ .id=t, .ne=n, .begin=(int)((int64_t)n*t/nt),
            .end=(int)((int64_t)n*(t+1)/nt), .excludePassive=excludePassive,
            .passive=passive, .centroid=centroid, .radius=*rmin,
            .maxNeighbors=*fnnzassumed, .counts=filternnzElems, .sums=sums };
        if (pthread_create(&threads[t], NULL, csr_worker, &args[t]) != 0) {
            fprintf(stderr, "pthread_create failed\n");
            for (int j = 0; j < started; ++j) pthread_join(threads[j], NULL);
            exit(EXIT_FAILURE);
        }
        ++started;
    }
    int failed = 0;
    for (int t = 0; t < started; ++t) {
        pthread_join(threads[t], NULL);
        failed |= args[t].failed;
    }
    fprintf(stderr, "\033[%d;1H\n", nt + 2);
    fflush(stderr);
    
    if (failed) { fprintf(stderr, "CSR assembly failed; shards retained for diagnosis\n"); exit(EXIT_FAILURE); }
    rowptr[0] = 0;
    for (int i = 0; i < n; ++i) {
        if ((int64_t)filternnzElems[i] > INT64_MAX - rowptr[i]) {
            fprintf(stderr, "CSR nnz overflow\n"); exit(EXIT_FAILURE);
        }
        rowptr[i+1] = rowptr[i] + (int64_t)filternnzElems[i];
    }
    /* The CSR files, not the legacy ITG output, are authoritative for nnz.
     * Never narrow a value that exceeds the current ITG range.
     * CSR consumers must read row_ptr[ne] and must not use filternnz as nnz.
     */
    if (sizeof(ITG) == 4 && rowptr[n] > INT32_MAX) {
        *filternnz = (ITG)0;  /* Sentinel: nnz must be read from row_ptr.bin. */
        fprintf(stderr,
                "CSR nnz=%" PRId64 " exceeds 32-bit ITG; "
                "legacy filternnz set to 0 (CSR readers use row_ptr.bin)\n",
                rowptr[n]);
    } else {
        *filternnz = (ITG)rowptr[n];
    }
    FILE *fr = fopen("row_ptr.bin", "wb");
    FILE *fc = fopen("col_idx.bin", "wb");
    FILE *fv = fopen("dval.bin", "wb");
    FILE *fs = fopen("dsum.bin", "wb");
    FILE *fn = fopen("dnnz.bin", "wb");
    if (!fr || !fc || !fv || !fs || !fn) { perror("open CSR output"); exit(EXIT_FAILURE); }
    if (fwrite(rowptr, sizeof *rowptr, (size_t)n+1, fr) != (size_t)n+1 ||
        fwrite(sums, sizeof *sums, (size_t)n, fs) != (size_t)n) {
        perror("write CSR metadata"); exit(EXIT_FAILURE);
    }
    /* Preserve the original dnnz.bin 32-bit integer layout. */
    for (int i = 0; i < n; ++i) {
        int32_t count = (int32_t)filternnzElems[i];
        if (fwrite(&count, sizeof count, 1, fn) != 1) { perror("write dnnz.bin"); exit(EXIT_FAILURE); }
    }
    for (int t = 0; t < nt; ++t) {
        char name[80];
        snprintf(name, sizeof name, "csr_col_%d.tmp", t);
        if (copy_shard(fc, name)) exit(EXIT_FAILURE);
        snprintf(name, sizeof name, "csr_val_%d.tmp", t);
        if (copy_shard(fv, name)) exit(EXIT_FAILURE);
    }
    int close_error = 0;
    close_error |= fclose(fr); close_error |= fclose(fc); close_error |= fclose(fv);
    close_error |= fclose(fs); close_error |= fclose(fn);
    if (close_error) {
        perror("close CSR output"); exit(EXIT_FAILURE);
    }
    /* Final output files have been closed successfully; remove temporary shards. */
    for (int t = 0; t < nt; ++t) {
        char name[80];
        snprintf(name, sizeof name, "csr_col_%d.tmp", t);
        if (remove(name) != 0) perror("remove CSR column shard");
        snprintf(name, sizeof name, "csr_val_%d.tmp", t);
        if (remove(name) != 0) perror("remove CSR value shard");
    }
    printf("CSR filter assembled: ne=%d nnz=%" PRId64 " threads=%d\n", n, rowptr[n], nt);
    printf("Files: row_ptr.bin col_idx.bin dval.bin dsum.bin dnnz.bin\n");
    free(passive); free(sums); free(rowptr); free(threads); free(args);
    SFREE(centroid);
}

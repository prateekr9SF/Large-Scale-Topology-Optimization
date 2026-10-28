/* CalTop: bounded-memory, pthread-parallel CSR transpose sensitivity filtering.
 * CSR files (all zero-based): row_ptr.bin: int64_t[ne+1],
 * col_idx.bin: int64_t[nnz], dval.bin: double[nnz],
 * dsum.bin: double[ne].
 *
 * Compile as a separate source alongside the existing COO implementation.
 * Uses POSIX fseeko for offsets larger than 2 GiB.
 */
#define _FILE_OFFSET_BITS 64
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#ifdef PROFILING_ON
#include <TAU.h>
#endif

#define CSR_CHUNK_ENTRIES 65536U

typedef struct {
    int id, begin, end, ne;
    const int64_t *rowptr;
    const double *input;
    double *output;
    int failed;
} CSRFilterWorker;

static void csr_fatal(const char *msg)
{
    fprintf(stderr, "CSR sensitivity filter: %s\n", msg);
    exit(EXIT_FAILURE);
}

static int64_t csr_file_items(const char *name, size_t item_size)
{
    struct stat st;
    if (stat(name, &st) != 0) { perror(name); exit(EXIT_FAILURE); }
    if (st.st_size < 0 || (uint64_t)st.st_size % item_size != 0)
        csr_fatal("invalid binary file size");
    return (int64_t)((uint64_t)st.st_size / item_size);
}

static void csr_read_exact(const char *name, void *dst, size_t item_size, size_t n)
{
    FILE *f = fopen(name, "rb");
    if (!f) { perror(name); exit(EXIT_FAILURE); }
    if (fread(dst, item_size, n, f) != n) {
        fprintf(stderr, "CSR sensitivity filter: short read of %s\n", name);
        fclose(f); exit(EXIT_FAILURE);
    }
    if (fclose(f) != 0) { perror(name); exit(EXIT_FAILURE); }
}

static void *csr_filter_worker(void *opaque)
{
    CSRFilterWorker *a = (CSRFilterWorker *)opaque;
#ifdef PROFILING_ON
    TAU_REGISTER_THREAD();
    TAU_PROFILE_TIMER(t_csr_worker, "CSR Sensitivity Filter Worker", "", TAU_USER);
    TAU_PROFILE_START(t_csr_worker);
#endif
    FILE *fc = NULL, *fv = NULL;
    int64_t *cols = NULL;
    double *weights = NULL;

    fc = fopen("col_idx.bin", "rb");
    fv = fopen("dval.bin", "rb");
    cols = (int64_t *)malloc(CSR_CHUNK_ENTRIES * sizeof(*cols));
    weights = (double *)malloc(CSR_CHUNK_ENTRIES * sizeof(*weights));
    if (!fc || !fv || !cols || !weights) {
        fprintf(stderr, "CSR worker %d: unable to open files or allocate buffers\n", a->id);
        a->failed = 1;
        goto cleanup;
    }

    /* Each worker reads a contiguous region of nonzeros and owns its output rows. */
    int64_t first = a->rowptr[a->begin];
    if ((uint64_t)first > (uint64_t)INT64_MAX / sizeof(int64_t)) {
        a->failed = 1; goto cleanup;
    }
    off_t offset = (off_t)(first * (int64_t)sizeof(int64_t));
    if ((int64_t)offset != first * (int64_t)sizeof(int64_t) ||
        fseeko(fc, offset, SEEK_SET) != 0 ||
        fseeko(fv, offset, SEEK_SET) != 0) {
        fprintf(stderr, "CSR worker %d: seek failed\n", a->id);
        a->failed = 1; goto cleanup;
    }

    size_t buffered = 0, pos = 0;
    int64_t cursor = first;
    const int64_t last = a->rowptr[a->end];

    for (int i = a->begin; i < a->end; ++i) {
        double sum = 0.0;
        const int64_t row_end = a->rowptr[i+1];
        while (cursor < row_end) {
            if (pos == buffered) {
                int64_t remain = last - cursor;
                size_t take = (remain < (int64_t)CSR_CHUNK_ENTRIES)
                              ? (size_t)remain : CSR_CHUNK_ENTRIES;
                if (take == 0 || fread(cols, sizeof(*cols), take, fc) != take ||
                    fread(weights, sizeof(*weights), take, fv) != take) {
                    fprintf(stderr, "CSR worker %d: truncated CSR data near nnz=%" PRId64 "\n",
                            a->id, cursor);
                    a->failed = 1; goto cleanup;
                }
                buffered = take;
                pos = 0;
            }
            size_t take = buffered - pos;
            int64_t row_remaining = row_end - cursor;
            if ((int64_t)take > row_remaining) take = (size_t)row_remaining;
            for (size_t k = 0; k < take; ++k) {
                int64_t j = cols[pos+k];
                if (j < 0 || j >= a->ne || !isfinite(weights[pos+k])) {
                    fprintf(stderr, "CSR worker %d: invalid column/weight at nnz=%" PRId64 "\n",
                            a->id, cursor+(int64_t)k);
                    a->failed = 1; goto cleanup;
                }
                sum += weights[pos+k] * a->input[(size_t)j];
            }
            pos += take;
            cursor += (int64_t)take;
        }
        a->output[i] = sum;
    }
cleanup:
    if (fc && fclose(fc) != 0) a->failed = 1;
    if (fv && fclose(fv) != 0) a->failed = 1;
    free(cols);
    free(weights);
#ifdef PROFILING_ON
    TAU_PROFILE_STOP(t_csr_worker);
#endif
    return NULL;
}

/* Compute H^T D^{-1} g = H D^{-1} g.
 * Valid only because the assembled unnormalized H is symmetric.
 * Keep the COO driver's int-based signature for initial validation.
 */
void filterSensitivity_buffered_bin_mt_csr(double *Vector, double *VectorFiltered,
                                       int *filternnzElems, int *ne_ptr,
                                       int *fnnzassumed_ptr, double *q_ptr,
                                       int filternnz_total)
{
#ifdef PROFILING_ON
    TAU_PROFILE_TIMER(t_csr_total, "CSR Sensitivity Filter: Total", "", TAU_USER);
    TAU_PROFILE_START(t_csr_total);
#endif
    (void)filternnzElems;
    (void)fnnzassumed_ptr;
    (void)q_ptr; /* Use stored dval directly. */
    (void)filternnz_total; /* Use row_ptr.bin's actual nnz count instead. */
    if (!Vector || !VectorFiltered || !ne_ptr || *ne_ptr <= 0)
        csr_fatal("invalid input arrays or ne");
    if (Vector == VectorFiltered)
        csr_fatal("input and output must not alias");
    const int ne = *ne_ptr;
    if (csr_file_items("row_ptr.bin", sizeof(int64_t)) != (int64_t)ne + 1 ||
        csr_file_items("dsum.bin", sizeof(double)) != ne)
        csr_fatal("row_ptr.bin or dsum.bin dimensions do not match ne");

    int64_t *rowptr = malloc(((size_t)ne+1) * sizeof(*rowptr));
    double *dsum = malloc((size_t)ne * sizeof(*dsum));
    if (!rowptr || !dsum) csr_fatal("metadata allocation failed");
    csr_read_exact("row_ptr.bin", rowptr, sizeof(*rowptr), (size_t)ne+1);
    csr_read_exact("dsum.bin", dsum, sizeof(*dsum), (size_t)ne);

    /* H is symmetric for the assembled centroid-distance filter.
     * Apply D^{-1} to the input BEFORE multiplying by H.
     */
    double *scaled = malloc((size_t)ne * sizeof(*scaled));
    if (!scaled) csr_fatal("scaled sensitivity allocation failed");
    for (int i = 0; i < ne; ++i) {
        if (!(dsum[i] > 0.0) || !isfinite(dsum[i]) || !isfinite(Vector[i]))
            csr_fatal("invalid normalization or input sensitivity");
        scaled[i] = Vector[i] / dsum[i];
    }

    if (rowptr[0] != 0) csr_fatal("row_ptr[0] is not zero");
    for (int i = 0; i < ne; ++i) {
        if (rowptr[i+1] <= rowptr[i])
            csr_fatal("row pointers not strictly increasing (empty row)");
    }
    int64_t nnz = rowptr[ne];
    if (csr_file_items("col_idx.bin", sizeof(int64_t)) != nnz ||
        csr_file_items("dval.bin", sizeof(double)) != nnz)
        csr_fatal("CSR column/weight file lengths do not match row_ptr[ne]");

    int nt = 4;
    const char *env = getenv("OMP_NUM_THREADS");
    if (env && *env) {
        char *end = NULL;
        long requested = strtol(env, &end, 10);
        if (end != env && *end == '\0' && requested > 0 && requested <= INT_MAX)
            nt = (int)requested;
    }
    if (nt > ne) nt = ne;

    pthread_t *threads = calloc((size_t)nt, sizeof(*threads));
    CSRFilterWorker *args = calloc((size_t)nt, sizeof(*args));
    if (!threads || !args) csr_fatal("thread allocation failed");

    /* Partition by nonzero count, not row count, to balance uneven sparsity. */
    int *bounds = malloc(((size_t)nt+1) * sizeof(*bounds));
    if (!bounds) csr_fatal("partition allocation failed");
    bounds[0] = 0;
    bounds[nt] = ne;
    for (int t = 1; t < nt; ++t) {
        int64_t target = (nnz / nt) * t + (nnz % nt) * t / nt;
        int lo = bounds[t-1], hi = ne;
        while (lo < hi) {
            int mid = lo + (hi-lo)/2;
            if (rowptr[mid] < target) lo = mid+1;
            else hi = mid;
        }
        bounds[t] = lo;
    }

    int started = 0;
    for (int t = 0; t < nt; ++t) {
        args[t] = (CSRFilterWorker){.id=t, .begin=bounds[t], .end=bounds[t+1],
            .ne=ne, .rowptr=rowptr,
            .input=scaled, .output=VectorFiltered, .failed=0};
        int rc = pthread_create(&threads[t], NULL, csr_filter_worker, &args[t]);
        if (rc != 0) {
            fprintf(stderr, "CSR pthread_create: %s\n", strerror(rc));
            for (int j = 0; j < started; ++j) pthread_join(threads[j], NULL);
            csr_fatal("failed to create worker threads");
        }
        ++started;
    }
    int failed = 0;
    for (int t = 0; t < started; ++t) {
        int rc = pthread_join(threads[t], NULL);
        if (rc != 0) { fprintf(stderr, "CSR pthread_join: %s\n", strerror(rc)); failed = 1; }
        failed |= args[t].failed;
    }
    free(bounds);
    free(args);
    free(threads);
    free(dsum);
    free(scaled);
    free(rowptr);
    if (failed) csr_fatal("worker failure");
#ifdef PROFILING_ON
    TAU_PROFILE_STOP(t_csr_total);
#endif
}

/* Drop-in CSR replacement for the legacy COO stress sensitivity filter.
 * SensIn  = dPnorm_drho (derivative w.r.t. physical densities)
 * SensOut = dPnorm_drhoFiltered (derivative w.r.t. design densities)
 * nnz_total is retained for compatibility; row_ptr[ne] is authoritative.
 */
void filterSensitivity_bin_buffered_mts_csr(const double *SensIn,
                                            double *SensOut,
                                            int ne,
                                            long long nnz_total)
{
    int unused_nnz_assumed = 0;
    double qfilter = 1.0;
    (void)nnz_total;
    filterSensitivity_buffered_bin_mt_csr((double *)SensIn, SensOut,
                                          NULL, &ne, &unused_nnz_assumed,
                                          &qfilter, 0);
}

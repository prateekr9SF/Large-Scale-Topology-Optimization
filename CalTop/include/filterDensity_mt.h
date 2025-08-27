#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

#define BLOCK_SIZE 100000000

/**
 * @brief Struct to pass arguments to each pthread for vector filtering using atomic updates.
 *
 * Each thread processes a portion of the filter triplets, accumulating results
 * into shared global arrays using atomic additions.
 */
typedef struct {
    int thread_id;            ///< ID of the current thread (0 to num_threads-1)
    int num_threads;          ///< Total number of threads
    int64_t *drow;                ///< Row indices of filter matrix
    int64_t *dcol;                ///< Column indices of filter matrix
    double *dval;             ///< Filter weights
    double *Vector;           ///< Input vector to be filtered
    double *VectorFiltered;   ///< Shared output vector (filtered result)
    double *weight_sum;       ///< Shared sum of weights for normalization
    double q;                 ///< Filter exponent (e.g., 1.0 or 2.0)
    int ne;                   ///< Number of elements
    int block_read;           ///< Number of entries read in the current block
} ThreadArgs;

/**
 * @brief Pthread worker that performs atomic accumulation of filter contributions.
 *
 * Each thread handles a subset of entries in the filter triplets and adds their contributions
 * atomically to the shared output and weight arrays.
 *
 * @param args_ptr Pointer to ThreadArgs struct
 * @return NULL
 */
void *thread_filter_worker_atomic(void *args_ptr);

/**
 * @brief Main filtering function that streams filter data from disk and applies the filter using pthreads.
 *
 * Reads filter triplets (row, col, val) in blocks from disk, launches pthreads to process them,
 * accumulates results using atomic updates, and finally normalizes the filtered output vector.
 *
 * @param Vector           Input vector to be filtered
 * @param VectorFiltered   Output vector (must be allocated by caller, will be normalized in-place)
 * @param filternnzElems   Not used (can be NULL)
 * @param ne               Number of elements in the vector
 * @param fnnzassumed      Estimated nonzeros per element (unused)
 * @param q                Filter power/exponent
 * @param filternnz_total  Total number of filter triplets to process
 */
void filterDensity_buffered_mt(double *Vector, double *VectorFiltered,
                              int *filternnzElems,
                              int *ne_ptr, int *fnnzassumed_ptr,
                              double *q_ptr, int filternnz_total);

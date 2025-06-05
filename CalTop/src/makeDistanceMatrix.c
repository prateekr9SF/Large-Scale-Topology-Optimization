/**
 * Constructs a sparse filter matrix based on distances between element centroids.
 *
 * Each element `i` is connected to its neighbors `j` (where `j >= i`) if they lie
 * within a distance `rmin`. The filter matrix entry is weighted by `(rmin - distance)`.
 *
 * This function avoids high memory use by assuming a fixed maximum number of neighbors (`fnnzassumed`)
 * and storing the result in flat arrays: row-wise in C-style.
 *
 * @param ne                Total number of elements
 * @param ttime, time       Not used (included for compatibility)
 * @param ne0               Offset to first element index (not used)
 * @param nea, neb          Element range assigned to this process/thread
 * @param elCentroid        Array of shape [ne][3], holding (x,y,z) coordinates for each element's centroid
 * @param rmin              Filter radius
 * @param filternnz         Output: total number of nonzeros (updated during loop)
 * @param FilterMatrixs     Output: filter weights (size = fnnzassumed * ne)
 * @param rowFilters        Output: row indices corresponding to filtered entries (size = fnnzassumed * ne)
 * @param colFilters        Output: column indices (neighbor elements) (size = fnnzassumed * ne)
 * @param filternnzElems    Output: number of neighbors for each row element (size = ne)
 * @param elarr             Mapping from thread-local ii index to global element index (size ≥ neb)
 * @param fnnzassumed       Maximum allowed nonzeros per row
 */
void mafillsm_filter2(int ne, double ttime, double time,
                      int ne0, int nea, int neb,
                      double (*elCentroid)[3], double rmin, int *filternnz,
                      double *FilterMatrixs, int *rowFilters, int *colFilters,
                      int *filternnzElems, int *elarr, int fnnzassumed)
{
    int ii, i, j, offset;
    double xi, yi, zi, xj, yj, zj;
    double d_squared, rmind;
    int dummy1 = fnnzassumed / 3;

    for (ii = nea; ii <= neb; ++ii)
    {
        i = elarr[ii];  // i is the global element index

        // Centroid of element i
        xi = elCentroid[i][0];
        yi = elCentroid[i][1];
        zi = elCentroid[i][2];

        filternnzElems[i] = 0;  // initialize number of neighbors for element i

        for (j = i; j < ne; ++j)
        {
            // Centroid of element j
            xj = elCentroid[j][0];
            yj = elCentroid[j][1];
            zj = elCentroid[j][2];

            // Squared Euclidean distance
            d_squared = (xj - xi) * (xj - xi)
                      + (yj - yi) * (yj - yi)
                      + (zj - zi) * (zj - zi);

            // Check if within filter radius squared
            rmind = rmin * rmin - d_squared;

            if (rmind >= 0.0)
            {
                int nz = filternnzElems[i];
                if (nz >= dummy1)
                    break;  // prevent memory overflow

                offset = i * fnnzassumed + nz;

                rowFilters[offset]    = i;                  // zero-based index
                colFilters[offset]    = j;                  // zero-based index
                FilterMatrixs[offset] = rmin - sqrt(d_squared);  // linear weight

                filternnzElems[i]++;
                (*filternnz)++;
            }
        }
    }
}

#include <math.h>
#include <stdio.h>

void mafillsm_filter2(int ne, double ttime, double time,
                      int ne0, int nea, int neb,
                      double *elCentroid, double rmin, int *filternnz,
                      double *FilterMatrixs, int *rowFilters, int *colFilters,
                      int *filternnzElems, int *elarr, int fnnzassumed)
{
    int ii, i, j, row_idx, offset;
    double xi, yi, zi, xj, yj, zj;
    double dist_sq, rmind, weight;
    int dummy1 = fnnzassumed / 3;

    for (ii = nea; ii <= neb; ++ii)
    {
        i = elarr[ii];  // 0-based indexing

        xi = elCentroid[3 * i + 0];
        yi = elCentroid[3 * i + 1];
        zi = elCentroid[3 * i + 2];

        filternnzElems[i] = 0;

        for (j = i; j < ne; ++j)
        {
            xj = elCentroid[3 * j + 0];
            yj = elCentroid[3 * j + 1];
            zj = elCentroid[3 * j + 2];

            dist_sq = (xj - xi) * (xj - xi) +
                      (yj - yi) * (yj - yi) +
                      (zj - zi) * (zj - zi);

            rmind = rmin * rmin - dist_sq;

            if (rmind >= 0.0)
            {
                row_idx = filternnzElems[i];
                if (row_idx >= dummy1) break;

                weight = rmin - sqrt(dist_sq);
                offset = row_idx + fnnzassumed * i;

                rowFilters[offset]    = i;
                colFilters[offset]    = j;
                FilterMatrixs[offset] = weight;

                filternnzElems[i]++;
                (*filternnz)++;
            }
        }
    }
}

void mafillsm_filter2(int ne, double ttime, double time,
                      int ne0, int nea, int neb,
                      double *elCentroid, double rmin, int *filternnz,
                      double *FilterMatrixs, int *rowFilters, int *colFilters,
                      int *filternnzElems, int *elarr, int fnnzassumed)
{
    int ii, i, j, offset;
    double xi, yi, zi, xj, yj, zj;
    double d_squared, d, rmind;
    int dummy1 = fnnzassumed / 3;

    for (ii = nea; ii <= neb; ++ii)
    {
        i = elarr[ii];  // already 0-based

        xi = elCentroid[3 * i + 0];
        yi = elCentroid[3 * i + 1];
        zi = elCentroid[3 * i + 2];

        filternnzElems[i] = 0;

        for (j = i; j < ne; ++j)
        {
            xj = elCentroid[3 * j + 0];
            yj = elCentroid[3 * j + 1];
            zj = elCentroid[3 * j + 2];

            d_squared = (xj - xi) * (xj - xi)
                      + (yj - yi) * (yj - yi)
                      + (zj - zi) * (zj - zi);

            rmind = rmin * rmin - d_squared;

            if (rmind >= 0.0)
            {
                int row_index = filternnzElems[i];

                if (row_index >= dummy1)
                    break;

                d = sqrt(d_squared);
                offset = i * fnnzassumed + row_index;

                rowFilters[offset]    = i;
                colFilters[offset]    = j;
                FilterMatrixs[offset] = (rmin - d) * (rmin - d);

                filternnzElems[i]++;
                (*filternnz)++;
            }
        }
    }
}

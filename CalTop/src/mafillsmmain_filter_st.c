#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"



// Main filter routine (single-threaded)
void mafillsmmain_filter2_st(int *ipkon, double *rmin, int *filternnz,
                          int *ne, double *ttime, double *time,
                          int *ne0, double *elCentroid, double *FilterMatrixs,
                          int *rowFilters, int *colFilters,
                          int *filternnzElems, int *fnnzassumed)
{
    // Element array (copy of all indices)
    int *elarr;
    NNEW(elarr, ITG, *ne);
    for (ITG i = 0; i < *ne; ++i)
        elarr[i] = i;

    // Call filter computation over all elements (no threads)
    mafillsm_filter2(*ne, *ttime, *time,
                     *ne0, 0, *ne - 1,  // Full range
                     elCentroid, *rmin, 0, filternnz,
                     FilterMatrixs, rowFilters, colFilters,
                     filternnzElems, elarr, *fnnzassumed);

    SFREE(elarr);
}

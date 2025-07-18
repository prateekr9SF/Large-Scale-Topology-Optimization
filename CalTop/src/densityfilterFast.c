#include <stdio.h>
#include <math.h>
#include "Calculix.h"
#include <unistd.h>



void densityfilterFast(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
                   ITG *ne, double *ttime, double *timepar,
                   ITG *mortar, double *rmin, ITG *filternnz,
                   double *FilterMatrixs, ITG *rowFilters, ITG *colFilters,
                   ITG *filternnzElems, ITG itertop, ITG *fnnzassumed)

    {

        char *lakon = NULL;
        ITG i, j, ne0, *kon = NULL, *ipkon = NULL;
        double dtime, time, *tper;


        ipkon = *ipkonp; lakon = *lakonp; kon = *konp;
        tper = &timepar[1];
        time = *tper; dtime = *tper;
        ne0 = *ne;

        // Initilaize buid filter flag
        int build_filter =0;

        // See if filter files already exit on path.
        if (access("drow.dat", F_OK) != 0 || access("dcol.dat", F_OK) != 0 || access("dval.dat", F_OK) != 0)
        {   
            // Filter files not found, should build filter. 
            build_filter = 1;
        }

        if (build_filter == 1) 
        {
            printf("Filter files not found => building filter matrix\n");

            // Allocate memory for element centroids
            double *eleCentroid = NULL;
            NNEW(elCentroid, double, 3 * ne0);

            printf("Computing element centroids...");
            mafillsmmain_filter(co, nk, kon, ipkon, lakon, ne, ttime, &time, mortar, &ne0, elCentroid);





    }
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
            printf("done!\n");



            printf("Building and writing symmetric distance matrix...");

            FILE *frow = fopen("drow.dat", "w"); // Row indices
            FILE *fcol = fopen("dcol.dat", "w"); // Col indices
            FILE *fval = fopen("dval.dat", "w"); // Kernel weight
            FILE *fdnnz = fopen("dnnz.dat", "w"); // non-zeros

            if (!frow || !fcol || !fval || !fdnnz) 
            {
                fprintf(stderr, "Error: Failed to open output files.\n");
                exit(EXIT_FAILURE);
            }

            double rmin_local = *rmin;
            double sum = 0.0;
            ITG nnz_total = 0;

            // Begin loop over all elements
            for (i = 0; i < ne0; ++i) 
            {
                int count = 0;

                double xi = elCentroid[3 * i + 0];
                double yi = elCentroid[3 * i + 1];
                double zi = elCentroid[3 * i + 2];

                for (j = 0; j < ne0; ++j) 
                {
                    if (i == j) continue;

                    double xj = eleCentroid[3 * j + 0];
                    double yj = eleCentroid[3 * j + 1];
                    double zj = eleCentroid[3 * j + 2];

                    double dx = xi - xj;
                    double dy = yi - yj;
                    double dz = zi - zj;

                    double dist = sqrt(dx * dx + dy * dy + dz * dz);

                    if (dist <= rmin_local)
                    {
                        double w = rmin_local - dist;

                        // write (i, j, w)
                        fprintf(frow, "%d\n", i);
                        fprintf(fcol, "%d\n", j);
                        fprintf(fval, "%.6f\n", w);

                        // Mirror (j, i w)
                        fprintf(frow, "%d\n", j);
                        fprintf(fcol, "%d\n", i);
                        fprintf(fval, "%.6f\n", w);

                        sum += 2 * w;
                        count++;

                    }
                }

                filternnzElems[i] = count;
                fprintf(fdnnz, "%d\n", count);

                if (count > *fnnzassumed) 
                {
                    printf("WARNING: Too many neighbours (%d) for element %d. Increase fnnzassumed or reduce rmin. \n", count , i);
                    exit(EXIT_FAILURE);
                }

                nnz_total += 2 * count;
            }

            // Close all files
            fclose(frow);
            fclose(fcol);
            fclose(fval);
            fclose(fdnnz);

            SFREE(elCentroid);

            *filternnz = nnz_total;


            printf("done!\n");
            printf("Total filter nonzeros (including symmetry): %d\n", nnz_total);
            printf("DVAL SUM FOR FULL MATRIX: %f\n", sum);


        }

    }


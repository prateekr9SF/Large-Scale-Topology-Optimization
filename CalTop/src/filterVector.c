/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

#define BLOCK_SIZE 1000000

void filterVector(ITG **ipkonp,double *Vector, double *VectorFiltered,double *FilterMatrix,
ITG *filternnzElem,ITG *rowFilter,ITG *colFilter,ITG *ne,double *ttime, double *timepar, ITG *fnnzassumed, double *q, int fnnz)
{




  ITG i,j,ne0,*ipkon=NULL;

  


  double *tper;
  
 

  double dtime,time;

  ipkon=*ipkonp;
  tper=&timepar[1];
  time=*tper;
  dtime=*tper;
  ne0=*ne;
   

/*      FILE *elCentroid_file;
       elCentroid_file=fopen("Centroidindensityfilter.dat","w"); //open in write mode
      for(int iii=0;iii<3*ne0;iii++){
        fprintf(elCentroid_file,"%.3f\n",elCentroid[iii]);
        }
        fclose(elCentroid_file);
*/


// Legacy method: calculate Apply filter to a vector
/*  mafillsmmain_Vectorfilter(ipkon,Vector,VectorFiltered,FilterMatrix,filternnzElem,
            rowFilter,colFilter,ne,ttime,&time,&ne0, fnnzassumed,q);
*/


  // Dynamically determine fnnz
 //   int fnnz = count_lines("drow.dat");

// Direct file buffering method (trial)
mafillsmvectorfilter_buffered_filtering(Vector, VectorFiltered,
                                        filternnzElem,  // corrected
                                        *ne, *fnnzassumed,
                                        *q, fnnz);     // corrected q and filternnz


/*      FILE *filter_file;
       filter_file=fopen("Filterfilte.dat","w"); //open in write mode
      for(int iii=0;iii<100*ne0;iii++)
        fprintf(filter_file,"%d , %d , %.3f\n",rowFilter[iii],colFilter[iii],FilterMatrix[iii]);
        }
        fclose(filter_file);

*/




  (*ttime)+=(*tper);

  return;
}



void mafillsmvectorfilter_buffered_filtering(double *Vector, double *VectorFiltered,
                                             int *filternnzElems,
                                             int ne, int fnnzassumed,
                                             double q, int filternnz_total)
{
    FILE *frow = NULL, *fcol = NULL, *fval = NULL;
    int *drow_block = NULL, *dcol_block = NULL;
    double *dval_block = NULL;

    printf("Number of non-zeros beofre reading from disk: %d \n", filternnz_total);
    printf("Streaming filter matrix from disk...\n");

    frow = fopen("drow.dat", "r");
    fcol = fopen("dcol.dat", "r");
    fval = fopen("dval.dat", "r");

    if (!frow || !fcol || !fval) {
        perror("Error opening filter input files");
        exit(EXIT_FAILURE);
    }

    drow_block = (int *)malloc(BLOCK_SIZE * sizeof(int));
    dcol_block = (int *)malloc(BLOCK_SIZE * sizeof(int));
    dval_block = (double *)malloc(BLOCK_SIZE * sizeof(double));

    if (!drow_block || !dcol_block || !dval_block) {
        fprintf(stderr, "Memory allocation failed for block buffers.\n");
        exit(EXIT_FAILURE);
    }

    double *weight_sum = (double *)calloc(ne, sizeof(double));  // stores denominator per row

    int total_read = 0;
    while (total_read < filternnz_total) {
        int remaining = filternnz_total - total_read;
        int block_read = (remaining < BLOCK_SIZE) ? remaining : BLOCK_SIZE;

        // Read one block
        for (int i = 0; i < block_read; ++i) {
            if (fscanf(frow, "%d", &drow_block[i]) != 1 ||
                fscanf(fcol, "%d", &dcol_block[i]) != 1 ||
                fscanf(fval, "%lf", &dval_block[i]) != 1) {
                fprintf(stderr, "Error reading triplet at %d\n", total_read + i);
                exit(EXIT_FAILURE);
            }
        }

        // Apply filter contribution
        for (int i = 0; i < block_read; ++i) {
            int row = drow_block[i] - 1;
            int col = dcol_block[i] - 1;
            double w = pow(dval_block[i], q);

            VectorFiltered[row] += w * Vector[col];
            weight_sum[row] += w;
        }

        total_read += block_read;
    }

    // Normalize filtered values
    for (int i = 0; i < ne; ++i) {
        VectorFiltered[i] = (weight_sum[i] > 0.0) ? VectorFiltered[i] / weight_sum[i] : 0.0;
    }

    // Cleanup
    fclose(frow); fclose(fcol); fclose(fval);
    free(drow_block); free(dcol_block); free(dval_block);
    free(weight_sum);
}


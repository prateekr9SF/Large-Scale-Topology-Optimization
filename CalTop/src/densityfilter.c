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
#include <unistd.h>

void densityfilter(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
	                ITG *ne,
                  double *ttime, double *timepar,
	                ITG *mortar,double *rmin,ITG *filternnz,
                  double *FilterMatrixs,ITG *rowFilters,ITG *colFilters,ITG *filternnzElems, ITG itertop,ITG *fnnzassumed){

  char *lakon=NULL;

  ITG i,j,ne0,*kon=NULL,*ipkon=NULL;

  double *tper;
  

  double dtime,time;

  ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;

  tper=&timepar[1];

  time=*tper;
  dtime=*tper;

  ne0=*ne;
  
  /* Define counter for filter matrix eval */
  int build_filter = 0;

  /* Check if filter files exist */
  if (access("drow.dat", F_OK) != 0 || access("dcol.dat", F_OK) != 0 || access("dval.dat", F_OK) != 0) 
  {
    build_filter = 1;
  }


  
  /* Filter files not found, build the density filter */
  if(build_filter==1)
  {
   
    printf("Filter files not found => building filter matrix \n");

    double *elCentroid=NULL; //pointer to store Centroid of elements
    NNEW(elCentroid,double,3*ne0);  //allocate memory to element CG, initialize to 0 

    /* Calculate centeroid of elements */
    printf("Computing element centeroids...");
    mafillsmmain_filter(co,nk,kon,ipkon,lakon,ne,ttime,&time,mortar,&ne0,elCentroid);


    /* Calculate density filter */
    printf("\nComputing distance matrix...");
    mafillsmmain_filter2(ipkon,rmin,filternnz,ne,ttime,&time,&ne0,elCentroid,
                                    FilterMatrixs,rowFilters,colFilters,filternnzElems,fnnzassumed);
    
    /* Free elCentroid as no longer required */
    SFREE(elCentroid);

    
    /* Pointers for filter matrix  */
    FILE *drow; FILE *dcol; FILE *dnnz; FILE *dval;			

    /* Go through each nnz and copy to respective other half, must be in serial */
    printf("Constructing filter matrix...");
    FORTRAN(mafillsm_expandfilter,(FilterMatrixs,filternnzElems,rowFilters,colFilters,ne,ttime,&time,&ne0,fnnzassumed));
    printf("done! \n");
    printf("Writing row indices to file...");

    /* Write non zero row values for density filter */
    drow=fopen("drow.dat","w"); //open in write mode
    printf("done!\n");
    
    for(int iii=0;iii< (*fnnzassumed)*(ne0);iii++)
    {
      if(FilterMatrixs[iii]>0)
      {
        fprintf(drow,"%d\n",rowFilters[iii]);
      }
    }
    fclose(drow);

    printf("Writing column indices to file...");
    /* Write non zero col values for density filter */
    dcol=fopen("dcol.dat","w"); //open in write mode
    printf("done!\n");

    for(int iii=0;iii<(*fnnzassumed)*ne0;iii++)
    {
      if(FilterMatrixs[iii]>0)
      {
        fprintf(dcol,"%d\n",colFilters[iii]);
      }
    }
    fclose(dcol);

    printf("Writing element values to file...");
    /* Write non zero filter values for density filter */
    dval=fopen("dval.dat","w"); //open in write mode
    printf("done!\n");

    for(int iii=0;iii<*fnnzassumed * ne0;iii++)
    {
      if(FilterMatrixs[iii]>0)
      {
        fprintf(dval,"%.6f\n",FilterMatrixs[iii]);
      }
    }
    fclose(dval);

    printf("Writing non-zero values to file...");
    /* Write number of non zero filter values for each element */
    dnnz=fopen("dnnz.dat","w"); //open in write mode
    printf("done!\n");

    for(int iii=0;iii<ne0;iii++)
    {
      fprintf(dnnz,"%d\n",filternnzElems[iii]);

      if(filternnzElems[iii]>*fnnzassumed)
      {
        printf("WARNING: Number of elements,%d inside filter radius for element %d exceeds %d", filternnzElems[iii],iii,*fnnzassumed);
        printf("CAUTION: ADJUST FILTER SETTING \n");
        exit(1);
      }
    }
    fclose(dnnz);
  }
  else
  {
    printf("Filter files found => assembling filter matrix\n");

    /* Read non zeros in each row from dnnz.dat and calculate total nnzs */
    *filternnz=0; //initialize
    int iii=0;
		FILE *dnnzw;

    dnnzw=fopen("dnnz.dat","r"); //open in read mode

    if (dnnzw!=NULL)
    {
      for (iii=0;iii<ne0;iii++)
      {
        fscanf(dnnzw,"%d",&filternnzElems[iii]);
        *filternnz+=filternnzElems[iii];

        if(filternnzElems[iii]>*fnnzassumed)
        {
          printf("WARNING during read: Number of elements,%d inside filter radius for element %d exceeds %d",filternnzElems[iii],iii,*fnnzassumed);
          printf("CAUTION: ADJUST FILTER SETTING \n");
          exit(1);
        }
      }

      fclose(dnnzw);
    }
    else
    {
      perror("Error reading dval.dat");
      exit(1);
    }

    
    printf("Allocating memory for filter files...");
    double *dval=NULL; //pointer to store density filter values
    NNEW(dval,double,*filternnz);  //allocate memory to dval 

    ITG *drow=NULL; //pointer to store density filter rows 
    NNEW(drow,ITG,*filternnz);  

    ITG *dcol=NULL; //pointer to store density filter cols
    NNEW(dcol,ITG,*filternnz);  

    printf("Done \n");
    FILE *dcolw;

    /*
    dcolw=fopen("dcol.dat","r"); //open in read mode

    if (dcolw!=NULL)
    {
      for (iii=0;iii<*filternnz;iii++)
      {
        fscanf(dcolw,"%d",&dcol[iii]);
      }

      fclose(dcolw);
    }
    else
    {
      perror("Error reading dcol.dat");
    }

    FILE *droww;

    droww=fopen("drow.dat","r"); //open in read mode

    if (droww!=NULL)
    {
      for (iii=0;iii<*filternnz;iii++)
      {
        fscanf(droww,"%d",&drow[iii]);
      }

      fclose(droww);

    }
    else
    {
      perror("Error reading drow.dat");
    }

    FILE *dvalw;

    dvalw=fopen("dval.dat","r"); //open in read mode

    if (dvalw!=NULL)
    {
      for (iii=0;iii<*filternnz;iii++)
      {
        fscanf(dvalw,"%lf",&dval[iii]);
      }

      fclose(dvalw);
    }
    else
    {
      perror("Error reading dval.dat");
    }
    */
    printf("Assembling density filter \n");
    /* Legacy method  */
    //FORTRAN(readfilter,(FilterMatrixs,filternnzElems,rowFilters,colFilters,ne,ttime,&time,&ne0,filternnz,drow,dcol,dval,fnnzassumed));
    
    /* c-based method with in-memory drow, dcol and dval handling */
    //assembleFilter(FilterMatrixs, rowFilters, colFilters,filternnzElems, drow, dcol, dval, ne, ne0, filternnz,fnnzassumed);

    /* c-based method with I/O-based drow, dcol and dval handling */
    //assembleFilter_beta(FilterMatrixs, rowFilters, colFilters,filternnzElems, ne, ne0, filternnz,fnnzassumed); 

    /* c-based method with buffered I/O-based drow, dcol and dval handling */
    //assembleFilter_beta_buffer(FilterMatrixs, rowFilters, colFilters,filternnzElems, ne, ne0, filternnz,fnnzassumed); 

    /* c-based method with I/O-based drow, dcol and dval handling and filter.bin I/O */
    assembleFilter_beta_to_binary("filter.bin",filternnz,fnnzassumed);

    double val_0 = FilterMatrixs[0];
    double val_1 = FilterMatrixs[1];
    double val_2 = FilterMatrixs[2];

    printf("\nFirst value: %f", val_0);
    printf("\nSecond value: %f", val_1);
    printf("\nThird value: %f", val_2);

    /* Density filter build, free up memory */
    SFREE(dcol);
    SFREE(drow);
    SFREE(dval);
  }

  (*ttime)+=(*tper);
  return;
}

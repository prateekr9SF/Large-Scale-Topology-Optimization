#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

/*!
 \file
 \brief Reads or initializes the density (rho) values from a file and stores them in a vector.

 \details
 This function reads rho values from a file named `density.dat` and stores them in the provided vector `design`. 
 If the file does not exist, it initializes all rho values to 1.0 and creates the `density.dat` file with these values.
 The function is designed to work with the finite element method by handling density data for `ne` elements.

 \param[in, out] design A pointer to the vector where rho values will be stored. 
 If the file does not exist, this vector will be initialized to 1.0.
 \param[in] ne The number of elements, which determines the size of the `design` vector.

 \details
 - If `density.dat` exists, it reads the rho values from the file.
 - If `density.dat` does not exist, it initializes the `design` vector to 1.0 for all elements and creates the file.

 Example usage:
 \code
 double *design = (double *)malloc(ne * sizeof(double));
 rho(design, ne);
 \endcode

 \note
 The file `density.dat` is either created or overwritten if it does not exist.

 \author G.Das, P.Ranjan
 \date 2019-2024
*/

void rho(double *design,int ne)
{

    FILE *rhoFile;
    
    /* Open the file in read mode */
    rhoFile=fopen("density.dat","r");
    int i;  //counter

    /* Initialize design array with ones */
    for(i=0; i<ne; i++)
    {
        design[i]=1.0000;
    }

    /* If the density file does not exist, create one */
    if(rhoFile==NULL)
    {
        printf("\n...density.dat not found, initialized to 1");

        rhoFile=fopen("density.dat","w");

        for (i=0;i<ne;i++)
        {
            fprintf(rhoFile,"%.15f \n",design[i]);
        }
    }
    
    /* desinty file exists, read and populate design[] */
    else
    {
        for (i=0;i<ne;i++)
        {
            fscanf(rhoFile,"%lf",&design[i]);
        }
    }
    
    /* all operations done, close the file */
    fclose(rhoFile);
}

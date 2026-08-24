
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


/**
 * Computes the total compliance of the structure
 * @param ne                Number of elements
 * @param elComp            Array of element compliance values
 */

 void getCompliance(int ne,const double *elComp,double *compliance_sum)
{
    /* Loop over all elements and write their sensitivities to file */
    for (int i = 0; i < ne; i++)
    {
        /* compute total element compliance */
        *compliance_sum += elComp[i];
    }
}
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/**
 * Writes volume sensitivities to a CSV file.
 *
 * @param ne                Number of elements.
 * @param eleVol            Array of original (geometric) element volumes.
 * @param rhoPhys           Array of filtered element densities.
 */
void write_objectives(int ne,
                                const double *eleVol,
                                const double *rhoPhys,
                                const double * compliance_sum)
{
    const char *filename = "objectives.csv";

    /* total material volume with rho = 1 */
    double initialVol_sum = 0;

    /* total material volume with optimized rho */
    double  designVol_sum = 0;

    /* design domain discreteness */
    double discreteness = 0.0;

    /* Check for old sensitivity file and delete it */
    if (access(filename, F_OK) == 0)
    {
        if (remove(filename) != 0)
        {
            perror("Error deleting existing objectives file");
            exit(EXIT_FAILURE);
        }
    }

    /* Open new sensitivity file */
    FILE *obj_file = fopen(filename, "w");
    if (!obj_file) {
        perror("Error opening objectives.csv");
        exit(EXIT_FAILURE);
    }

    /* NOTE:
        ORIGINAL VOLUME := Volume of the structure with rho = 1 
        DESIGN VOLUME   := Volume of the structure with filtered rho
    */


    /* Write file header */
    fprintf(obj_file, "COMPLIANCE, ORIGINAL VOLUME, DESIGN VOLUME, VOLUME_FRACTION, DISCRETENESS\n");

    /* Loop over all elements and compute the initial and current volume*/
    for (int i = 0; i < ne; i++)
    {
        initialVol_sum+= eleVol[i];
        designVol_sum+= (eleVol[i]*rhoPhys[i]);
        discreteness += rhoPhys[i] * (1.0 - rhoPhys[i]);

        //fprintf(sens_file, "%.15f,%.15f,%.15f\n", eleVol[i], eleVol[i] * rhoPhys[i], eleVolFiltered[i]);
    }

     /* Compute volume fraction */
    double volume_fraction = designVol_sum / initialVol_sum;

    /* Write structure compliance and volume to file */
    fprintf(obj_file, "%.15f, %.15f, %.15f, %.15f, %.15f \n", *compliance_sum, initialVol_sum, designVol_sum, volume_fraction, discreteness);
    
    fclose(obj_file);
}

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/**
 * Writes volume sensitivities to a CSV file.
 *
 * @param ne                Number of elements.
 * @param eleVol            Array of original (geometric) element volumes.
 * @param rhoPhys           Array of element densities updated by the optimizer.
 * @param eleVolFiltered    Array of element volume sensitivities (filtered).
 */
void write_volume_sensitivities(int ne,
                                const double *eleVol,
                                const double *rhoPhys,
                                const double *eleVolFiltered)
{
    const char *filename = "volume_sens.csv";

    /* Check for old sensitivity file and delete it */
    if (access(filename, F_OK) == 0)
    {
        if (remove(filename) != 0)
        {
            perror("Error deleting existing volume sensitivity file");
            exit(EXIT_FAILURE);
        }
    }

    /* Open new sensitivity file */
    FILE *sens_file = fopen(filename, "w");
    if (!sens_file) {
        perror("Error opening volume_sens.csv");
        exit(EXIT_FAILURE);
    }

    /* Write file header */
    fprintf(sens_file, "ELEMENT VOLUME, CURRENT ELEMENT VOLUME, VOLUME GRADIENT\n");

    /* Loop over all elements and write their sensitivities to file */
    for (int i = 0; i < ne; i++)
    {
        fprintf(sens_file, "%.15f,%.15f,%.15f\n", eleVol[i], eleVol[i] * rhoPhys[i], eleVolFiltered[i]);
    }

    fclose(sens_file);
}

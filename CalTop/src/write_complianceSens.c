
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Function to open a file and handle errors */
FILE *open_file(const char *filename, const char *mode) 
{
    FILE *file = fopen(filename, mode);
    if (!file) {
        perror(filename);
        exit(EXIT_FAILURE);
    }
    return file;
}

/**
 * Writes compliance sensitivities to a csv file and returns the 
 * total compliance of the structure
 * @param ne            Number of elements
 * @param gradCompl     Array of un-filtered compliance sensitivities
 * @param gradComplFiltered Array of filtered compliance sensitivities
 * @param elComp            Array of element compliance values
 */

 void write_compliance_sensitivities(int ne,
                                     const double *gradCompl,
                                     const double *gradComplFiltered,
                                     const double *elComp,
                                     double *compliance_sum)
{
    const char *filename = "compliance_sens.csv";

    /* Check for old sensitivity file and delete it */
    if (access(filename, F_OK) == 0)
    {
        if (remove(filename) != 0)
        {
            perror("Error deleting existing compliance sensitibity file! \n");
            exit(EXIT_FAILURE);
        }
    }

    /* Old file purged, write new sensitivity file */
    FILE *sens_file = open_file(filename, "w");

    /* Write file header */
    fprintf(sens_file, "Compliance GRADIENT UNFILTERED, Compliance GRADIENT FILTERED \n");
    

    /* Loop over all elements and write their sensitivities to file */
    for (int i = 0; i < ne; i++)
    {
        fprintf(sens_file, "%.15f,%.15f\n", gradCompl[i], gradComplFiltered[i]);

        /* compute total element compliance */
        *compliance_sum += elComp[i];
    }

    fclose(sens_file);

}
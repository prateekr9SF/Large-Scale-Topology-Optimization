#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  // For access() function

// Function to open a file and handle errors
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
 * Writes compliance and volume sensitivities to a CSV file and returns
 * compliance sum, initial volume sum, and design volume sum.
 * 
 * @param ne               Number of elements.
 * @param gradCompl        Array of un-filtered compliance gradients.
 * @param gradComplFiltered Array of filtered compliance gradients.
 * @param elCompl          Array of element compliance values.
 * @param eleVol           Array of element volumes.
 * @param rhoPhys          Array of physical density values.
 * @param eleVolFiltered   Array of filtered volume values.
 * @param compliance_sum   Pointer to store the compliance sum.
 * @param initialVol_sum   Pointer to store the initial volume sum.
 * @param designVol_sum    Pointer to store the design volume sum.
 */
void write_elastic_sensitivities(int ne, 
                                 const double *gradCompl, 
                                 const double *gradComplFiltered, 
                                 const double *elCompl, 
                                 const double *eleVol, 
                                 const double *rhoPhys, 
                                 const double *eleVolFiltered,
                                 double *compliance_sum,
                                 double *initialVol_sum,
                                 double *designVol_sum) 
{
    const char *filename = "elastic_sens.csv";

    // Check if the file exists and delete it
    if (access(filename, F_OK) == 0) 
    {
        if (remove(filename) != 0) 
        {
            perror("Error deleting existing CSV file");
            exit(EXIT_FAILURE);
        }
    }

    // Open new CSV file
    FILE *sens_file = open_file(filename, "w");

    // Write CSV header
    fprintf(sens_file, "Element,Compliance,Compliance_Filtered,"
                       "Volume,Volume_Current,Volume_Filtered\n");

    // Initialize summation variables
    *compliance_sum = 0;
    *initialVol_sum = 0;
    *designVol_sum = 0;

    // Loop through all elements
    for (int i = 0; i < ne; i++) 
    {
        *compliance_sum += elCompl[i];
        *initialVol_sum += eleVol[i];
        *designVol_sum += eleVol[i] * rhoPhys[i];

        // Write element-wise sensitivities
        fprintf(sens_file, "%d,%.15f,%.15f,%.15f,%.15f,%.15f\n",
                i, gradCompl[i], gradComplFiltered[i], 
                eleVol[i], eleVol[i] * rhoPhys[i], eleVolFiltered[i]);
    }
    fclose(sens_file);
}

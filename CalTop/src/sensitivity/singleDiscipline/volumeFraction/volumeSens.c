#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

/**
 * Evalates element volume fraction sensitivities
 *
 * @param ne                Number of elements.
 * @param eleVol            Array of original (geometric) element volumes.
 * @param rhoPhys           Array of filtered element densities.
 */

 void volumeSens(int ne, const double *eleVol, const int *passiveIDs,
                int numPassive, double *volFracSens)
{

    /* Allocate element volume sensitivity array */
    //double *volFracSens = calloc(ne, sizeof(double));

    //if (!volFracSens)
   // {
   //     perror("calloc failed for volFracSens");
   //     exit (EXIT_FAILURE);
   // }

    /* Pre-allocate skin element indices here */
    int *is_skin = calloc(ne, sizeof(int));

    if (!is_skin)
    {
        perror("calloc failed for is_skin");
        exit(EXIT_FAILURE);
    }

    for (int k = 0; k <numPassive; k++)
    {
        int id = passiveIDs[k] -1;
        if (id >= 0 && id < ne)
        {
            is_skin[id] =1;
        }
    }

    int nthreads = omp_get_max_threads();

    printf("Using %d CPU(s) for volume fraction sensitivity calculation", nthreads);
    fflush(stdout);

    /* Evaluate full solid volume (note: should exclude passive volume if applicable) */
    
    double fullSolidVol_sum = 0.0;

    #pragma omp parallel for reduction(+:fullSolidVol_sum)
    for (int i = 0; i <ne; i++)
    {
        if (!is_skin[i])
        {
            fullSolidVol_sum += eleVol[i];
        }
    }

    /* sanity check for full solud volume */
    if (fullSolidVol_sum <= 0.0)
    {
        fprintf(stderr, "Error: Full solid volume is zero or negative. Check element volumes.\n");
        free(is_skin);
        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for
    for (int i = 0; i < ne; i++)
    {
        if (is_skin[i])
        {
            volFracSens[i] = 0.0;  // Passive elements have zero sensitivity
        }
        else
        {
            volFracSens[i] = eleVol[i] / fullSolidVol_sum;
        }
    }

    free(is_skin);
}


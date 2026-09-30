#include <stdio.h>
#include <math.h>

#include "CalFSI.h"

void write_converged_rhs(const char *filename,
                                const double *b,
                                ITG neq)
{
    FILE *fp = fopen(filename, "w");

    if (fp == NULL)
    {
        fprintf(stderr,
                "ERROR: Could not open %s for writing.\n",
                filename);
        perror("fopen");
        FORTRAN(stop,());
    }

    double rhs_norm_sq = 0.0;

    for (ITG i = 0; i < neq; ++i)
    {
        /*
         * Fixed-point notation:
         * one equation-space RHS value per line.
         */
        if (fprintf(fp, "%.10f\n", b[i]) < 0)
        {
            fprintf(stderr,
                    "ERROR: Failed while writing %s at equation %d.\n",
                    filename,
                    (int)(i + 1));

            fclose(fp);
            FORTRAN(stop,());
        }

        rhs_norm_sq += b[i] * b[i];
    }

    if (fclose(fp) != 0)
    {
        fprintf(stderr,
                "ERROR: Could not properly close %s.\n",
                filename);
        FORTRAN(stop,());
    }

    printf("Structural RHS written to %s\n", filename);
    printf("Number of equations : %d\n", (int)neq);
    printf("RHS L2 norm         : %.10f\n", sqrt(rhs_norm_sq));
    fflush(stdout);
}
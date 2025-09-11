#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "CalculiX.h"

/*--------------------------------------------------------------
  write_cg_sens_csv

  Purpose:
    Write CG sensitivity arrays to a single CSV file with columns:
      dCGx, dCGy, dCGz

  Parameters:
    path  - output CSV file path
    ne    - number of elements (length of each array)
    dCGx  - [ne] sensitivities of CGx w.r.t. rho
    dCGy  - [ne] sensitivities of CGy w.r.t. rho
    dCGz  - [ne] sensitivities of CGz w.r.t. rho

  Returns:
    0 on success, nonzero on error (also prints a message to stderr).

  CSV format:
    dCGx,dCGy,dCGz
    val_x0,val_y0,val_z0
    val_x1,val_y1,val_z1
    ...
----------------------------------------------------------------*/
int write_cg_sens(const char *path,
                      size_t ne,
                      const double *dCGx,
                      const double *dCGy,
                      const double *dCGz)
{
    /* Basic argument checks */
    if (!path || !dCGx || !dCGy || !dCGz) {
        fprintf(stderr, "cg_sens.csv: null argument\n");
        return 1;
    }

    /* Open file for writing (overwrite if exists) */
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cg_sens.csv: cannot open '%s': %s\n",
                path, strerror(errno));
        return 2;
    }

    /* Write header */
    if (fprintf(fp, "CGx GRADIENT,CGy GRADIENT,CGz GRADIENT\n") < 0) {
        fprintf(stderr, "cg_sens.csv: write header failed\n");
        fclose(fp);
        return 3;
    }

    /* Write rows (use high precision; CSV uses '.' as decimal separator) */
    for (size_t i = 0; i < ne; ++i) {
        if (fprintf(fp, "%.15g,%.15g,%.15g\n",
                    dCGx[i], dCGy[i], dCGz[i]) < 0) {
            fprintf(stderr, "write_cg_sens_csv: write row %zu failed\n", i);
            fclose(fp);
            return 4;
        }
    }

    /* Flush and close */
    if (fclose(fp) != 0) {
        fprintf(stderr, "write_cg_sens_csv: close failed for '%s': %s\n",
                path, strerror(errno));
        return 5;
    }

    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>


#include "CalculiX.h"



void addPassiveComponent(const char *filename,
                         const char *description,
                         int **passiveIDs,
                         int *numPassive,
                         int **componentIDs,
                         int *numComponent)
{
    struct stat buffer;

    /* --------------------------------------------------------
       Check whether component element list exists
       -------------------------------------------------------- */

    if (stat(filename, &buffer) != 0)
    {
        printf("File '%s' not found -> "
               "%s will not be defined.\n",
               filename, description);

        *componentIDs = NULL;
        *numComponent = 0;

        return;
    }

    printf("File '%s' found -> "
           "Setting %s definition.\n",
           filename, description);

    fflush(stdout);

    /* --------------------------------------------------------
       Read component element IDs

       This array is retained so that it can be used later
       for component-specific VTU output.
       -------------------------------------------------------- */

    *componentIDs = passiveElements(filename,
                                    numComponent);

    if (*componentIDs == NULL)
    {
        fprintf(stderr,
                "ERROR: Could not read %s\n",
                filename);

        FORTRAN(stop,());
    }

    printf("Read %d %s elements.\n",
           *numComponent, description);

    /* --------------------------------------------------------
       Expand global passive array for worst case
       -------------------------------------------------------- */

    int *tmp = realloc(*passiveIDs,
                       (*numPassive + *numComponent)
                       * sizeof(int));

    if (tmp == NULL)
    {
        fprintf(stderr,
                "ERROR: Could not reallocate passiveIDs "
                "while adding %s.\n",
                description);

        FORTRAN(stop,());
    }

    *passiveIDs = tmp;

    /* --------------------------------------------------------
       Add only unique IDs to global passive domain
       -------------------------------------------------------- */

    int numAdded = 0;
    int numDuplicate = 0;

    for (int i = 0; i < *numComponent; i++)
    {
        int id = (*componentIDs)[i];
        int alreadyPassive = 0;

        for (int j = 0; j < *numPassive; j++)
        {
            if ((*passiveIDs)[j] == id)
            {
                alreadyPassive = 1;
                break;
            }
        }

        if (!alreadyPassive)
        {
            (*passiveIDs)[*numPassive] = id;
            (*numPassive)++;

            numAdded++;
        }
        else
        {
            numDuplicate++;
        }
    }

    /* --------------------------------------------------------
       Summary
       -------------------------------------------------------- */

    printf("%s elements added to passive : %d\n",
           description, numAdded);

    printf("Already passive / duplicates  : %d\n",
           numDuplicate);

    printf("Total passive elements        : %d\n",
           *numPassive);

    fflush(stdout);
}
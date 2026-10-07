#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>


#include "CalculiX.h"




void addPassiveElements(const char *filename,
                        const char *description,
                        int **passiveIDs,
                        int *numPassive)
{
    struct stat buffer;

    int *newPassiveIDs = NULL;
    int numNewPassive = 0;

    /* --------------------------------------------------------
       Check whether element list exists
       -------------------------------------------------------- */

    if (stat(filename, &buffer) != 0)
    {
        printf("File '%s' not found -> "
               "%s will not be defined.\n",
               filename, description);

        return;
    }

    printf("File '%s' found -> "
           "Setting %s definition.\n",
           filename, description);

    fflush(stdout);

    /* --------------------------------------------------------
       Read passive element IDs
       -------------------------------------------------------- */

    newPassiveIDs = passiveElements(filename,
                                    &numNewPassive);

    if (newPassiveIDs == NULL)
    {
        fprintf(stderr,
                "ERROR: Could not read %s\n",
                filename);

        FORTRAN(stop,());
    }

    printf("Read %d %s elements.\n",
           numNewPassive, description);

    /* --------------------------------------------------------
       Allocate enough space for worst case:
       no overlap with existing passive elements
       -------------------------------------------------------- */

    int *tmp = realloc(*passiveIDs,
                       (*numPassive + numNewPassive)
                       * sizeof(int));

    if (tmp == NULL)
    {
        fprintf(stderr,
                "ERROR: Could not reallocate passiveIDs "
                "while adding %s.\n",
                description);

        free(newPassiveIDs);
        free(*passiveIDs);

        FORTRAN(stop,());
    }

    *passiveIDs = tmp;

    /* --------------------------------------------------------
       Add only IDs that are not already passive
       -------------------------------------------------------- */

    int numAdded = 0;
    int numDuplicate = 0;

    for (int i = 0; i < numNewPassive; i++)
    {
        int id = newPassiveIDs[i];
        int alreadyPassive = 0;

        /*
         * Search existing passive IDs.
         *
         * numPassive increases when a new ID is added,
         * so duplicates within the new list are also detected.
         */
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

    printf("%s elements added           : %d\n",
           description, numAdded);

    printf("Already passive / duplicates   : %d\n",
           numDuplicate);

    printf("Total passive elements         : %d\n",
           *numPassive);

    fflush(stdout);

    free(newPassiveIDs);
}
#include <stdio.h>
#include <stdlib.h>

#define BLOCK_SIZE 1000  // Tune this for performance vs memory

/**
 * @brief Populate the flattened filter matrix and associated row/column index arrays
 *        from sparse data given in coordinate format (drow, dcol, dval).
 *
 * This function assumes that the 2D matrices FilterMatrixs, rowFilters, and colFilters
 * are stored in a flattened 1D array using column-major storage, i.e., (i, j) maps to
 * i + fnnzassumed * j.
 *
 * @param FilterMatrixs  Flattened 2D array of filter weights [fnnzassumed x ne]
 * @param rowFilters     Flattened 2D array of row indices [fnnzassumed x ne]
 * @param colFilters     Flattened 2D array of column indices [fnnzassumed x ne]
 * @param filternnzElems Array storing the number of nonzeros expected per row (size ne)
 * @param drow           Row indices of nonzero entries in COO format (size filternnz)
 * @param dcol           Column indices of nonzero entries in COO format (size filternnz)
 * @param dval           Values of nonzero entries in COO format (size filternnz)
 * @param ne             Number of elements (rows in the filter matrix)
 * @param ttime          Total time (not used here, kept for compatibility)
 * @param time           Current time (not used here, kept for compatibility)
 * @param ne0            Starting element index (not used)
 * @param filternnz      Total number of nonzeros in the filter
 * @param fnnzassumed    Assumed max number of nonzeros per row (used to size flattened arrays)
 */
void assembleFilter(double *FilterMatrixs, int *rowFilters, int *colFilters,
                int *filternnzElems, int *drow, int *dcol, double *dval,
                int ne, int ne0,
                int *filternnz, int *fnnzassumed) 
{

    int i;              // Loop counter for nonzeros
    int rowval, colval; // Row and column indices for the current nonzero
    double value;       // Value of the current nonzero

    int index = 1;      // Shared counter for writing into flattened arrays
                        // (shared across all rows — may cause overwrite)
    
    printf("Starting loop \n");
    for (i = 0; i < *filternnz; i++) 
    {
        rowval = drow[i];    // Target row (0-based index)
        colval = dcol[i];    // Target column
        value  = dval[i];    // Corresponding value

        // Compute flattened array offset for (index, rowval)
        //int offset = index + (*fnnzassumed) * rowval;

        int offset = (index - 1) + (*fnnzassumed) * (rowval - 1);
        //FilterMatrixs[offset] = value;

        // Store data in the filter structure
        rowFilters[offset]    = rowval;
        colFilters[offset]    = colval;
        FilterMatrixs[offset] = value;

        // Update index (shared for all rows — may lead to memory overwrites)
        if (index < filternnzElems[rowval]) 
        {
            index++;
        }
        else 
        {
            index = 1;
        }
    }

    return;
}



void assembleFilter_beta(double *FilterMatrixs, int *rowFilters, int *colFilters,
                    int *filternnzElems,
                    int ne, int ne0,
                    int *filternnz, int *fnnzassumed)
{
    FILE *frow = NULL, *fcol = NULL, *fval = NULL;
    int i, rowval, colval, index = 1;
    double value;

    // Open the input files
    frow = fopen("drow.dat", "r");
    fcol = fopen("dcol.dat", "r");
    fval = fopen("dval.dat", "r");

    if (!frow || !fcol || !fval) {
        perror("Error opening one or more input files");
        exit(EXIT_FAILURE);
    }

    printf("Streaming filter data from disk...\n");

    for ( i = 0; i < *filternnz; ++i)
    {
        if (fscanf(frow, "%d", &rowval) != 1 ||
            fscanf(fcol, "%d", &colval) != 1 ||
            fscanf(fval, "%lf", &value) != 1)
            {
                fprintf(stderr, "Error reading line %d in filter data\n", i);
                exit(EXIT_FAILURE);
            }
    
        // Compute flattened index (column-major: index + fnnzassumed * row)
        int offset = (index - 1) + (*fnnzassumed) * (rowval - 1);

        rowFilters[offset]    = rowval;
        colFilters[offset]    = colval;
        FilterMatrixs[offset] = value;

        // Update index (shared for all rows — may lead to memory overwrites)
        if (index < filternnzElems[rowval])
        {
            index++;
        }
        else
        {
            index = 1;
        }
    }

    /* Free memeory */
    fclose(frow);
    fclose(fcol);
    fclose(fval);
}

void assembleFilter_beta_buffer(double *FilterMatrixs, int *rowFilters, int *colFilters,
                    int *filternnzElems,
                    int ne, int ne0,
                    int *filternnz, int *fnnzassumed)
{
    FILE *frow = NULL, *fcol = NULL, *fval = NULL;
    int i, block_read, index = 1;
    int *drow_block = NULL, *dcol_block = NULL;
    double *dval_block = NULL;

    frow = fopen("drow.dat", "r");
    fcol = fopen("dcol.dat", "r");
    fval = fopen("dval.dat", "r");

    if (!frow || !fcol || !fval) 
    {
        perror("Error opening filter input files");
        exit(EXIT_FAILURE);
    }

    // Allocate buffer for block reads
    drow_block = (int *)malloc(BLOCK_SIZE * sizeof(int));
    dcol_block = (int *)malloc(BLOCK_SIZE * sizeof(int));
    dval_block = (double *)malloc(BLOCK_SIZE * sizeof(double));

    if (!drow_block || !dcol_block || !dval_block) 
    {
        fprintf(stderr, "Memory allocation failed for block buffers.\n");
        exit(EXIT_FAILURE);
    }

    printf("Reading filter triplets in blocks of %d...\n", BLOCK_SIZE);

    int total_read = 0;
    while (total_read < *filternnz) 
    {
        int remaining = *filternnz - total_read;
        block_read = (remaining < BLOCK_SIZE) ? remaining : BLOCK_SIZE;

        // Read one block into buffers
        for (i = 0; i < block_read; ++i) 
        {
            if (fscanf(frow, "%d", &drow_block[i]) != 1 ||
                fscanf(fcol, "%d", &dcol_block[i]) != 1 ||
                fscanf(fval, "%lf", &dval_block[i]) != 1)
            {
                fprintf(stderr, "Error reading triplet at %d\n", total_read + i);
                exit(EXIT_FAILURE);
            }
        }

        // Process the block
        for (i = 0; i < block_read; ++i) 
        {
            int rowval = drow_block[i];
            int colval = dcol_block[i];
            double value = dval_block[i];

            int offset = (index - 1) + (*fnnzassumed) * (rowval - 1);

            rowFilters[offset]    = rowval;
            colFilters[offset]    = colval;
            FilterMatrixs[offset] = value;

            if (index < filternnzElems[rowval]) 
            {
                index++;
            } 
            else 
            {
                index = 1;
            }
        }

        total_read += block_read;
    }

    // Cleanup
    fclose(frow);
    fclose(fcol);
    fclose(fval);
    free(drow_block);
    free(dcol_block);
    free(dval_block);
}


void assembleFilter_beta_to_binary(const char* outfile,
                                   int* filternnz,
                                   int *fnnzassumed)
{
              
    FILE *frow = fopen("drow.dat", "r");
    FILE *fcol = fopen("dcol.dat", "r");
    FILE *fval = fopen("dval.dat", "r");
    FILE *fbin = fopen(outfile, "wb");

    if (!frow || !fcol || !fval || !fbin) {
        perror("Error opening one or more input/output files");
        exit(EXIT_FAILURE);
    }

    printf("Streaming filter data from disk into binary...\n");

    int rowval, colval;
    double value;
    for (int i = 0; i < *filternnz; ++i) {
        if (fscanf(frow, "%d", &rowval) != 1 ||
            fscanf(fcol, "%d", &colval) != 1 ||
            fscanf(fval, "%lf", &value) != 1)
        {
            fprintf(stderr, "Error reading line %d in filter data\n", i);
            exit(EXIT_FAILURE);
        }

        int row0 = rowval - 1;  // Convert to 0-based
        int col0 = colval - 1;

        fwrite(&row0, sizeof(int), 1, fbin);
        fwrite(&col0, sizeof(int), 1, fbin);
        fwrite(&value, sizeof(double), 1, fbin);
    }

    fclose(frow);
    fclose(fcol);
    fclose(fval);
    fclose(fbin);
}

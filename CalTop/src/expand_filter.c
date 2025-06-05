#include <stdio.h>

void mafillsm_expandfilter(double *FilterMatrixs, int *filternnzElems,
                           int *rowFilters, int *colFilters,
                           int ne, int fnnzassumed)
{
    for (int i = 0; i < ne; ++i) {
        int row_nnz = filternnzElems[i];

        if (row_nnz > fnnzassumed) {
            printf("[ERROR] filternnzElems[%d] = %d exceeds fnnzassumed = %d\n", i, row_nnz, fnnzassumed);
            continue;
        }

        for (int j = 0; j < row_nnz; ++j) {
            int offset = j + fnnzassumed * i;

            // Defensive checks
            if (offset >= fnnzassumed * ne) {
                printf("[ERROR] offset %d out-of-bounds (max %d)\n", offset, fnnzassumed * ne - 1);
                continue;
            }

            int rowval = rowFilters[offset];
            int colval = colFilters[offset];
            double value = FilterMatrixs[offset];

            if (colval < 0 || colval >= ne || rowval < 0 || rowval >= ne) {
                printf("[ERROR] Invalid colval=%d or rowval=%d at i=%d, j=%d\n", colval, rowval, i, j);
                continue;
            }

            if (colval > rowval) {
                int dest_row = colval;
                int dest_col = rowval;

                int dest_nnz = filternnzElems[dest_row];
                if (dest_nnz >= fnnzassumed) {
                    printf("[WARNING] Overflow at row %d (nnz = %d)\n", dest_row, dest_nnz);
                    continue;
                }

                int new_offset = dest_nnz + fnnzassumed * dest_row;
                if (new_offset >= fnnzassumed * ne) {
                    printf("[ERROR] new_offset %d out-of-bounds\n", new_offset);
                    continue;
                }

                rowFilters[new_offset]    = dest_row;
                colFilters[new_offset]    = dest_col;
                FilterMatrixs[new_offset] = value;

                filternnzElems[dest_row]++;
            }
        }
    }
}

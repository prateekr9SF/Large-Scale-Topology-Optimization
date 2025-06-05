void mafillsm_expandfilter(double *FilterMatrixs, int *filternnzElems,
                           int *rowFilters, int *colFilters,
                           int ne, int fnnzassumed)
{
    for (int i = 0; i < ne; i++) {
        int row_nnz = filternnzElems[i];

        for (int j = 0; j < row_nnz; j++) {
            int offset = j + fnnzassumed * i;

            int rowval = rowFilters[offset];
            int colval = colFilters[offset];
            double value = FilterMatrixs[offset];

            // Expand only if colval > i (upper triangle → reflect to lower)
            if (colval > i) {

                printf("In the upper traingular matrix \n");
                int dest_row = colval;
                int dest_col = i;
                int dest_nnz = filternnzElems[dest_row];

                if (dest_nnz >= fnnzassumed) {
                    printf("[WARNING] Overflow at row %d in expansion\n", dest_row);
                    continue;
                }

                int new_offset = dest_nnz + fnnzassumed * dest_row;

                rowFilters[new_offset] = dest_row;
                colFilters[new_offset] = dest_col;
                FilterMatrixs[new_offset] = value;

                filternnzElems[dest_row]++;
            }
        }
    }
}

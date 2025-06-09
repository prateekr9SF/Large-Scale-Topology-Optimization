#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Counts the number of lines in a given text file.
 *
 * This function opens the specified file in read mode and counts the number
 * of newline characters (`'\n'`) to determine how many lines it contains.
 * It is commonly used to determine the number of entries in a text-based dataset.
 *
 * @param filename Path to the input file (e.g., "drow.dat")
 * @return int Number of lines in the file
 *
 * @note If the file cannot be opened, the function prints an error and exits.
 */
int count_lines(const char *filename) 
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("[count_lines] Cannot open file");
        exit(EXIT_FAILURE);
    }

    int count = 0;
    char ch;
    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\n') count++;
    }

    fclose(fp);
    return count;
}

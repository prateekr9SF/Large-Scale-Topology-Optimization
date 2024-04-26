/**
 * @file main.c
 * @brief Provides memory reallocation utility for CalculiX programs.
 *
 * This file includes a custom memory reallocation function that handles memory resizing
 * with additional error handling and optional debugging and logging features.
 * It is part of the CalculiX, a 3-dimensional finite element program.
 * 
 * CalculiX is licensed under the GNU General Public License version 2.
 */

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
extern int log_realloc; ///< External global variable to control memory allocation logging.

/**
 * @brief Reallocates memory and includes detailed logging and error handling.
 *
 * This function attempts to reallocate memory block to a new size. It provides
 * extensive error reporting and can log successful operations based on environment settings.
 *
 * @param ptr Pointer to the memory previously allocated (may be NULL).
 * @param size New size in bytes to which the memory block is to be resized.
 * @param file Filename from which the reallocation request is made (for logging).
 * @param line Line number in the source file where the reallocation request occurs (for logging).
 * @param ptr_name Name of the variable being reallocated (for logging).
 * @return Pointer to the reallocated memory block. Exits program if reallocation fails and conditions apply.
 */
void *u_realloc(void* ptr, size_t size, const char *file, const int line, const char* ptr_name) {
    void *a; // Pointer to hold the new memory block address.
    char *env; // Variable to hold the environment variable value.

    // Attempt to reallocate the memory block to the new size.
    a = realloc(ptr, size);

    // Check if reallocation failed but the original pointer was not NULL and size was not zero.
    if (a == NULL && ptr != NULL && size != 0) {
        // Print an error message with details about the reallocation attempt.
        printf("*ERROR in u_realloc: error allocating memory\n");
        printf("variable=%s, file=%s, line=%d, size(bytes)=%ld, oldaddress=%ld\n",
               ptr_name, file, line, size, (long int)ptr);
        exit(16); // Exit the program with an error code.
    } else {
        // Check if logging for memory reallocations has been initialized.
        if (log_realloc == -1) {
            log_realloc = 0; // Set flag indicating logging status is checked.
            env = getenv("CCX_LOG_ALLOC"); // Get the environment variable for logging.
            if (env) {
                log_realloc = atoi(env); // Convert environment variable to integer for logging flag.
            }
        }
        
        // If logging is enabled, print reallocation details.
        if (log_realloc == 1) {
            printf("REALLOCATION of variable %s, file %s, line=%d: size(bytes)=%ld, oldaddress= %ld, address= %ld\n",
                   ptr_name, file, line, size, (long int)ptr, (long int)a);
        }
        
        return a; // Return the new memory block pointer.
    }
}

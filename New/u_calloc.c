#include <stdio.h>
#include <stdlib.h>
#include "main.h"

int log_realloc = -1; // Initialize log_realloc to -1 to signal that logging is not set up yet.

/**
 * A customized calloc function to allocate memory with additional error handling and logging.
 *
 * @param num Number of elements to allocate.
 * @param size Size of each element.
 * @param file Name of the source file where the allocation request is made.
 * @param line Line number in the source file of the allocation request.
 * @param ptr_name The name of the variable that will receive the allocated memory.
 * @return Pointer to the allocated block of memory. If the allocation fails, the program will exit with an error code.
 */
void *u_calloc(size_t num, size_t size, const char *file, const int line, const char* ptr_name) {
    void *a; // Pointer to hold the address of the allocated memory.
    char *env; // Variable to hold the environment value for logging.

    // Check if the requested number of elements is zero.
    if (num == 0) {
        a = NULL; // If no elements are requested, return NULL immediately.
        return a;
    }

    a = calloc(num, size); // Attempt to allocate memory for num elements of size bytes each.
    if (a == NULL) {
        // If memory allocation fails, print an error message with details about the request.
        printf("*ERROR in u_calloc: error allocating memory\n");
        printf("variable=%s, file=%s, line=%d, num=%ld, size=%ld\n", ptr_name, file, line, num, size);

        // Additional check for negative numbers which are not valid for memory allocation sizes.
        if (num < 0) {
            printf("\n It looks like you may need the i8 (integer*8) version of CalculiX\n");
        }
        exit(16); // Exit the program with an error code of 16.
    } else {
        // Check if logging of memory allocations is enabled.
        if (log_realloc == -1) {
            log_realloc = 0; // Set to 0 indicating logging checked but not enabled.
            env = getenv("CCX_LOG_ALLOC"); // Retrieve the environment variable for log control.
            if (env) {
                log_realloc = atoi(env); // Convert the environment variable to an integer for logging control.
            }
        }
        
        // If logging is enabled, print details of the memory allocation.
        if (log_realloc == 1) {
            printf("ALLOCATION of variable %s, file %s, line=%d, num=%ld, size=%ld, address= %ld\n", ptr_name, file, line, num, size, (long int)a);
        }
        
        return a; // Return the pointer to the allocated memory.
    }
}

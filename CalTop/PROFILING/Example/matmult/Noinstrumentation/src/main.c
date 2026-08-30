#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"


int main(void)
{
    double *A;
    double *B;
    double *C;

    double norm;

    printf("Allocating memory ...");
    fflush(stdout);
    /* Allocate matrices */
    allocateMatrices(&A, &B, &C);
    printf("done!\n");
    fflush(stdout);

    printf("Populating matrices...");
    fflush(stdout);
    /* Populate matrices */
    populateMatrices(A, B);
    printf("done!\n");

    printf("Multiplying matrices...");
    fflush(stdout);
    /* Matrix multiplication */
    matrixMultiply(A, B, C);
    printf("done!\n");

    printf("Computing norm...");
    fflush(stdout);
    /* Compute matrix norm */
    norm = matrixNorm(C);
    printf("done!\n");
    fflush(stdout);

    printf("Matrix size: %d x %d\n", N, N);
    printf("Frobenius norm of C = %.12e\n", norm);

    /* Free memory */
    free(A);
    free(B);
    free(C);

    return 0;
}
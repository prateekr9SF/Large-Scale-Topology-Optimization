#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
//include <TAU.h>



int main(void)
{
    double *A;
    double *B;
    double *C;
    double norm;

    //TAU_PROFILE("main()", "", TAU_DEFAULT);

    /* 
    Define TAU timers
    *
    * The first argument is the timer handle.
    * The second argument is the name that will appear in TAU.
    * The third argument is the type string.
    * The fourth argument specifies the timer group.
    */
    
    printf("Allocating memory ...");
    fflush(stdout);
    allocateMatrices(&A, &B, &C);
    printf("done!\n");
    fflush(stdout);

    printf("Populating matrices...");
    fflush(stdout);
    populateMatrices(A, B);
    printf("done!\n");


    printf("Multiplying matrices...");
    fflush(stdout);
    
    matrixMultiply(A, B, C);

    printf("done!\n");
    
    printf("Computing norm...");
    fflush(stdout);
    norm = matrixNorm(C);
    printf("done!\n");
    fflush(stdout);

    printf("Matrix size: %d x %d\n", N, N);
    printf("Frobenius norm of C = %.12e\n", norm);


    printf("Adding matrices...");
    fflush(stdout);
    matrixAdd(A, B, C);
    printf("done!\n");
    fflush(stdout);

    printf("Trace calculation...");
    fflush(stdout);
    matrixTrace(A);

    printf("done!\n");
    fflush(stdout);



    /* Free memory */
    free(A);
    free(B);
    free(C);

    return 0;
}
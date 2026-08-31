#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include <TAU.h>



int main(void)
{
    double *A;
    double *B;
    double *C;
    double norm;

    TAU_PROFILE("main()", "", TAU_DEFAULT);

    /* 
    Define TAU timers
    *
    * The first argument is the timer handle.
    * The second argument is the name that will appear in TAU.
    * The third argument is the type string.
    * The fourth argument specifies the timer group.
    */
    
    TAU_PROFILE_TIMER(timer_allocate, "allocateMatrices", "", TAU_USER);
    TAU_PROFILE_TIMER(timer_multiply, "matrixMultiply", "", TAU_USER);
    TAU_PROFILE_TIMER(timer_populate, "populateMatrices", "", TAU_USER);
    TAU_PROFILE_TIMER(timer_norm, "matrixNorm", "", TAU_USER);      

    printf("Allocating memory ...");
    fflush(stdout);

    TAU_PROFILE_START(timer_allocate);
    allocateMatrices(&A, &B, &C);
    TAU_PROFILE_STOP(timer_allocate);
    printf("done!\n");
    fflush(stdout);

    printf("Populating matrices...");
    fflush(stdout);
    
    TAU_PROFILE_START(timer_populate);
    populateMatrices(A, B);
    TAU_PROFILE_STOP(timer_populate);
    printf("done!\n");


    printf("Multiplying matrices...");
    fflush(stdout);
    
    TAU_PROFILE_START(timer_multiply);
    matrixMultiply(A, B, C);
    TAU_PROFILE_STOP(timer_multiply);

    printf("done!\n");
    
    printf("Computing norm...");
    fflush(stdout);

    TAU_PROFILE_START(timer_norm);
    norm = matrixNorm(C);
    TAU_PROFILE_STOP(timer_norm);

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
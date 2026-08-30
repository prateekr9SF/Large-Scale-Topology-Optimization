#ifndef MATRIX_H
#define MATRIX_H

#define N 2000

/* Allocate matrices A, B, and C */
void allocateMatrices(double **A, double **B, double **C);

/* Populate matrices A and B */
void populateMatrices(double *A, double *B);

/* Compute C = A * B */
void matrixMultiply(double *A, double *B, double *C);

/* Compute Frobenius norm of C */
double matrixNorm(double *C);

#endif
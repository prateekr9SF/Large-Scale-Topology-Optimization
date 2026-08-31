#ifndef MATRIX_H
#define MATRIX_H

#define N 2000

/* Allocate matrices A, B, and C */
void allocateMatrices(double **A, double **B, double **C);

/* Populate matrices A and B */
void populateMatrices(double *A, double *B);

/* Compute C = A * B */
void matrixMultiplyOMP(double *A, double *B, double *C);

/* Compute Frobenius norm of C */
double matrixNorm(double *C);

void matrixMultiplyTransposed(double *A, double *BT, double *C);

double matrixFrobeniusProduct(double *A, double *B);

double matrixTrace(double *A);

void matrixVectorMultiply(double *A, double *x, double *y);

void matrixTranspose(double *A, double *AT);

void matrixScale(double *A, double *B, double alpha);

void matrixAdd(double *A, double *B, double *C);

void matrixMultiply(double *A, double *B, double *C);

#endif
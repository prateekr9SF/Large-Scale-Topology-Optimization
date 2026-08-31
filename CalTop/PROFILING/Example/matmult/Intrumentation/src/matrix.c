#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "matrix.h"


/* ------------------------------------------------------------
   Allocate memory for matrices A, B, and C
   ------------------------------------------------------------ */
void allocateMatrices(double **A, double **B, double **C)
{
    *A = (double *)malloc(N * N * sizeof(double));
    *B = (double *)malloc(N * N * sizeof(double));
    *C = (double *)malloc(N * N * sizeof(double));

    if (*A == NULL || *B == NULL || *C == NULL)
    {
        printf("Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }
}


/* ------------------------------------------------------------
   Populate matrices A and B
   ------------------------------------------------------------ */
void populateMatrices(double *A, double *B)
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i * N + j] = (double)(i + j + 1);
            B[i * N + j] = (double)(i - j + 1);
        }
    }
}


/* ------------------------------------------------------------
   Compute C = A * B
   ------------------------------------------------------------ */
void matrixMultiply(double *A, double *B, double *C)
{
    int i, j, k;

    #pragma omp parallel for private(j, k)
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            C[i * N + j] = 0.0;

            for (k = 0; k < N; k++)
            {
                C[i * N + j] +=
                    A[i * N + k] * B[k * N + j];
            }
        }
    }
}


/* ------------------------------------------------------------
   Compute Frobenius norm of matrix C
   ------------------------------------------------------------ */
double matrixNorm(double *C)
{
    double sum = 0.0;
    int i, j;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            double value = C[i * N + j];
            sum += value * value;
        }
    }

    return sqrt(sum);
}
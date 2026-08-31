#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#ifdef PROFILING_ON
#include <TAU.h>
#else

#define TAU_PROFILE_TIMER(timer, name, type, group)
#define TAU_PROFILE_START(timer)
#define TAU_PROFILE_STOP(timer)

#define TAU_USER 0

#endif

#include <pthread.h>
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

    #pragma omp parallel private(i, j)
    {
        TAU_PROFILE_TIMER(timer_populate,
                          "populateMatrices_thread_work",
                          "",
                          TAU_USER);

        TAU_PROFILE_START(timer_populate);

        #pragma omp for
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                A[i * N + j] = (double)(i + j + 1);
                B[i * N + j] = (double)(i - j + 1);
            }
        }

        TAU_PROFILE_STOP(timer_populate);
    }
}


/* ------------------------------------------------------------
   Compute C = A * B
   ------------------------------------------------------------ */
void matrixMultiplyOMP(double *A, double *B, double *C)
{
    int i, j, k;

    #pragma omp parallel private(i, j, k)
    {
        TAU_PROFILE_TIMER(timer_multiply,
                          "matrixMultiply_thread_work",
                          "",
                          TAU_USER);

        TAU_PROFILE_START(timer_multiply);

        #pragma omp for
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

        TAU_PROFILE_STOP(timer_multiply);
    }
}


typedef struct
{
    double *A;
    double *B;
    double *C;

    int start_row;
    int end_row;

} MatrixMultiplyThreadData;


/* ------------------------------------------------------------
   Pthread worker for matrix multiplication
   ------------------------------------------------------------ */
void *matrixMultiplyWorker(void *arg)
{
    MatrixMultiplyThreadData *data =
        (MatrixMultiplyThreadData *)arg;

    double *A = data->A;
    double *B = data->B;
    double *C = data->C;

    int i, j, k;

    TAU_PROFILE_TIMER(timer_multiply,
                      "matrixMultiply_thread_work",
                      "",
                      TAU_USER);

    TAU_PROFILE_START(timer_multiply);

    for (i = data->start_row; i < data->end_row; i++)
    {
        for (j = 0; j < N; j++)
        {
            C[i * N + j] = 0.0;

            for (k = 0; k < N; k++)
            {
                C[i * N + j] +=
                    A[i * N + k] *
                    B[k * N + j];
            }
        }
    }

    TAU_PROFILE_STOP(timer_multiply);

    return NULL;
}

/* ------------------------------------------------------------
   Compute C = A * B using Pthreads
   ------------------------------------------------------------ */
void matrixMultiply(double *A, double *B, double *C)
{
    char *env;
    int num_threads;
    int rows_per_thread;
    int remaining_rows;
    int start_row;
    int rows;
    int t;

    pthread_t *threads;
    MatrixMultiplyThreadData *threadData;


    /* Read number of Pthreads from environment variable */
    env = getenv("PTHREAD_NUM_THREADS");

    if (env != NULL)
    {
        num_threads = atoi(env);
    }
    else
    {
        num_threads = 1;
    }

    if (num_threads < 1)
    {
        num_threads = 1;
    }


    /* Allocate thread information */
    threads = (pthread_t *)
        malloc(num_threads * sizeof(pthread_t));

    threadData = (MatrixMultiplyThreadData *)
        malloc(num_threads *
               sizeof(MatrixMultiplyThreadData));

    if (threads == NULL || threadData == NULL)
    {
        printf("Pthread memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }


    /* Determine row distribution */
    rows_per_thread = N / num_threads;
    remaining_rows  = N % num_threads;

    start_row = 0;


    /* Create threads */
    for (t = 0; t < num_threads; t++)
    {
        rows = rows_per_thread;

        if (t < remaining_rows)
        {
            rows++;
        }

        threadData[t].A = A;
        threadData[t].B = B;
        threadData[t].C = C;

        threadData[t].start_row = start_row;
        threadData[t].end_row   = start_row + rows;

        pthread_create(&threads[t],
                       NULL,
                       matrixMultiplyWorker,
                       &threadData[t]);

        start_row += rows;
    }


    /* Wait for all threads to finish */
    for (t = 0; t < num_threads; t++)
    {
        pthread_join(threads[t], NULL);
    }


    /* Free temporary thread data */
    free(threads);
    free(threadData);
}



/* ------------------------------------------------------------
   Compute Frobenius norm of matrix C
   ------------------------------------------------------------ */
double matrixNorm(double *C)
{
    double sum = 0.0;
    int i, j;

    #pragma omp parallel private(i, j)
    {
        TAU_PROFILE_TIMER(timer_norm,
                          "matrixNorm_thread_work",
                          "",
                          TAU_USER);

        TAU_PROFILE_START(timer_norm);

        #pragma omp for reduction(+:sum)
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                double value = C[i * N + j];
                sum += value * value;
            }
        }

        TAU_PROFILE_STOP(timer_norm);
    }

    return sqrt(sum);
}

void matrixAdd(double *A, double *B, double *C)
{
    int i, j;

    #pragma omp parallel private(i, j)
    {
        TAU_PROFILE_TIMER(timer_add,
                          "matrixAdd_thread_work",
                          "",
                          TAU_USER);

        TAU_PROFILE_START(timer_add);

        #pragma omp for
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                C[i * N + j] =
                    A[i * N + j] + B[i * N + j];
            }
        }

        TAU_PROFILE_STOP(timer_add);
    }
}

void matrixScale(double *A, double *B, double alpha)
{
    int i, j;

    #pragma omp parallel for private(j)
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            B[i * N + j] =
                alpha * A[i * N + j];
        }
    }
}

void matrixTranspose(double *A, double *AT)
{
    int i, j;

    #pragma omp parallel for private(j)
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            AT[j * N + i] =
                A[i * N + j];
        }
    }
}

void matrixVectorMultiply(double *A, double *x, double *y)
{
    int i, j;

    #pragma omp parallel private(i, j)
    {
        TAU_PROFILE_TIMER(timer_matvec,
                          "matrixVectorMultiply_thread_work",
                          "",
                          TAU_USER);

        TAU_PROFILE_START(timer_matvec);

        #pragma omp for
        for (i = 0; i < N; i++)
        {
            double sum = 0.0;

            for (j = 0; j < N; j++)
            {
                sum += A[i * N + j] * x[j];
            }

            y[i] = sum;
        }

        TAU_PROFILE_STOP(timer_matvec);
    }
}


double matrixTrace(double *A)
{
    double trace = 0.0;
    int i;

    #pragma omp parallel private(i)
    {
        TAU_PROFILE_TIMER(timer_trace,
                          "matrixTrace_thread_work",
                          "",
                          TAU_USER);

        TAU_PROFILE_START(timer_trace);

        #pragma omp for reduction(+:trace)
        for (i = 0; i < N; i++)
        {
            trace += A[i * N + i];
        }

        TAU_PROFILE_STOP(timer_trace);
    }

    return trace;
}

double matrixFrobeniusProduct(double *A, double *B)
{
    double sum = 0.0;
    int i, j;

    #pragma omp parallel for private(j) reduction(+:sum)
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            sum += A[i * N + j] * B[i * N + j];
        }
    }

    return sum;
}

void matrixMultiplyTransposed(double *A, double *BT, double *C)
{
    int i, j, k;

    #pragma omp parallel for private(j, k)
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            double sum = 0.0;

            for (k = 0; k < N; k++)
            {
                sum +=
                    A[i * N + k] *
                    BT[j * N + k];
            }

            C[i * N + j] = sum;
        }
    }
}
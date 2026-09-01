# OpenMP Matrix Multiplication Profiling with TAU

This example demonstrates how to profile an OpenMP-based matrix multiplication program using the TAU Performance System.

The profiling workflow uses:

- OpenMP for shared-memory parallelism
- TAU's OpenMP Tools Interface (OMPT) support
- Event-Based Sampling (EBS)
- `pprof` for command-line inspection
- ParaProf for graphical visualization

---

## 1. OpenMP Matrix Multiplication

The program computes the product of two dense matrices,

\[
\mathbf{C} = \mathbf{A}\mathbf{B},
\]

where

\[
C_{ij}=\sum_{k=0}^{N-1} A_{ik} B_{kj}.\]

The matrices are dynamically allocated as one-dimensional arrays and accessed using row-major indexing. For example,

```c
A[i * N + j]
```

corresponds to the matrix entry \(A_{ij}\).

The matrix multiplication is parallelized using OpenMP by distributing iterations of the outer row loop among the available threads:

```c
#include <omp.h>

void matrixMultiply(double *A, double *B, double *C)
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
                sum += A[i * N + k] * B[k * N + j];
            }

            C[i * N + j] = sum;
        }
    }
}
```

Each OpenMP thread operates on a different subset of rows of the output matrix `C`.

Because different threads write to different rows, no synchronization is required during the matrix multiplication. The matrices `A` and `B` are shared among all threads and are only read during the computation.

The computational complexity remains

\[
\mathcal{O}(N^3),
\]

but the outer loop can now execute concurrently across multiple CPU threads.

A relatively large matrix, such as

```c
#define N 2000
```

is useful for profiling because the application must execute long enough for TAU to collect a representative number of samples.

---

## 2. Compiling the OpenMP Application

The application must be compiled with OpenMP support.

For GCC, the required compiler option is

```text
-fopenmp
```

Debug information should also be retained so that TAU can associate sampled instruction addresses with source-code locations.

Recommended compilation flags are

```makefile
CFLAGS = -O2 -g -Wall -Wextra -fopenmp
```

and the link flags should include

```makefile
LDFLAGS = -fopenmp -lm
```

For example:

```makefile
CC = g++-12

CFLAGS  = -O2 -g -Wall -Wextra -fopenmp
LDFLAGS = -fopenmp -lm
```

Rebuild the program using

```bash
make clean
make
```

---

## 3. Selecting the Number of OpenMP Threads

The number of OpenMP threads is controlled using the `OMP_NUM_THREADS` environment variable.

For example, to use four threads:

```bash
export OMP_NUM_THREADS=4
```

The OpenMP runtime configuration can be displayed using

```bash
export OMP_DISPLAY_ENV=TRUE
```

followed by

```bash
./bin/matrix_multiply.exe
```

The output should indicate that four OpenMP threads have been requested.

For profiling runs, the display option can be disabled afterward:

```bash
unset OMP_DISPLAY_ENV
```

---

## 4. TAU OpenMP Profiling with OMPT

The TAU installation used for this example contains OpenMP and OMPT support.

The relevant TAU runtime binding is

```text
shared-caltop-openmp-ompt-pdt-openmp
```

OMPT allows TAU to interact with the OpenMP runtime and identify OpenMP execution regions and worker threads.

Before profiling, enable full OMPT support:

```bash
export TAU_OMPT_SUPPORT_LEVEL=full
```

For example, for a four-thread run:

```bash
export OMP_NUM_THREADS=4
export TAU_OMPT_SUPPORT_LEVEL=full
```

---

## 5. Event-Based Sampling

TAU Event-Based Sampling is enabled using the `-ebs` option.

For the OpenMP application, EBS is used together with OMPT:

```bash
tau_exec -T ompt -ompt -ebs ./bin/matrix_multiply.exe
```

The options have the following roles:

### `-T ompt`

Selects a TAU configuration containing OMPT support.

### `-ompt`

Activates OpenMP runtime profiling through OMPT.

### `-ebs`

Enables Event-Based Sampling.

Together, these options allow TAU to sample the application while also identifying the OpenMP execution context.

---

## 6. Running a Profiling Experiment

Before each profiling run, remove profile files from previous experiments:

```bash
rm -f profile.*
```

Select the OpenMP thread count:

```bash
export OMP_NUM_THREADS=4
```

Enable full OMPT support:

```bash
export TAU_OMPT_SUPPORT_LEVEL=full
```

Then run:

```bash
tau_exec -T ompt -ompt -ebs ./bin/matrix_multiply.exe
```

The complete profiling sequence is therefore:

```bash
export OMP_NUM_THREADS=4
export TAU_OMPT_SUPPORT_LEVEL=full

rm -f profile.*

tau_exec -T ompt -ompt -ebs ./bin/matrix_multiply.exe
```

---

## 7. Verbose TAU Output

During setup or debugging, TAU can be run in verbose mode:

```bash
tau_exec -T ompt -ompt -ebs -v ./bin/matrix_multiply.exe
```

The verbose output shows which TAU libraries and bindings are loaded.

For this installation, the important line is

```text
Using shared-caltop-openmp-ompt-pdt-openmp
```

which confirms that TAU selected the OpenMP/OMPT-enabled configuration.

TAU may also report that no MPI binding is available and then fall back to the non-MPI configuration. This is expected because this matrix multiplication example does not use MPI.

---

## 8. TAU Profile Files

After execution, TAU generates profile files in the current working directory.

For an OpenMP run, multiple thread profiles may be produced, for example:

```text
profile.0.0.0
profile.0.0.1
profile.0.0.2
profile.0.0.3
```

TAU profile filenames generally follow the convention

```text
profile.<node>.<context>.<thread>
```

These files contain measurements for the execution contexts registered by TAU.

---

## 9. Inspecting the Profile with `pprof`

TAU provides the command-line utility `pprof`.

From the directory containing the generated `profile.*` files, run

```bash
pprof
```

This displays a text summary of the measured execution.

For this example, the matrix multiplication and associated OpenMP parallel region should dominate the application runtime.

`pprof` is useful for quickly checking that TAU generated sensible profiling data before opening the graphical interface.

---

## 10. Visualizing the Results with ParaProf

The profile can be visualized using TAU's graphical profiler:

```bash
paraprof
```

ParaProf reads the `profile.*` files in the current directory.

For the OpenMP application, ParaProf can be used to examine:

- individual OpenMP threads
- execution time per thread
- OpenMP parallel regions
- sampled code locations
- inclusive and exclusive time
- thread load balance
- runtime overhead

The matrix multiplication should appear as the principal computational hotspot.

Conceptually, the execution structure is

```text
main
 |
 +-- matrixMultiply()
       |
       +-- OpenMP parallel region
             |
             +-- Thread 0
             +-- Thread 1
             +-- Thread 2
             +-- Thread 3
```

---

## 11. Optional Call-Stack Sampling

If TAU was built with `libunwind` support, EBS can also collect call-stack information.

Enable stack unwinding with

```bash
export TAU_EBS_UNWIND=1
```

and optionally set the maximum unwind depth:

```bash
export TAU_EBS_UNWIND_DEPTH=10
```

Then run:

```bash
tau_exec -T ompt -ompt -ebs ./bin/matrix_multiply.exe
```

This can provide additional information about the calling context associated with sampled instructions.

---

## 12. Recommended Profiling Procedure

A typical OpenMP profiling experiment can be performed using:

```bash
# Build the application
make clean
make

# Select the number of OpenMP threads
export OMP_NUM_THREADS=4

# Enable TAU OMPT support
export TAU_OMPT_SUPPORT_LEVEL=full

# Remove previous profiling results
rm -f profile.*

# Run TAU with OMPT and Event-Based Sampling
tau_exec -T ompt -ompt -ebs ./bin/matrix_multiply.exe

# Inspect the profile from the terminal
pprof

# Visualize the profile
paraprof
```

For configuration debugging, use:

```bash
tau_exec -T ompt -ompt -ebs -v ./bin/matrix_multiply.exe
```

---

## 13. Profiling Different Thread Counts

The same procedure can be repeated for different numbers of OpenMP threads:

```bash
export OMP_NUM_THREADS=1
```

```bash
export OMP_NUM_THREADS=2
```

```bash
export OMP_NUM_THREADS=4
```

```bash
export OMP_NUM_THREADS=8
```

Before each run, remove or archive the previous TAU profile files:

```bash
rm -f profile.*
```

The resulting profiles can be compared to study:

- parallel speedup
- thread load balance
- OpenMP runtime overhead
- changes in sampled hotspots
- scaling behavior

For \(p\) threads, the parallel speedup is

\[
S_p = \frac{T_1}{T_p},
\]

where \(T_1\) is the runtime using one OpenMP thread and \(T_p\) is the runtime using \(p\) threads.

The corresponding parallel efficiency is

\[
E_p = \frac{S_p}{p}.
\]

---

## 14. Summary

The OpenMP profiling workflow is:

```text
OpenMP matrix multiplication
          |
          v
Compile with
-g -fopenmp
          |
          v
Set OMP_NUM_THREADS
          |
          v
Enable TAU OMPT support
TAU_OMPT_SUPPORT_LEVEL=full
          |
          v
tau_exec -T ompt -ompt -ebs
          |
          v
       profile.*
          |
     +----+----+
     |         |
     v         v
   pprof    ParaProf
```

This example provides a simple test case for understanding TAU profiling of OpenMP applications using OMPT and Event-Based Sampling.
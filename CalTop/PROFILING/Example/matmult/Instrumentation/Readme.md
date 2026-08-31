# TAU Profiling for OpenMP and Pthreads

This repository provides a simple example of profiling shared-memory parallel
applications using the **TAU Performance System**. The example demonstrates
manual TAU instrumentation for code containing both **OpenMP** and **POSIX
Threads (Pthreads)**.

The numerical operations used in the example are not important. The primary
purpose of the repository is to demonstrate how TAU can be incorporated into
a C application and used to measure parallel regions executed by individual
threads.

## Requirements

The following software is required:

- GCC
- OpenMP
- Pthreads
- TAU Performance System
- ParaProf (optional, for graphical profile visualization)

The TAU installation used by the Makefile is specified through:

```makefile
TAU_ROOT = /path/to/TAU_install
```

The appropriate TAU makefile must also be available in the environment. For
example:

```bash
export TAU_MAKEFILE=/path/to/TAU_install/x86_64/lib/Makefile.tau-...
```

## Enabling and Disabling Profiling

Profiling is controlled using the `PROFILING` Makefile variable.

A normal build is performed using:

```bash
make clean
make
```

In this case, TAU manual instrumentation is disabled.

To enable TAU profiling:

```bash
make clean
make PROFILING=1
```

When profiling is enabled, the Makefile defines:

```text
PROFILING_ON
```

which allows the source code to conditionally include TAU:

```c
#ifdef PROFILING_ON
#include <TAU.h>
#endif
```

TAU profiling calls can similarly be disabled when `PROFILING_ON` is not
defined. This allows the same source code to be compiled either with or
without TAU instrumentation.

## Manual TAU Instrumentation

A manually instrumented region is defined using:

```c
TAU_PROFILE_TIMER(timer,
                  "region_name",
                  "",
                  TAU_USER);
```

The region to be measured is enclosed by:

```c
TAU_PROFILE_START(timer);

/* Code being profiled */

TAU_PROFILE_STOP(timer);
```

This allows TAU to measure selected computational regions rather than relying
only on automatic function-level instrumentation.

## Profiling OpenMP Regions

For OpenMP code, a TAU timer can be placed inside the OpenMP parallel region:

```c
#pragma omp parallel
{
    TAU_PROFILE_TIMER(timer,
                      "parallel_thread_work",
                      "",
                      TAU_USER);

    TAU_PROFILE_START(timer);

    #pragma omp for
    for (int i = 0; i < N; i++)
    {
        /* Parallel work */
    }

    TAU_PROFILE_STOP(timer);
}
```

Because every OpenMP thread enters the parallel region, the timer measures
the work performed by each participating thread.

The number of OpenMP threads is controlled using:

```bash
export OMP_NUM_THREADS=6
```

## Profiling Pthreads

Pthreads require explicit thread creation and work assignment. A TAU timer
can therefore be placed inside the Pthread worker function:

```c
void *worker(void *arg)
{
    TAU_PROFILE_TIMER(timer,
                      "pthread_thread_work",
                      "",
                      TAU_USER);

    TAU_PROFILE_START(timer);

    /* Work assigned to this Pthread */

    TAU_PROFILE_STOP(timer);

    return NULL;
}
```

Each Pthread executes the worker function independently, allowing TAU to
record the execution time associated with each thread.

Unlike OpenMP, Pthreads does not define a standard environment variable for
the number of threads. This example uses the user-defined variable:

```bash
export PTHREAD_NUM_THREADS=6
```

The application reads this variable and explicitly creates the requested
number of Pthreads.

## Using OpenMP and Pthreads Together

The application may contain both OpenMP and Pthread parallel regions. Their
thread counts can be controlled independently:

```bash
export OMP_NUM_THREADS=6
export PTHREAD_NUM_THREADS=6
```

`OMP_NUM_THREADS` is interpreted by the OpenMP runtime, whereas
`PTHREAD_NUM_THREADS` is read explicitly by the application.

The executable must be linked with support for both threading models:

```text
-fopenmp -pthread
```

## Running a Profile

After building with profiling enabled:

```bash
make clean
make PROFILING=1
```

set the desired thread counts:

```bash
export OMP_NUM_THREADS=6
export PTHREAD_NUM_THREADS=6
```

and run the application:

```bash
./bin/matrix_multiply.exe
```

TAU will generate profile files of the form:

```text
profile.*
```

## Viewing Results

A text summary of the TAU profile can be viewed using:

```bash
pprof
```

For graphical analysis, use:

```bash
paraprof
```

ParaProf can be used to examine the execution time associated with individual
threads and manually instrumented computational regions.

## Typical Workflow

A complete profiling run therefore consists of:

```bash
export TAU_MAKEFILE=/path/to/TAU/Makefile.tau-...

make clean
make PROFILING=1

export OMP_NUM_THREADS=6
export PTHREAD_NUM_THREADS=6

./bin/matrix_multiply.exe

pprof
```

or:

```bash
paraprof
```

To return to a non-profiled build:

```bash
make clean
make
```

This provides a simple mechanism for maintaining a single code base that can
be compiled either as a normal application or as a manually instrumented TAU
profiling build.
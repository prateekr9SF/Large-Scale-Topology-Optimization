# TAU Profiling Setup for CalTop

This directory contains installation and configuration notes for using the **TAU Performance System** with **CalTop**.

CalTop is a mixed-language application containing **C and Fortran** source code. Its parallel execution model uses a combination of **OpenMP** and **POSIX threads (pthreads)**. TAU is configured here with OpenMP support, OMPT runtime instrumentation, and source-level instrumentation through the **Program Database Toolkit (PDT)**.

The configuration documented below was tested on an Ubuntu workstation using the GCC 12 toolchain.

---

## Software Configuration

| Package     | Version / Configuration |
| ----------- | ----------------------- |
| TAU         | 2.35.2                  |
| PDT         | 3.25.2                  |
| GCC         | 12.5.0                  |
| G++         | 12.5.0                  |
| GNU Fortran | GCC 12 toolchain        |
| OpenMP      | Enabled                 |
| pthreads    | Used by CalTop          |
| OMPT        | Enabled                 |
| BFD         | Enabled                 |
| libunwind   | Enabled                 |

The installation uses the following directory structure:

```text
~/Documents/Software/
├── PDT/
│   ├── pdtoolkit-3.25.2/
│   └── PDT_install/
│
└── TAU/
    ├── tau-2.35.2/
    └── TAU_install/
```

The exact installation paths can be changed as needed.

---

# 1. Installing PDT

PDT is used by TAU for automatic source-code instrumentation.

Download and extract PDT 3.25.2 and enter the source directory:

```bash
cd ~/Documents/Software/PDT/pdtoolkit-3.25.2
```

The installation prefix used here is:

```text
/home/prateek/Documents/Software/PDT/PDT_install
```

Configure PDT with:

```bash
./configure \
    -prefix=/home/prateek/Documents/Software/PDT/PDT_install/
```

## 1.1 Building PDT with GCC 12

On systems containing multiple GCC versions, the PDT-generated Makefiles may use `/usr/bin/g++` directly rather than respecting the shell `CXX` variable.

For example, the generated `ductape/Makefile` may contain:

```make
PDT_GXX=/usr/bin/g++
```

To guarantee that PDT is built with GCC 12, explicitly override `PDT_GXX`:

```bash
make PDT_GXX=/usr/bin/g++-12 install
```

This ensures that PDT components, including `tau_instrumentor`, are compiled using G++ 12.

The resulting installation is located under:

```text
/home/prateek/Documents/Software/PDT/PDT_install/x86_64/
```

Set the PDT environment:

```bash
export PDT_HOME=/home/prateek/Documents/Software/PDT/PDT_install
export PATH=$PDT_HOME/x86_64/bin:$PATH
```

Verify the installation:

```bash
which tau_instrumentor
which pdbconv
which pdbtree
```

---

# 2. Installing TAU

Download and extract TAU 2.35.2 and enter the source directory:

```bash
cd ~/Documents/Software/TAU/tau-2.35.2
```

The installation prefix used here is:

```text
/home/prateek/Documents/Software/TAU/TAU_install
```

## 2.1 TAU Configuration for CalTop

CalTop contains C and Fortran source code and uses both OpenMP and pthreads.

For the present installation, TAU uses **OpenMP as its thread package**. OMPT is enabled to provide OpenMP runtime-level instrumentation. PDT is enabled for source-level instrumentation of CalTop routines.

Configure TAU using:

```bash
./configure \
    -prefix=/home/prateek/Documents/Software/TAU/TAU_install \
    -cc=gcc-12 \
    -c++=g++-12 \
    -fortran=gnu \
    -pdt=/home/prateek/Documents/Software/PDT/PDT_install \
    -pdt_c++=g++-12 \
    -openmp \
    -ompt=download \
    -bfd=download \
    -unwind=download \
    -tag=caltop-openmp
```

The resulting configuration provides:

```text
CalTop
│
├── C
│   └── GCC 12
│
├── Fortran
│   └── GNU Fortran
│
├── Parallelism
│   ├── OpenMP
│   │   └── OMPT instrumentation
│   │
│   └── POSIX threads (pthreads)
│       └── Application-level threading
│
└── TAU
    ├── PDT
    │   └── Source instrumentation
    ├── BFD
    │   └── Symbol/address resolution
    └── libunwind
        └── Stack unwinding
```

### OpenMP and pthreads

The `-openmp` option selects OpenMP as the threading package used by this TAU configuration.

CalTop may still contain and execute pthread-based code. The presence of pthreads in the application does not, by itself, require replacing the TAU OpenMP configuration with a `-pthread` configuration.

If detailed TAU instrumentation of pthread creation and synchronization is required in addition to OpenMP/OMPT profiling, a separate TAU pthread configuration or additional TAU runtime instrumentation can be investigated.

For the primary CalTop profiling configuration, OpenMP + OMPT is used because OpenMP runtime behavior is an important component of CalTop's parallel execution.

---

# 3. Third-Party Packages

Several optional packages are downloaded automatically during TAU configuration.

## libunwind

```bash
-unwind=download
```

enables stack unwinding.

TAU should report:

```text
NOTE: Enabling Stack Unwinding support
NOTE: Using libunwind Unwinder.
```

When built with the GCC 12 configuration, TAU should also report:

```text
Using unwind_c_compiler=gcc-12
Using unwind_cxx_compiler=g++-12
```

## BFD

```bash
-bfd=download
```

downloads/builds GNU BFD support.

TAU should report:

```text
NOTE: Using BFD support
NOTE: Using ELF support in BFD
NOTE: Using DEMANGLE support
```

BFD provides binary address and symbol resolution.

## OMPT

```bash
-ompt=download
```

enables the OpenMP Tools Interface.

If the system OpenMP runtime does not provide the required OMPT functionality, TAU may download and build an LLVM OpenMP runtime.

A successful configuration should report:

```text
Using OMPT 5
```

and:

```text
NOTE: Using OMPT OpenMP options for OMPT 5.0 specification
```

OMPT and PDT serve complementary purposes:

```text
                 CalTop
                   │
        ┌──────────┴──────────┐
        │                     │
        ▼                     ▼
       PDT                   OMPT
        │                     │
        ▼                     ▼
Source instrumentation    OpenMP runtime
        │                 instrumentation
        │                     │
        ▼                     ▼
CalTop functions         Parallel regions
Solver routines          Threads
Sensitivities            Barriers
FEA routines             Synchronization
Optimization routines    Tasks
```

---

# 4. Verify the TAU Configuration

Before compiling TAU, inspect the output from `configure`.

The C/C++ compiler detection should report GCC 12:

```text
Default C++ compiler will be g++ version 12.5.0
Default C compiler will be gcc version 12.5.0

Using GNU lib dir as /usr/lib/gcc/x86_64-linux-gnu/12/
Using GNU stdc++ lib dir as /usr/lib/gcc/x86_64-linux-gnu/12/
```

Fortran support should also be enabled:

```text
Setting F90 compiler based on requested: gfortran
Default Fortran compiler will be GNU gfortran
```

The final configuration summary should contain entries similar to:

```text
NOTE: GNU gfortran compiler specific options used
NOTE: Using the specified Fortran compiler instead of the default
NOTE: Using PDT for TAU Source Code Instrumentation ***
NOTE: Using OpenMP Threads ***
NOTE: Using GNU's OpenMP options ***
NOTE: Using OMPT OpenMP options for OMPT 5.0 specification
NOTE: Using BFD support
NOTE: Enabling Stack Unwinding support
NOTE: Using libunwind Unwinder.
```

PDT should additionally report:

```text
PDT supports Fortran Loop Level information
```

Do not proceed if TAU unexpectedly detects a different GCC version.

---

# 5. Build and Install TAU

Once configuration completes successfully:

```bash
make install
```

The TAU installation is located at:

```text
/home/prateek/Documents/Software/TAU/TAU_install
```

with architecture-specific executables under:

```text
/home/prateek/Documents/Software/TAU/TAU_install/x86_64/bin
```

---

# 6. Environment Setup

Add PDT and TAU to the shell environment:

```bash
export PDT_HOME=/home/prateek/Documents/Software/PDT/PDT_install
export TAU_ROOT=/home/prateek/Documents/Software/TAU/TAU_install

export PATH=$PDT_HOME/x86_64/bin:$PATH
export PATH=$TAU_ROOT/x86_64/bin:$PATH
```

These definitions may be added to `~/.bashrc` or `~/.zshrc`.

For example, for Zsh:

```bash
source ~/.zshrc
```

---

# 7. Verify the Installation

Verify the TAU compiler wrappers and runtime tools:

```bash
which tau_cc.sh
which tau_cxx.sh
which tau_f90.sh
which tau_exec
```

The returned paths should point to:

```text
/home/prateek/Documents/Software/TAU/TAU_install/x86_64/bin/
```

Inspect the installed TAU configurations with:

```bash
ls $TAU_ROOT/x86_64/lib/Makefile.tau*
```

The CalTop configuration should contain the `caltop-openmp` tag together with the enabled PDT/OpenMP/OMPT configuration.

---

# 8. Compiler Consistency

The workstation may contain multiple GCC versions. Verify the compiler toolchain before rebuilding PDT or TAU:

```bash
gcc-12 --version
g++-12 --version
gfortran-12 --version
```

The TAU configuration should select GCC 12 libraries:

```text
/usr/lib/gcc/x86_64-linux-gnu/12/
```

This is particularly important on systems where `/usr/bin/gcc`, `/usr/bin/g++`, or `/usr/bin/gfortran` refer to a newer system compiler.

---

# 9. CalTop Profiling Workflow

The intended profiling workflow is:

```text
CalTop C / Fortran source
          │
          ▼
     PDT instrumentation
          │
          ▼
   TAU compiler wrappers
          │
          ├── C       → GCC 12
          └── Fortran → GNU Fortran
          │
          ▼
  Instrumented CalTop executable
          │
          ├── OpenMP execution
          ├── OMPT runtime events
          └── pthread-based execution
          │
          ▼
       TAU profiles
```

The CalTop build system can subsequently be configured to switch between the normal GCC 12 build and a TAU-instrumented build without modifying the underlying source files.

---

# References

* TAU Performance System: https://www.cs.uoregon.edu/research/tau/
* TAU Installation Guide: https://www.cs.uoregon.edu/research/tau/docs/html-docs/latest/installguide/installguide.html
* Program Database Toolkit (PDT): https://www.cs.uoregon.edu/research/pdt/

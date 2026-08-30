# TAU Profiling Setup for CalTop

This directory contains notes and configuration information for building and using the **TAU Performance System** with **CalTop**.

CalTop contains both **C and Fortran** source code and uses **OpenMP** for shared-memory parallelism. The TAU installation described here enables source-level instrumentation through the **Program Database Toolkit (PDT)** and OpenMP runtime instrumentation through **OMPT**.

The configuration documented below was tested on an Ubuntu workstation using the GCC 12 toolchain.

---

## 1. Software Configuration

The following software versions were used:

| Package     | Version / Configuration |
| ----------- | ----------------------- |
| TAU         | 2.35.2                  |
| PDT         | 3.25.2                  |
| GCC         | 12.5.0                  |
| G++         | 12.5.0                  |
| GNU Fortran | GCC 12 toolchain        |
| OpenMP      | Enabled                 |
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

The exact paths can be changed as needed.

---

# 2. Installing PDT

PDT is used by TAU for automatic source-code instrumentation.

Download and extract PDT 3.25.2, then enter the source directory:

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

## 2.1 Building PDT with GCC 12

On systems where the default `/usr/bin/g++` points to a newer compiler, the PDT-generated Makefiles may ignore the shell `CXX` environment variable and use `/usr/bin/g++` directly.

For example, the generated `ductape/Makefile` may contain:

```make
PDT_GXX=/usr/bin/g++
```

Therefore, explicitly override `PDT_GXX` when building:

```bash
make PDT_GXX=/usr/bin/g++-12 install
```

This ensures that PDT and `tau_instrumentor` are compiled using G++ 12.

The resulting installation is located under:

```text
/home/prateek/Documents/Software/PDT/PDT_install/x86_64/
```

with executables under:

```text
/home/prateek/Documents/Software/PDT/PDT_install/x86_64/bin
```

Verify the installation:

```bash
export PDT_HOME=/home/prateek/Documents/Software/PDT/PDT_install
export PATH=$PDT_HOME/x86_64/bin:$PATH

which tau_instrumentor
which pdbconv
which pdbtree
```

---

# 3. Installing TAU

Download and extract TAU 2.35.2, then enter the TAU source directory:

```bash
cd ~/Documents/Software/TAU/tau-2.35.2
```

The installation prefix used here is:

```text
/home/prateek/Documents/Software/TAU/TAU_install
```

## 3.1 TAU Configuration for CalTop

CalTop requires profiling support for:

* C
* Fortran
* OpenMP
* PDT source instrumentation
* OMPT OpenMP instrumentation

Additional TAU functionality is enabled through BFD and libunwind.

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

### Important

TAU's `-cc` and `-c++` options expect compiler names such as:

```text
gcc-12
g++-12
```

rather than absolute compiler paths such as:

```text
/usr/bin/gcc-12
/usr/bin/g++-12
```

Using the latter may cause TAU to reject the compiler selection and silently fall back to the system-default compiler.

---

# 4. Verify the TAU Configuration

Before compiling TAU, inspect the output from `configure`.

For this installation, TAU correctly detected:

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

PDT should also report Fortran support:

```text
PDT supports Fortran Loop Level information
```

Do not proceed with the TAU build if the configuration unexpectedly reports GCC/G++ from another compiler version.

---

# 5. Third-Party Packages

Several optional packages are downloaded automatically during the TAU configuration.

## libunwind

The option

```bash
-unwind=download
```

downloads and builds `libunwind`.

For the GCC 12 configuration, TAU should report something similar to:

```text
Using unwind_c_compiler=gcc-12
Using unwind_cxx_compiler=g++-12
```

libunwind provides stack-unwinding functionality useful for call-stack and sampling-based performance analysis.

## BFD

The option

```bash
-bfd=download
```

downloads GNU binutils and builds BFD support.

TAU should report:

```text
NOTE: Using BFD support
NOTE: Using ELF support in BFD
NOTE: Using DEMANGLE support
```

BFD allows TAU to resolve binary addresses and symbols, which is useful for obtaining meaningful function-level profiling information.

## OMPT

The option

```bash
-ompt=download
```

enables the OpenMP Tools Interface.

If the system OpenMP runtime does not provide the required OMPT interface, TAU automatically downloads and builds an LLVM OpenMP runtime.

A successful configuration should report:

```text
Using OMPT 5
```

and later:

```text
NOTE: Using OMPT OpenMP options for OMPT 5.0 specification
```

OMPT complements PDT instrumentation. PDT instruments CalTop source routines, while OMPT exposes events from the OpenMP runtime such as parallel regions, threads, synchronization, barriers, and tasks.

---

# 6. Build and Install TAU

Once configuration has completed successfully:

```bash
make install
```

The resulting TAU installation is located under:

```text
/home/prateek/Documents/Software/TAU/TAU_install
```

For the `x86_64` architecture, TAU executables are installed under:

```text
/home/prateek/Documents/Software/TAU/TAU_install/x86_64/bin
```

---

# 7. Environment Setup

Add PDT and TAU to the shell environment:

```bash
export PDT_HOME=/home/prateek/Documents/Software/PDT/PDT_install
export TAU_ROOT=/home/prateek/Documents/Software/TAU/TAU_install

export PATH=$PDT_HOME/x86_64/bin:$PATH
export PATH=$TAU_ROOT/x86_64/bin:$PATH
```

These commands can optionally be added to `~/.bashrc` or `~/.zshrc`.

Reload the shell configuration if necessary:

```bash
source ~/.bashrc
```

or:

```bash
source ~/.zshrc
```

---

# 8. Verify the TAU Installation

Verify that the TAU compiler wrappers and utilities are available:

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

The available TAU configurations can be inspected using:

```bash
ls $TAU_ROOT/x86_64/lib/Makefile.tau*
```

The CalTop-specific configuration should contain the `caltop-openmp` tag and the enabled OpenMP/OMPT/PDT options.

---

# 9. TAU Configuration Used for CalTop

The complete configuration is summarized below:

```text
CalTop
│
├── C
│   └── GCC 12
│
├── Fortran
│   └── GNU Fortran / GCC 12
│
├── OpenMP
│   ├── TAU OpenMP thread support
│   └── OMPT 5 runtime instrumentation
│
└── TAU 2.35.2
    ├── PDT 3.25.2
    │   └── Source-code instrumentation
    ├── BFD
    │   └── Symbol/address resolution
    └── libunwind
        └── Stack unwinding
```

The corresponding TAU configuration command is:

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

---

# 10. Notes

### Compiler consistency

The workstation may contain multiple GCC versions. Always verify the compiler toolchain before rebuilding TAU or PDT:

```bash
gcc-12 --version
g++-12 --version
gfortran-12 --version
```

TAU configuration output should also be inspected to ensure that GCC 12 libraries are being selected:

```text
/usr/lib/gcc/x86_64-linux-gnu/12/
```

### Reconfiguring TAU

When changing important TAU configuration options, such as enabling Fortran support or changing compilers, perform a clean rebuild before reconfiguring.

### Multiple TAU configurations

TAU supports multiple configurations within the same installation. The

```bash
-tag=caltop-openmp
```

option is used to distinguish the CalTop OpenMP configuration from other TAU configurations that may be built later.

---

# 11. References

* TAU Performance System: https://www.cs.uoregon.edu/research/tau/
* TAU Installation Guide: https://www.cs.uoregon.edu/research/tau/docs/html-docs/latest/installguide/installguide.html
* PDT: https://www.cs.uoregon.edu/research/pdt/

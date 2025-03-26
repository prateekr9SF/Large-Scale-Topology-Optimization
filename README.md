# CalTop: A CalculiX-based Topology Optimization Framework

## Overview
CalTop is a high-performance topology optimization framework built upon **CalculiX 2.15**, featuring advanced solvers such as **SPOOLES, PARDISO, SUPERLU, and PASTIX**. It implements a **density-based finite element framework** interfaced with **IPOPT** for gradient-based optimization. Additionally, CalTop integrates with **SU2_CFD** to solve **static aero-elastic topology optimization problems**.

## Features
- **Finite Element Analysis (FEA):** Uses CalculiX for **linear** structural analysis.
- **Fast Solvers:** Supports SPOOLES, PARDISO, SUPERLU, and PASTIX for fast stiffness matrix factorizations.
- **Density-Based Topology Optimization:** Tetrahedral elements parameterized with densities.
- **IPOPT Integration:** Utilizes the Interior Point OPTimizer for optimization.
- **Sensitivities:** Computes analytical sensitivities for mass, compliance, center of gravity and material stress
- **SU2_CFD Integration:** Utilized preCICE coupling adapter for force-displacement and adjoint sensitivties  tranfer to/from SU2.


## File Structure
```
├── CalTop/                # Density-based Finite Element Analysis     
    ├── src/               # Source code density-based CalculiX 2.15
│       ├── ccx_2.15.c     # Main CalculiX driver
│       ├── add_file_1     # Optimization routine
│       ├── add_file_2     # Density filtering and sensitivity analysis
│       └── add_file_3     # Utility functions
    ├── include/
         ├── ccx_2.15.h     # Main CalculiX driver header
├── CalGeo/                 # Mesh passive element identifier
├── Deps/                   # Dependencies (ARPACK,SPOOLES,YAML)
├── TestCases/              # Test cases for topology optimization
    ├── RAE2822/            # RAE2822 3D wing section

└── README.md               # This file
```
## Build Instructions

1. Clone the repository:
   ```sh
   git clone https://github.com/prateekr9SF/Large-Scale-Topology-Optimization.git
   ```
   This is the `ROOT` directory

2. Install dependency ARPACK:
   ``` sh
   wget https://web.archive.org/web/20220526222500fw_/https://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz
   wget https://web.archive.org/web/20220526222500fw_/https://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz
   ```
   Open ARmake.inc and make the following changes:
   - **Line 28:** Change `home = $(HOME)/ARPACK` to the path where ARPACK is extracted.
   - **Line 115:** Change `MAKE=/bin/make` to `MAKE=make`
   - **Line 120:** Change `SHELL =/bin/sh` to `SHELL=sh`
   - **Lines 104-105:** Set Fortran compiler flags:
   ```sh
   FC = gfortran
   ```
   ***Line 35:*** Set platform to INTEL (if applicable)
   ```sh
   PLAT = INTEL
   ```
   Open `UTIL/second.f' and comment out line 24:
   ```sh
   * EXTERNAL  ETIME
   ```
   Now build ARPACK using:
   ```sh
   make lib
   ```
3. Install dependency yamp-cpp:
   Get the latest verion of yamp-cp and build as a shared library:
   ```sh
   wget https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.2.zip
   unzip yaml-cpp-0.6.2.zip
   cd yaml-cpp-yaml-cpp-0.6.2
   mkdir build
   cd build
   cmake -DBUILD_SHARED_LIBS=ON ..
   make 
   ```
   After building, make sure to set `LD_BIBRARY_PATH` to the installation directory

4. CalTop currently supports **SPOOLES** and **INTEL MKL PARDISO** for matrix factorization.
   #### Option A: SPOOLES (Single-thread build)
   
   ```sh
   wget http://www.netlib.org/linalg/spooles/spooles.2.2.tgz
   mkdir SPOOLES.2.2
   tar zxvf spooles.2.2.tgz -C SPOOLES.2.2
   cd SPOOLES.2.2
   ```
   Edit `Make.inc` to set compiler version:
   ```sh
   CC=gcc
   ```
   Build SPOOLES
   ```sh
   make lib
   ```

   Navigate to the `SPOOLES_MAKE` directory and move the Makefile to `ROOT`.
   Open the `Makefile` and edit paths:
   ```sh
   SPOOLES_PATH = <spooles_installation_dir/src>
   ARPACK_PATH = <ARPACK installation_dir>
   ```
   #### Option B: INTEL MKL PARDISO (Multi-thread build)

   Install the Intel oneAPI Base Toolkit
   https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html

   Install the Intel oneAPI HPC Toolkit
   https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html

   Navigate to the `PARDISO_MAKE` directory and move the Makefile to `ROOT`.
   Open the `Makefile` and edit paths:
   ```sh
   ARPACK_PATH = <ARPACK installation_dir>
   MKL_LIB = <oneAPI_installation_path/intel/oneapi/mkl/year/lib/intel64>
   MKL_INCLUDE = <oneAPI_installation_path/intel/oneapi/mkl/year/include>
   MKL_INCLUDE = <oneAPI_installation_path/intel/oneapi/compiler/year/bin>
   ```
   `NOTE: Defaulty instalation directory for Intel MKL library is opt/`

6. Build and install    `CalTop`

   If installing in current directory:
   ```sh
   make 
   make install
   ```

   If using custom installation directory:
   ```sh
   make install PREFIX=$HOME/<installation_dir>
   ```

   Note: Run `make` or `make -j N` to build with `N` CPUs.

6. Set CalTop path
    In your `.bashrc` or `.profile` set:
   ```sh
   CALTOP_PATH=<installation_dir>
   export PATH=$CALTOP_PATH:$PATH
   ```
   source `.bashrc`

7. Uninstalling CalTop:
   
   If installed in default (current) directory:
   ```sh 
   make uninstall
   ```

   If installed in custom directory:
   ```sh
   make install PREFIX=$HOME/<installation_dir>
   ```


## Usage
TODO

## License
This project is licensed under [MIT License](LICENSE).

## Contact
For inquiries, please reach out to **Prateek Ranjan** at `prateekr@mit.edu`.

## Acknowledgement
This work was supported by NASA under award number **80NSSC19M0125** as part of the **C**enter for **H**igh-**E**fficiency **E**lectrical **T**echnologies for **A**ircraft **(CHEETA)**.
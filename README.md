# CalTop: CalculiX-based Topology Optimization Framework

CalTop is a density-based topology optimization framework built on CalculiX 2.15. It combines finite element analysis, density filtering, and adjoint sensitivities with gradient-based optimization. The CalFSI and CalADJ components extend the framework to coupled aeroelastic analysis and derivatives with SU2.

The repository contains three main components:

- **CalGeo** processes SU2 solid meshes and generates mesh and CalculiX node-set files for loads, supports, and designated surface regions.
- **CalFilt** constructs the density filter used during optimization, with parallel assembly for large tetrahedral meshes.
- **CalTop** evaluates structural responses and design sensitivities using a multithreaded CalculiX-based solver. External optimization drivers can use these results to update the design.
- **CalFSI** uses [preCICE](https://precice.org/) to couple CalTop with [SU2](https://github.com/prateekr9SF/SU2/) for fluid–structure interaction and static aeroelastic analysis.
- **CalADJ** couples CalTop with SU2 to compute coupled aeroelastic derivatives for gradient-based optimization.

## Repository layout

- `CalGeo/`: mesh preprocessing utilities.
- `CalTop/`: CalculiX-based analysis and topology optimization source.
- `CalFSI/`: preCICE-based CalTop–SU2 coupling for aeroelastic analysis.
- `CalADJ/`: CalTop–SU2 coupling for aeroelastic derivatives.
- `Deps/`: bundled or referenced dependencies.
- `TestCases/`: example problems, including the RAE 2822 wing section.
- `SPOOLES_MAKE/`, `PARDISO_MAKE/`: build configurations for the corresponding solvers.

## Build Instructions

1. Clone the repository:

   ```sh
   git clone https://github.com/prateekr9SF/Large-Scale-Topology-Optimization.git
   ```

   This is the `ROOT` directory.

2. Install dependency ARPACK:

   ```sh
   wget https://web.archive.org/web/20220526222500fw_/https://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz
   wget https://web.archive.org/web/20220526222500fw_/https://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz
   ```

   Open `ARmake.inc` and make the following changes:

   - **Line 28:** Change `home = $(HOME)/ARPACK` to the path where ARPACK is extracted.
   - **Line 115:** Change `MAKE=/bin/make` to `MAKE=make`.
   - **Line 120:** Change `SHELL =/bin/sh` to `SHELL=sh`.
   - **Lines 104–105:** Set the Fortran compiler:

     ```make
     FC = gfortran
     ```

   - **Line 35:** Set the platform to INTEL, if applicable:

     ```make
     PLAT = INTEL
     ```

   Open `UTIL/second.f` and comment out line 24:

   ```fortran
   * EXTERNAL  ETIME
   ```

   Build ARPACK:

   ```sh
   make lib
   ```

3. Install dependency yaml-cpp:

   Get yaml-cpp and build it as a shared library:

   ```sh
   wget https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.2.zip
   unzip yaml-cpp-0.6.2.zip
   cd yaml-cpp-yaml-cpp-0.6.2
   mkdir build
   cd build
   cmake -DBUILD_SHARED_LIBS=ON ..
   make
   ```

   After building, set `LD_BIBRARY_PATH` to the installation directory.

4. Choose a matrix factorization solver. CalTop supports **SPOOLES** and **Intel MKL PARDISO**.

   **Option A: SPOOLES (single-thread build)**

   ```sh
   wget http://www.netlib.org/linalg/spooles/spooles.2.2.tgz
   mkdir SPOOLES.2.2
   tar zxvf spooles.2.2.tgz -C SPOOLES.2.2
   cd SPOOLES.2.2
   ```

   Edit `Make.inc` to set the compiler:

   ```make
   CC=gcc
   ```

   Build SPOOLES:

   ```sh
   make lib
   ```

   Navigate to the `SPOOLES_MAKE` directory and move the Makefile to `ROOT`. Edit the paths in the Makefile:

   ```make
   SPOOLES_PATH = <spooles_installation_dir/src>
   ARPACK_PATH = <ARPACK installation_dir>
   ```

   **Option B: Intel MKL PARDISO (multi-thread build)**

   Install the [Intel oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html) and the [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html).

   Navigate to the `PARDISO_MAKE` directory and move the Makefile to `ROOT`. Edit the paths in the Makefile:

   ```make
   ARPACK_PATH = <ARPACK installation_dir>
   MKL_LIB = <oneAPI_installation_path/intel/oneapi/mkl/year/lib/intel64>
   MKL_INCLUDE = <oneAPI_installation_path/intel/oneapi/mkl/year/include>
   MKL_INCLUDE = <oneAPI_installation_path/intel/oneapi/compiler/year/bin>
   ```

   Note: The default installation directory for Intel oneAPI is `opt/`.

5. Build and install `CalTop`.

   For installation in the current directory:

   ```sh
   make
   make install
   ```

   For a custom installation directory:

   ```sh
   make install PREFIX=$HOME/<installation_dir>
   ```

   Run `make` or `make -j N` to build with `N` CPUs.

6. Set the CalTop path. In your `.bashrc` or `.profile`, set:

   ```sh
   CALTOP_PATH=<installation_dir>
   export PATH=$CALTOP_PATH:$PATH
   ```

   Source `.bashrc`.

7. Uninstall CalTop.

   If installed in the default (current) directory:

   ```sh
   make uninstall
   ```

   If installed in a custom directory:

   ```sh
   make install PREFIX=$HOME/<installation_dir>
   ```

## Usage: calFilt

`calFilt` assembles a density filter for CalTop. Set the number of OpenMP threads, then run it with the mesh filename stem, filter radius, and filter-kernel storage parameter:

```sh
export OMP_NUM_THREADS=N
calFilt.exe -i <filename_without_extension> -r <filter_radius> -f <number_of_nonzeros_in_filter_kernel>
```

Choose `N` according to the CPUs allocated to the process.

## Usage: calGeo

`calGeo` reads an SU2 mesh and prepares files used by the structural analysis. You can call the Python script directly or define a shell alias:

```sh
alias calGeo='python3 /path/to/Large-Scale-Topology-Optimization/CalGeo/calGeo.py'
```

After reloading your shell configuration, run it with the mesh and marker arguments appropriate to the case. For example, to identify several skin markers as passive regions:

```sh
calGeo mesh_name.su2 SkinMarkerList skin1 skin3 skin10
```

The marker names must match those in the SU2 mesh. The preprocessing step generates the `.nam` and `.msh` files required by CalTop.

## Usage: calTop

CalTop reads a CalculiX input case and element densities from `density.dat` (or the executable's default densities). Set the OpenMP thread count before running a shared-memory build:

```sh
export OMP_NUM_THREADS=<number_of_threads>
```

**Structural analysis:**

```sh
calTop.exe <filename_without_extension>
```

This evaluates the linear elastic response and writes `elastic_Field.vtu` for visualization of densities, stresses, and displacements.

**Structural analysis and sensitivities:**

```sh
calTop.exe <filename_without_extension> -p 2
```

Here, `-p` sets the density penalization parameter. This mode evaluates the response, adjoint sensitivities, and filtering. Its outputs include:

| File | Contents |
| --- | --- |
| `stress_sens.csv` | Aggregatedn p-norm stress sensitivities |
| `compliance_sens.csv` | Compliance sensitivities |
| `volume_sens.csv` | Volume fraction sensitivities |
| `center_of_gravity_sens.csv` | Element center-of-gravity sensitivities |
| `rhos.dat` | Filtered element densities |
| `objectives.csv` | Compliance, volume fraction, center-of-gravity and aggregated p-norm stress values |

## Coupled aeroelastic workflow

CalFSI links the structural analysis in CalTop with the aerodynamic analysis in SU2. [preCICE](https://precice.org/) coordinates the exchange of interface forces and displacements so that the fluid and structural solutions can be coupled. CalADJ provides the corresponding coupled aeroelastic derivatives used by gradient-based design optimization. The SU2 coupling requires a separately configured SU2 installation and preCICE setup.

## Optimization with FADO and IPOPT

CalTop can be used with [FADO_pyoptsparse](https://github.com/WabalabaKing/FADO_pyoptsparse) to drive an optimization. That fork provides an IPOPT driver and a pyOptSparse interface; it is a separate project. Install the [`ipyopt` package](https://pypi.org/project/ipyopt/) when using the IPOPT interface.

## Authors

- **Prateek Ranjan** — Project lead; CalTop development. Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
- **Wanzheng Zheng** — CalTop development. Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
- **Ghanendra Das** — CalTop development. School of Aerospace Engineering, Georgia Institute of Technology.

## Research Direction and Support

- **Professor Phillip J. Ansell** — Research direction and financial support. Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
- **Professor Kai A. James** — Research direction and financial support. School of Aerospace Engineering, Georgia Institute of Technology.

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

For inquiries, contact Prateek Ranjan at `prateekr@mit.edu`.

## Acknowledgment

This work was supported by NASA under award 80NSSC19M0125 as part of the Center for High-Efficiency Electrical Technologies for Aircraft (CHEETA).

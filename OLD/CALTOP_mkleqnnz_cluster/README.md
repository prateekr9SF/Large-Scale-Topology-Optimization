# CalTop: A Topology Optimization Framework

## Overview
CalTop is a high-performance topology optimization framework built upon **CalculiX 2.15**, featuring advanced solvers such as **SPOOLES, PARDISO, SUPERLU, and PASTIX**. It implements a **density-based finite element framework** interfaced with **IPOPT** for gradient-based optimization. Additionally, CalTop integrates with **SU2_CFD** to solve **static aero-elastic topology optimization problems**.

## Features
- **Finite Element Analysis (FEA):** Uses CalculiX for structural analysis.
- **Fast Solvers:** Supports SPOOLES, PARDISO, SUPERLU, and PASTIX for efficient computations.
- **Density-Based Topology Optimization:** Elements parameterized with densities.
- **IPOPT Integration:** Utilizes the Interior Point OPTimizer for optimization.
- **SU2_CFD Integration:** Solves static aero-elastic topology optimization problems.

## Installation
### Prerequisites
Ensure that the following dependencies are installed:
- **CalculiX 2.15**
- **IPOPT**
- **SU2_CFD**
- **C and Fortran Compilers** (e.g., `gcc`, `gfortran`)
- **BLAS and LAPACK** (for numerical computations)
- **CMake** (if required by dependencies)

### Build Instructions
1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/CalTop.git
   cd CalTop/src
   ```
2. Set up the build environment:
   ```sh
   cp Makefile.inc.example Makefile.inc
   ```
   Edit `Makefile.inc` to specify solver paths and compiler options.

3. Compile the source code:
   ```sh
   make
   ```
4. (Optional) Install the binaries:
   ```sh
   make install
   ```
   The installation path can be specified in `Makefile.inc`.

## Usage
To run a topology optimization problem:
```sh
./CalTop input_file.inp
```
Example input files can be found in the `examples/` directory.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue.

## License
This project is licensed under [MIT License](LICENSE).

## Contact
For inquiries, please reach out to **Prateek Ranjan** at `prateekr@mit.edu`.

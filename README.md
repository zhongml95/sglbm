# sglbm

Stochastic Galerkin Lattice Boltzmann Method (SG-LBM) — a C++20 implementation for intrusive uncertainty quantification in Lattice Boltzmann simulations.  
This repository provides a research-grade framework to reproduce the algorithm described in:

> **M. Zhong, S. Simonis, M. Krause, M. Frank (2024).**  
> *A Stochastic Galerkin Lattice Boltzmann Method for Uncertainty Quantification.*  
> Journal of Computational Physics, 113344.  
> [https://doi.org/10.1016/j.jcp.2024.113344](https://doi.org/10.1016/j.jcp.2024.113344)

---

## 📖 Overview

The **Stochastic Galerkin Lattice Boltzmann Method (SG-LBM)** extends the classical LBM by embedding uncertainty directly into the solver.  
Key features:

- **Intrusive polynomial chaos expansion (gPC)** within LBM collision and streaming operators  
- **Hermite & Legendre polynomial bases** with configurable order  
- **Smolyak sparse grids & tensor quadrature rules**  
- **High-performance C++20 implementation** with MPI parallelization  
- Compatible with **OpenLB** data structures and workflows  

Applications include:
- Fluid dynamics with uncertain boundary conditions or parameters  
- Propagation of stochastic inflow conditions in urban/porous flows  
- Benchmarking against non-intrusive UQ methods (MCS, SC-gPC, etc.)

---

## ⚙️ Installation

### Dependencies
- C++20 compiler (tested with GCC ≥ 11, Clang ≥ 14)  
- [OpenMPI](https://www.open-mpi.org/) (for distributed parallel runs)  
- [OpenLB](http://www.openlb.net/) (linked as submodule or external dependency)  
- CMake ≥ 3.20  

### Build
```bash
git clone https://github.com/zhongml95/sglbm.git
cd sglbm
mkdir build && cd build
cmake ..
make -j

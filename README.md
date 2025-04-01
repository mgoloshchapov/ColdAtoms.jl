![logo](https://github.com/mgoloshchapov/ColdAtoms.jl/blob/main/docs/assets/logo.png)

# ColdAtoms.jl

[![Build Status](https://github.com/mgoloshchapov/ColdAtoms.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mgoloshchapov/ColdAtoms.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://mgoloshchapov.github.io/ColdAtoms.jl/dev/)


This package provides a set of tools to simulate experiments with neutral atoms. 

## Installation

Paste the following line into the Julia REPL:
```
]add "https://github.com/mgoloshchapov/ColdAtoms.jl"
```
or
```julia
using Pkg; Pkg.add(PackageSpec(url="https://github.com/mgoloshchapov/ColdAtoms.jl"))
```


## Package features
- Simulation of two-photon Rydberg excitation with different sources of decoherence: atom motion, laser phase noise, intermediate state decay

- Simulation of release and recapture experiment
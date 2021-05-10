# ElectronicStructure

[![Build Status](https://github.com/jlapeyre/ElectronicStructure.jl/workflows/CI/badge.svg)](https://github.com/jlapeyre/ElectronicStructure.jl/actions)
[![Coverage](https://codecov.io/gh/jlapeyre/ElectronicStructure.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlapeyre/ElectronicStructure.jl)

This package provides an interface to electronic structure calculation packages.
The only package supported is `pyscf`.

Note that doc strings for Python objects are available at the Julia REPL. For example

    help?> ElectronicStructure.pyscf.gto.Mole

Data structures compatible with Qiskit and OpenFermion are supported.

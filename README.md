# FMM Library

Fast multiplication and matrix library

## Info

This library is a collection of fast multiplication algorithms and matrix operations.
It is written in C++ and is header only. It is designed to be used in other projects
and rewritten to fit the needs of the project.

The library is mostly written for educational purposes rather
than anything else. I highly recommend not using it in any production code.

A test program is included. It can be compiled and executed using the provided makefile.
The test program is not meant to be used as a benchmark, but rather as a proof of concept.

### Important Notes

The library only supports square matrices (for now).

Matrices are stored in a 1D array, which is then mapped to 2D indices using the given () operator.
(See the Matrix class for more information) Obviously, this can be changed to best fit your needs.

There might (will) be bugs in the code and a lot of space for optimisation.

## Algorithms (fmm)

- Naive multiplication
- Strassen's multiplication
- Schonhage Laser Method multiplication (TODO)

## Matrix Operations (SMatrix)

- Addition/Subtraction
- Multiplication by a scalar
- Transposition
- Determinant
- Gaussian Elimination (Row Echelon Form)
- Adjunct
- Rank
- split (Creates 4 submatrices from any square matrix)
- merge (Merges 4 submatrices into a single matrix)
- Reduction (remove \[i]\[j] row and column)
- Identity (n-th order)
- Inverse
- Base (TODO)
- Tensor Product (TODO)
- Kronecker Product (TODO)
- Some other stuff... (read the code)

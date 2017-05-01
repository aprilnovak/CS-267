# hacky Makefile for OpenMP version

ftn -c -O3 -qopenmp solvers.f90 element_matrices.f90 csr_storage.f90 mesh.f90 quadrature.f90 read_data.f90
ftn *.o -O3 -qopenmp fmain.f90 -o openmp.out

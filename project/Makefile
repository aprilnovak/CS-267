MODS = csr_storage.f90 element_matrices.f90 mesh.f90 quadrature.f90 read_data.f90 solvers.f90
OBJS = csr_storage.o element_matrices.o mesh.o quadrature.o read_data.o solvers.o
mpi.out : cmain.f90 $(OBJS)
	ftn -O3 $(OBJS) cmain.f90 -o mpi.out

%.o: %.mod

csr_storage.o : csr_storage.f90 element_matrices.f90 mesh.f90
	ftn -c -O3 csr_storage.f90

element_matrices.o : element_matrices.f90 mesh.f90 quadrature.f90 read_data.f90
	ftn -c -O3 element_matrices.f90

quadrature.o : quadrature.f90 read_data.f90
	ftn -c -O3 quadrature.f90

read_data.o : read_data.f90
	ftn -c -O3 read_data.f90

solvers.o : solvers.f90 csr_storage.f90 mesh.f90
	ftn -c -O3 solvers.f90





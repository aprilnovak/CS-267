CFLAGS = -O3

OBJS = csr_storage.o element_matrices.o mesh.o quadrature.o read_data.o solvers.o

mpi.out : cmain.f90 $(OBJS)
	ftn $(CFLAGS) $(OBJS) cmain.f90 -o mpi.out
serial.out : serial.f90 $(OBJS)
	ftn $(CFLAGS) $(OBJS) serial.f90 -o serial.out

%.o: %.mod

read_data.o : read_data.f90
	ftn -c $(CFLAGS) read_data.f90

quadrature.o : quadrature.f90 read_data.o
	ftn -c $(CFLAGS) quadrature.f90

mesh.o : mesh.f90 read_data.o
	ftn -c $(CFLAGS) mesh.f90

element_matrices.o : element_matrices.f90 mesh.o quadrature.o read_data.o
	ftn -c $(CFLAGS) element_matrices.f90

csr_storage.o : csr_storage.f90 element_matrices.o mesh.o
	ftn -c $(CFLAGS) csr_storage.f90

solvers.o : solvers.f90 csr_storage.o mesh.o
	ftn -c $(CFLAGS) solvers.f90

.PHONY: clean

clean:
	rm *.o *.mod mpi.out serial.out

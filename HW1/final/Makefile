# on Edision we will benchmark you against the default vendor-tuned BLAS. The compiler wrappers handle all the linking. If you wish to compare with other BLAS implementations, check the NERSC documentation.
# This makefile is intended for the GNU C compiler. To change compilers, you need to type something like: "module swap PrgEnv-pgi PrgEnv-gnu" See the NERSC documentation for available compilers.

CC = cc
OPT = -O3 -march=native -ffast-math -ffp-contract=fast -fno-signed-zeros -freciprocal-math -funroll-loops
CFLAGS = -Wall -std=gnu99 $(OPT) #-O1 -g
LDFLAGS = -Wall
# librt is needed for clock_gettime
LDLIBS = -lrt

targets = benchmark-blocked
objects = benchmark.o dgemm-blocked.o

.PHONY : default
default : all

.PHONY : all
#all : $(targets)
all : clean $(targets)

benchmark-blocked : benchmark.o dgemm-blocked.o
	$(CC) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects)

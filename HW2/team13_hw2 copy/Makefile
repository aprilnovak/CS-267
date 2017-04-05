#
# Edison - NERSC
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = CC
MPCC = CC
CFLAGS = -DTWO_COORD -O3
CFLAGS += -DNDEBUG        # Turn off assertions
#CFLAGS += -DSTATS_ENABLE  # Enable statistics
#CFLAGS += -DRAND_SEED=1  # Force deterministic runs
#CFLAGS += -O0 -g         # Debugging
LIBS =
OPENMP = -fopenmp
# Use -openmp for Intel compilers.  But don't, because it's deprecated

ifneq ($(NERSC_HOST),edison)
  CC = g++
  MPCC = g++
endif


TARGETS = openmp autograder serial

all:	$(TARGETS)

serial: serial.o xcommon.o
	$(CC) -o $@ $(LIBS) serial.o xcommon.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o
openmp: openmp.o xcommon.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o xcommon.o

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
openmp.o: openmp.cpp xcommon.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp xcommon.h
	$(CC) -c $(CFLAGS) serial.cpp
xcommon.o: xcommon.cpp xcommon.h
	$(CC) -c $(CFLAGS) xcommon.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt

#
# Hopper
#
CC = CC
MPCC = CC
OPENMP =  -mp
CFLAGS = -O3
LIBS =

#
# local (gcc)
#
#CC = g++
#MPCC = mpic++
#OPENMP = -fopenmp
#CFLAGS = -O3
#LIBS =


TARGETS = serial openmp mpi pthreads

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ -pg $(LIBS) serial.o common.o
pthreads: pthreads.o common.o
	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o
openmp: openmp.o common.o
	$(CC) -o $@ -pg $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o

openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

package:
	rm -f team15_hw2.tar*
	tar -cf team15_hw2.tar *.cpp *.h *.pdf Makefile
	tar -tvf team15_hw2.tar
	gzip team15_hw2.tar

clean:
	rm -f *.o $(TARGETS) team15_hw2.tar*

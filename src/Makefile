# MAC
CBLASDIR = /opt/local
CC = gcc-mp-4.8 -m64
CXX = g++-mp-4.8 -m64
FC = gfortran-mp-4.7 -m64
CFLAGS += -I. -O3 -I$(CBLASDIR)/include/
CXXFLAGS += -I. -O3 -I$(CBLASDIR)/include/
FFLAGS += -w -O3
LAPACKLIB = -L$(CBLASDIR)/lib -llapack -latlas -lcblas -lf77blas

# NEMO
#CC = gcc -DPARALLEL
#CXX = g++ -DPARALLEL
#FC = gfortran -DMPI
#CFLAGS += -I. -g -O3 -I/home/pgraff/include -I/usr/include -I/usr/include/mpi
#CXXFLAGS += -I. -g -O3 -I/home/pgraff/include -I/usr/include -I/usr/include/mpi
#FFLAGS += -w -O3 -g
#LAPACKLIB = -L/usr/lib -L/home/pgraff/lib -lcblas -llapack -lmpi -lmpi_cxx -lmpi_f77

LIBS =  -L/home/pgraff/BAMBI -lbambi-1.2 $(LAPACKLIB) -lstdc++ -Lnumerical_recipes_double -lnr -lm -lgsl

CFLAGS += -I/home/pgraff/BAMBI/src -Inumerical_recipes_double

export CC CXX FC CFLAGS CXXFLAGS FFLAGS
 
OBJFILES = bambi.o utils.o main.o

default: all

all: main test genpop

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c

%.o: %.cc
	$(CXX) $(CFLAGS) -c $*.cc

test: kstest.c test.c
	$(CC) -O3 -o test kstest.c test.c 

nr:
	make -C numerical_recipes_double

main: nr $(OBJFILES)
	$(FC) $(FFLAGS) -o ../Analysis $(OBJFILES) $(LIBS)

genpop: nr
	$(CXX) $(CFLAGS) -o GenPop genPopTest.cc -Lnumerical_recipes_double -lnr -L/usr/lib -lgsl -lgslcblas

clean:
	rm -f *.o test ../Analysis GenPop

cleanall: clean
	make -C numerical_recipes_double clean

tidy:
	rm -f *.o

tidyall: tidy
	make -C numerical_recipes_double tidy

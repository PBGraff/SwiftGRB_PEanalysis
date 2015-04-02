CBLASDIR = /opt/local

CC = gcc-mp-4.8 -m64
CXX = g++-mp-4.8 -m64
FC = gfortran-mp-4.7 -m64
CFLAGS += -I. -O3 -I$(CBLASDIR)/include/
CXXFLAGS += -I. -O3 -I$(CBLASDIR)/include/
FFLAGS += -w -O3
LAPACKLIB = -L$(CBLASDIR)/lib -llapack -latlas -lcblas -lf77blas

LIBS =  -L/usr/local/lib -lbambi-1.2 $(LAPACKLIB) -lstdc++ -Lnumerical_recipes_double -lnr -lm -lgsl

CFLAGS += -I/usr/local/include -Inumerical_recipes_double

#CFLAGS += -g -pg
#FFLAGS += -g -pg

export CC CXX FC CFLAGS CXXFLAGS FFLAGS
 
OBJFILES = bambi.o utils.o main.o

default: main

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
	$(CXX) $(CFLAGS) -o GenPop genPopTest.cc -Lnumerical_recipes_double -lnr -L/opt/local/lib -lgsl

clean:
	rm -f *.o test ../Analysis GenPop
	make -C numerical_recipes_double clean

tidy:
	rm -f *.o
	make -C numerical_recipes_double tidy

# MAC
MCMC = $(PWD)/ptmcmc
CBLASDIR = /opt/local
BAMBIDIR = ~/Projects/GWDA/LISA-type-response/flare/BAMBI
INCLIB = $(CBLASDIR) $(BAMBIDIR) $(MCMC) /usr/local
INCS=$(INCLIB:%=-I %/include) -Inumerical_recipes_double
CC = x86_64-apple-darwin15-gcc-mp-5
CXX = x86_64-apple-darwin15-g++-mp-5
FC =  x86_64-apple-darwin15-gfortran-mp-5
CFLAGS += -I. -O1 -g $(INCS)  $(INCLIB:%=-L %/lib)
CXXFLAGS += -I. -O1 -g -std=c++11 $(INCS)
FFLAGS += -w -O3 -g $(INCLIB:%=-L %/lib)
LAPACKLIB = -L$(CBLASDIR)/lib -llapack -latlas -lcblas -lf77blas 
LIBS =  -L/usr/local/lib -lbambi-1.2 $(LAPACKLIB) -lstdc++ -Lnumerical_recipes_double -lnr -lm -lgsl
LIBS2 = -L/usr/local/lib -lbambi-1.2 -L$(CBLASDIR)/lib -lstdc++ -Lnumerical_recipes_double -lnr -lm -lgsl 

# NEMO
#CC = gcc -DPARALLEL
#CXX = g++ -DPARALLEL
#FC = gfortran -DMPI
#CFLAGS += -I. -g -O3 -I/home/pgraff/include -I/usr/include -I/usr/include/mpi
#CXXFLAGS += -I. -g -O3 -I/home/pgraff/include -I/usr/include -I/usr/include/mpi
#FFLAGS += -w -O3 -g
#LAPACKLIB = -L/usr/lib -L/home/pgraff/lib -lcblas -llapack -lmpi -lmpi_cxx -lmpi_f77
#LIBS =  -L/home/pgraff/BAMBI -lbambi-1.2 $(LAPACKLIB) -lstdc++ -Lnumerical_recipes_double -lnr -lm -lgsl
#LIBS2 = $(LIBS)
#CFLAGS += -I/home/pgraff/BAMBI/src -Inumerical_recipes_double


export CC CXX FC CFLAGS CXXFLAGS FFLAGS
 
OBJFILES = bambi.o utils.o main.o
OBJFILESM = utils.o mcmc.o
OBJFILES2 = utils.o detection_fraction.o
OBJFILES3 = utils.o detection_fraction_RF.o
OBJFILES3B = utils.o detection_fraction_AB.o
OBJFILES3C = utils.o detection_fraction_flux.o
OBJFILES4 = utils.o detection_fraction_both.o

default: all

#all: main detfrac rfdetfrac abdetfrac detfracflux detfracboth
#all: main detfrac rfdetfrac abdetfrac detfracflux mcmc
all: main detfrac rfdetfrac abdetfrac detfracflux

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $*.cc

test: kstest.c test.c
	$(CC) -O3 -o test kstest.c test.c 

nr:
	make -C numerical_recipes_double

main: nr $(OBJFILES)
	$(FC) $(FFLAGS) -o ../Analysis $(OBJFILES) $(LIBS)

detfrac: nr $(OBJFILES2)
	$(CXX) $(CFLAGS) -o ../DetectionFraction $(OBJFILES2) $(LIBS2)

rfdetfrac: nr $(OBJFILES3)
	$(CXX) $(CFLAGS) -o ../RFDetectionFraction $(OBJFILES3) $(LIBS2)

abdetfrac: nr $(OBJFILES3B)
	$(CXX) $(CFLAGS) -o ../ABDetectionFraction $(OBJFILES3B) $(LIBS2)

detfracflux: nr $(OBJFILES3C)
	$(CXX) $(CFLAGS) -o ../DetectionFractionFlux $(OBJFILES3C) $(LIBS2)

detfracboth: nr $(OBJFILES4)
	$(CXX) $(CFLAGS) -o ../DetectionFractionBoth $(OBJFILES4) $(LIBS2)

genpop: nr
	$(CXX) $(CFLAGS) -o GenPop genPopTest.cc -Lnumerical_recipes_double -lnr -L/usr/lib -lgsl -lgslcblas

mcmc.o:GRB_bayesian_components.hh

mcmc: nr  ${MCMC}/lib/libptmcmc.a $(OBJFILESM)
	$(FC) $(FFLAGS) -fopenmp -o ../MCMC_Analysis $(OBJFILESM) $(LIBS) -lptmcmc -lprobdist

#cleaning

clean:
	rm -f *.o test ../Analysis ../DetectionFraction GenPop
	${MAKE} -C ptmcmc clean
	rm -f lib/*.a
	rm -f include/*.h*

cleanall: clean
	make -C numerical_recipes_double clean

tidy:
	rm -f *.o

tidyall: tidy
	make -C numerical_recipes_double tidy

#Stuff to build ptmcmc

.ptmcmc-version: ${MCMC}/lib/libptmcmc.a ${MCMC}/lib/libprobdist.a
	cd ptmcmc;git rev-parse HEAD > ../.ptmcmc-version;git status >> ../.ptmcmc-version;git diff >> ../.ptmcmc-version

ptmcmc:
	@echo "Need to check out ptmcmc from github:"
	@git clone https://github.com/JohnGBaker/ptmcmc.git

$(MCMC)/lib/libptmcmc.a: ptmcmc
	@echo "Descending to ptmcmc"
	${MAKE} -C ${MCMC}  ${MFLAGS} ${MCMC}/lib/libptmcmc.a

$(MCMC)/lib/libprobdist.a: ${LIB} ${INCLUDE} ptmcmc
	@echo "Descending into ptmcmc"
	${MAKE} -C ${MCMC}  ${MFLAGS} ${MCMC}/lib/libprobdist.a

clean: tidy

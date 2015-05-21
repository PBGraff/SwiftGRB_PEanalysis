#ifndef __UTILS_H__
#define __UTILS_H__ 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>

typedef enum {NEURALNET, RANDOMFOREST, ADABOOST, FLUXTHRESH} MLmethod;

typedef struct RunOptions {
	double n0;
	double n1;
	double n2;
	double z1;
	long int seed;
	int resume;
	long int popsize;
	long int datapopsize;
	int help;
	char datafile[200];
	int nlive;
	bool nstar;
	bool ntotal;
	bool flatn0;
	double tobs;
	bool testpop;
	int nbins;
	bool flatzpop;
	double zbin_max;
	double zbin_min;
	int zpts;
	bool zeroLogLike;
	char outfile[200];
	int verbose;
	MLmethod method;
	bool vary_z1;
} RunArgs;

// return an array of detected z from input array of z and detection prob, given a prob threshold
void detected(double *trigz, double *ptrig, long int ntrig, double pth, double *detz, long int *ndet);

// count lines in a file
long int countlines(char filename[]);

// read command-line options
void read_options(int argc, char **argv, RunArgs *args);

// log of Poisson probability
double logPoisson(double k, double lambda);
double logPoisson2(double k, double lambda);
double logPoisson3(double k, double lambda);

// bin detected GRB redshifts
void bindetections(double *zdet, long int ndet, double zmin, double zmax, int nbins, int *count);

double logFactorial(double x);
double factorial(double x);

#endif

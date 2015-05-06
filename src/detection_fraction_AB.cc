#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include "NeuralNetwork.h"

extern "C" {
	#include "utils.h"
	#include "mock_sample_functions.h"
}

//#define POPSIZE		5000
//#define DPOPSIZE	100

#define NINPUTS		15 			// do not change this!

#define Z1DATA			3.60
#define XDATA			-0.65
#define YDATA			-3.00
#define LOGLSTARDATA	52.05

#define ZMIN		0.0
#define ZMAX		10.0

double *zdata=NULL;
long int ndetdata=0, ndetpop;
double *population=NULL, *zpop=NULL, *ppop=NULL, *zdetpop=NULL;
float *sample=NULL;
double ksd, ksp, logpois, logpois0 = 0.0;
float prob[2];
RunArgs runargs;
int *detcount=NULL, *popcount=NULL;
double *detfrac=NULL;

extern "C" {
	#include "mock_sample_functions.c"
	#include "population_fixedz.c"
}

int main(int argc, char *argv[])
{
#ifdef PARALLEL
 	MPI_Init(&argc,&argv);
#endif

	int myid = 0, ncpus = 1;
#ifdef PARALLEL
 	MPI_Comm_size(MPI_COMM_WORLD,&ncpus);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

	int i,j,k;

	// initialize default options for arguments
	runargs.resume = 0;
	runargs.help = 0;
	runargs.n0 = 0.84;
	runargs.n1 = 2.07;
	runargs.n2 = -0.7;
	runargs.popsize = 1000;
	runargs.datapopsize = 10000;
	runargs.seed = 0;
	strcpy(runargs.datafile,"\0");
	runargs.nlive = 100;
	runargs.nstar = false;
	runargs.flatn0 = false;
	runargs.tobs = 1.0;
	runargs.nbins = 50;
	runargs.zpts = 1001;

	long int dataseed=0;

	// get command-line options
	read_options(argc, argv, &runargs);

	if ( runargs.help == 1 )
	{
		char helpstr[] = "\n\
This program runs BAMBI on the Swift GRB population-fitting problem. Here are available options:\n\
\n\
--help       Print this help and exit\n\
\n\
Data Settings\n\
-------------------------------------------------------------------------------------------------\n\
--seed       Initial seed for generating simulated data, (default=0 uses sys time)\n\
--dpop       Population size for simulated data (default=10000)\n\
--zpts       Number of points in z space to evaluate (default=1001 => dz=0.01)\n\
\n";
		printf("%s",helpstr);
		return 0;
	}

	if (runargs.seed == 0)
	{
		dataseed = -34712 * (long int) time(NULL);
	} else {
		dataseed = -34712 * runargs.seed;	
	}

	// read in values for the background possibilities of simulated GRBs
	read_background_values();

	// load the splines
	load_splines();

	// allocate memory
	sample = (float *) malloc(NINPUTS * sizeof(float));

	double *datapop=NULL, *dataz=NULL, *dataprob=NULL;
	datapop = (double *) malloc(runargs.datapopsize * NINPUTS * sizeof(double));
	dataz = (double *) malloc(runargs.datapopsize * sizeof(double));
	dataprob = (double *) malloc(runargs.datapopsize * sizeof(double));
	zdata = (double *) malloc(runargs.datapopsize * sizeof(double));
	float dprob[2];

	double dz = (ZMAX - ZMIN) / (double) (runargs.zpts - 1);
	double *z = (double *) malloc(runargs.zpts * sizeof(double));
	double *detfrac = (double *) malloc(runargs.zpts * sizeof(double));

	for (i = 0; i < runargs.zpts; i++)
	{
		z[i] = ZMIN + (double) i * dz;
		detfrac[i] = 0.0;
	}

	int Nz = (int) ceil((double)runargs.zpts / (double) ncpus);
	int zstart = myid * Nz;
	int zend = (int) fmin((double)(myid+1)*Nz, (double)runargs.zpts);
	printf("ID/Start/End = %d/%d/%d\n", myid, zstart, zend-1);

	for (k = zstart; k < zend; k++)
	{
		// make file names
		char outfilename[100], infilename[100];
		sprintf(outfilename, "population_data_%d.txt", k);
		sprintf(infilename, "population_AB_predictions_%d.txt", k);

		// simulate population
		GeneratePopulationFixZ(datapop, runargs.datapopsize, z[k], XDATA, YDATA, LOGLSTARDATA, dataz, &dataseed);

		// write to file
		FILE *outfileptr = fopen(outfilename, "w");
		for ( i=0; i<runargs.datapopsize; i++)
		{
			for ( j=0; j<NINPUTS; j++ )
			{
				fprintf(outfileptr, "%lf ", datapop[i*NINPUTS+j]);
			}
			fprintf(outfileptr, "\n");
		}
		fclose(outfileptr);

		// run python script for RF
		char command[200];
		sprintf(command, "python evalAB.py %s %s", outfilename, infilename);
		system(command);

		// collect results
		FILE *infileptr = fopen(infilename, "r");
		for ( i=0; i<runargs.datapopsize; i++)
		{
			fscanf(infileptr, "%lf\n", &dataprob[i]);
		}
		fclose(infileptr);

		// find detected GRBs
		detected(dataz, dataprob, runargs.datapopsize, 0.5, zdata, &ndetdata);

		detfrac[k] = (double) ndetdata / (double) runargs.datapopsize;
		
		printf("Simulated %ld GRBs (%ld detected) at z=%g : %g%%\n", runargs.datapopsize, ndetdata, z[k], detfrac[k]*100.0);

		// clear data files
		sprintf(command, "rm -f %s %s", outfilename, infilename);
		system(command);	
	}

	// collect all results at root node
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid != 0)
	{
		MPI_Send(&detfrac[zstart], zend - zstart, MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
	}
	else
	{
		for (i = 1; i<ncpus; i++)
		{
			int zstart_tmp = i * Nz;
			int zend_tmp = (int) fmin((double)(i+1)*Nz, (double)runargs.zpts);
			MPI_Status status;
			MPI_Recv(&detfrac[zstart_tmp], zend_tmp - zstart_tmp, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
		}
	}
#endif
	
	// write results to file
	if (myid == 0)
	{
		// open file for outputs
		FILE *outfile = fopen("support_data/splines_detection_fraction_z_AB.txt", "w");

		// write results
		for (i = 0; i < runargs.zpts; i++)
		{
			fprintf(outfile, "%lf %lf\n", z[i], detfrac[i]);
		}

		// close the file
		fclose(outfile);	
	}
	
	// clean up allocated variables
	free(z);
	free(detfrac);
	free(datapop);
	free(dataz);
	free(dataprob);
	free(zdata);
	free(sample);
	unload_splines();

#ifdef PARALLEL
 	MPI_Finalize();
#endif

	return 0;
}


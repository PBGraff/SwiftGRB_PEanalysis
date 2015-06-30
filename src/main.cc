#include "main.h"
#include "NeuralNetwork.h"
extern "C" {
	#include "utils.h"
	#include "mock_sample_functions.h"
}
#include <time.h>

//#define POPSIZE		5000
//#define DPOPSIZE	100

#define NINPUTS		15 			// do not change this!

#define Z1DATA			3.60
#define XDATA			-0.65
#define YDATA			-3.00
#define LOGLSTARDATA	52.05

#define ZMIN		0.0
#define ZMAX		10.0

#define LOG10_FLUX_THRESHOLD	-7.243

NeuralNetwork *GRBnn;
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
	#include "population.c"
}

double CubeToFlatPrior(double r, double xmin, double xmax)
{
	return r*(xmax-xmin)+xmin;
}

double CubeToLogPrior(double r, double xmin, double xmax)
{
	double lmin=log(xmin),lmax=log(xmax);
	return exp(r*(lmax-lmin)+lmin);
}

/******************************************** getphysparams routine ****************************************************/

void getphysparams(double *Cube, int &ndim, int &nPar, void *context)
{
	double n0, n1, n2, nstar, ntotal, z1;

	// n1 from Cube[1]
	n1 = CubeToFlatPrior(Cube[1], 0.00, 4.00);

	// n2 from Cube[2]
	n2 = CubeToFlatPrior(Cube[2], -6.00, 0.00);

	// z1 from Cube[3]
	if (runargs.vary_z1) {
		z1 = CubeToFlatPrior(Cube[3], 0.00, 10.0);
	} else {
		z1 = Z1DATA;
	}

	// normalization parameters from Cube[0]
	if (runargs.nstar) {
		// nstar
		nstar = CubeToLogPrior(Cube[0], 0.10, 10000.0);
		// n0
		n0 = nstar * pow(1.0 + z1, -n1);
		// ntotal
		ntotal = GRBNumberIntegral(n0, n1, n2, z1);
	} else if (runargs.ntotal) {
		// ntotal
		ntotal = CubeToLogPrior(Cube[0], 1.00, 1e5);
		double ntmp = GRBNumberIntegral(1.0, n1, n2, z1);
		// n0
		n0 = ntotal / ntmp;
		// nstar
		nstar = n0 * pow(1.0 + z1, n1);
	} else {
		// n0
		if (runargs.flatn0) {
			n0 = CubeToFlatPrior(Cube[0], 0.01, 2.00);
		} else {
			n0 = CubeToLogPrior(Cube[0], 0.01, 2.00);
		}
		// nstar
		nstar = n0 * pow(1.0 + z1, n1);
		// ntotal
		ntotal = GRBNumberIntegral(n0, n1, n2, z1);
	}

	Cube[0] = n0;
	Cube[1] = n1;
	Cube[2] = n2;
	Cube[3] = z1;
	Cube[4] = nstar;
	Cube[5] = ntotal;
}

/******************************************** getallparams routine ****************************************************/

void getallparams(double *Cube, int &ndim, int &nPar, void *context)
{
	getphysparams(Cube,ndim,nPar,context);

	// x
	Cube[6] = XDATA;

	// y
	Cube[7] = YDATA;

	// log10(L_star)
	Cube[8] = LOGLSTARDATA;
}

/******************************************** loglikelihood routine ****************************************************/

// Now an example, sample an egg box likelihood

// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//	 
// Output arguments
// lnew 						= loglikelihood

void getLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
	//clock_t t = clock();

	// Extract population parameters
	getallparams(Cube,ndim,npars,context);
	double n0, n1, n2, z1;
	n0 = Cube[0];
	n1 = Cube[1];
	n2 = Cube[2];
	z1 = Cube[3];

	// compute logL as in notes
	if (runargs.zeroLogLike)
	{
		lnew = 0.0;
	}
	else
	{
		lnew = -1.0 * GRBRateIntegral(n0, n1, n2, z1);
		int i;
		for (i = 0; i < ndetdata; i++)
		{
			lnew += log(GRBRate(zdata[i], n0, n1, n2, z1));
		}
		//printf("logL = %lf\n", lnew);
	}

	//t = clock() - t;
	//printf("%lf ms for logL\n", 1.0e3 * (double)t / CLOCKS_PER_SEC);
}


/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void *context)
{
	/*
	// convert the 2D Fortran arrays to C arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	int i, j;
	
	double postdist[*nSamples][*nPar + 2];
	for( i = 0; i < *nPar + 2; i++ )
		for( j = 0; j < *nSamples; j++ )
			postdist[j][i] = posterior[0][i * (*nSamples) + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[*nlive][*nPar + 1];
	for( i = 0; i < *nPar + 1; i++ )
		for( j = 0; j < *nlive; j++ )
			pLivePts[j][i] = physLive[0][i * (*nlive) + j];
	*/
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
#ifdef PARALLEL
 	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
	
	int i,j;

	// initialize default options for arguments
	runargs.resume = 0;
	runargs.help = 0;
	runargs.n0 = 0.42;
	runargs.n1 = 2.07;
	runargs.n2 = -0.7;
	runargs.z1 = Z1DATA;
	runargs.popsize = 1000;
	runargs.datapopsize = 100;
	runargs.seed = 0;
	strcpy(runargs.datafile,"\0");
	runargs.nlive = 100;
	runargs.nstar = false;
	runargs.flatn0 = false;
	runargs.tobs = 0.8;
	runargs.nbins = 50;
	runargs.zeroLogLike = false;
	strcpy(runargs.outfile,"");
	runargs.verbose = 1;
	runargs.method = NEURALNET;
	runargs.vary_z1 = false;

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
Run Settings\n\
-------------------------------------------------------------------------------------------------\n\
--resume         Resume the run from before (default=off)\n\
--nlive          Number of live points for BAMBI to use (default=100)\n\
--bins           [deprecated] Number of bins to use in redshift for likelihood (default=50)\n\
--zeroLogLike    Sample from prior with logL=0\n\
--outfile        Root for output files\n\
--silent         Suppress BAMBI updates to stdout\n\
--method         Choice of model to use for Swift detection:\n\
                 0 = Neural Network (default)\n\
                 1 = Random Forest\n\
                 2 = AdaBoost\n\
                 3 = Flux threshold\n\
\n\
Data Settings\n\
-------------------------------------------------------------------------------------------------\n\
--seed           Initial seed for generating simulated data, will use later defined n0, n1, and n2\n\
                 (default=0 reads in data from file)\n\
--dpop           [deprecated] population size for simulated data (optional, default=100)\n\
--n0             n0 for generated simulated data (optional, default=0.42)\n\
--n1             n1 for generated simulated data (optional, default=2.07)\n\
--n2             n2 for generated simulated data (optional, default=-0.7)\n\
--z1             z1 for generated simulated data (optional, default=3.60)\n\
--file           data file to be read in (optional for when seed=0)\n\
\n\
Model Settings\n\
-------------------------------------------------------------------------------------------------\n\
--pop            [deprecated] population size for simulated models (optional, default=1000)\n\
--nstar          use GRB rate at peak instead of rate at z=0\n\
--ntotal         use total intrinsic GRB pop size instead of rate at z=0\n\
--flatn0         use flat prior (instead of log) on n0\n\
--varyz1         allow z1 to vary from fixed 3.6\n\
\n";
		printf("%s",helpstr);
#ifdef PARALLEL
 	MPI_Finalize();
#endif
		return 0;
	}

	if ( runargs.seed==0 && strcmp(runargs.datafile,"")==0 )
	{
		fprintf(stderr, "You need to provide either a data seed (--seed) or an input data file (--file).\n");
#ifdef PARALLEL
 	MPI_Finalize();
#endif
		return 0;
	}

	dataseed = -347*runargs.seed;

	// read in values for the background possibilities of simulated GRBs
	read_background_values();

	// load the splines
	load_splines();

	// allocate memory
	sample = (float *) malloc(NINPUTS * sizeof(float));
	detcount = (int *) malloc(runargs.nbins * sizeof(int));
	popcount = (int *) malloc(runargs.nbins * sizeof(int));

	char outroot[100] = "";
	if ( runargs.seed == 0 )
	{
		if (strlen(runargs.outfile)==0)
		{
			sprintf(outroot, "chains/analysis_realdata_");
		}
		else
		{
			strcpy(outroot, runargs.outfile);
		}
		
		// Read in data
		ndetdata = countlines(runargs.datafile);
		zdata = (double *) malloc(ndetdata * sizeof(double));
		FILE *fptr = fopen(runargs.datafile, "r");
		for ( i=0; i<ndetdata; i++ )
		{
			fscanf(fptr, "%lf\n", &zdata[i]);
		}
		fclose(fptr);
		printf("Data read in from file with %d detected GRBs.\n", (int) ndetdata);
		logpois0 = -0.5 * log(2.0 * M_PI * (double) ndetdata);
	}
	else
	{
		if (strlen(runargs.outfile)==0)
		{
			sprintf(outroot, "chains/analysis_n0%d_n1%d_n2%d_seed%ld_", (int)round(fabs(runargs.n0*100)), 
					(int)round(fabs(runargs.n1*100)), (int)round(fabs(runargs.n2*100)), runargs.seed);
		}
		else
		{
			strcpy(outroot, runargs.outfile);
		}
		
		// simulate data
		// calculate population size
		double all_sky_rate = GRBNumberIntegral(runargs.n0, runargs.n1, runargs.n2, runargs.z1);
		printf("Computed an all-sky intrinsic rate of %lf GRBs/yr.\n", all_sky_rate);
		runargs.datapopsize = (long int) (all_sky_rate * runargs.tobs / 6.0);
		// allocate memory
		double *datapop=NULL, *dataz=NULL, *dataprob=NULL;
		datapop = (double *) malloc(runargs.datapopsize * NINPUTS * sizeof(double));
		dataz = (double *) malloc(runargs.datapopsize * sizeof(double));
		dataprob = (double *) malloc(runargs.datapopsize * sizeof(double));
		zdata = (double *) malloc(runargs.datapopsize * sizeof(double));
		float dprob[2];
		// simulate population
		GeneratePopulation(datapop, runargs.datapopsize, runargs.n0, runargs.n1, runargs.n2, runargs.z1, XDATA, YDATA, LOGLSTARDATA, dataz, &dataseed);
		
		// find detection probabilities
		if (runargs.method == NEURALNET)
		{
			// read in saved NN
			GRBnn = new FeedForwardClassNetwork();
			GRBnn->read("support_data/Swift_NN_all_nhid-100-50_act330_network.txt");
			
			// perform predictions
			for ( i=0; i<runargs.datapopsize; i++)
			{
				for ( j=0; j<NINPUTS; j++ )
				{
					sample[j] = (float) datapop[i*NINPUTS+j];
				}
				GRBnn->forwardOne(1, &sample[0], &dprob[0]);
				dataprob[i] = (double) dprob[1];
			}
		}
		else if (runargs.method == RANDOMFOREST || runargs.method == ADABOOST)
		{
			char outfilename[100], infilename[100], command[200];
			sprintf(outfilename, "population_data.txt");
			sprintf(infilename, "population_predictions.txt");

			if (myid == 0)
			{
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
				if (runargs.method == RANDOMFOREST)
				{
					sprintf(command, "python evalRF.py %s %s", outfilename, infilename);
				}
				else
				{
					sprintf(command, "python evalAB.py %s %s", outfilename, infilename);
				}
				system(command);
			}

#ifdef PARALLEL
 			MPI_Barrier(MPI_COMM_WORLD);
#endif

			FILE *infileptr = fopen(infilename, "r");
			for ( i=0; i<runargs.datapopsize; i++)
			{
				fscanf(infileptr, "%lf\n", &dataprob[i]);
			}
			fclose(infileptr);

#ifdef PARALLEL
 			MPI_Barrier(MPI_COMM_WORLD);
#endif

			if (myid == 0)
			{
				sprintf(command, "rm -f %s %s", outfilename, infilename);
				system(command);
			}
		}
		else if (runargs.method == FLUXTHRESH)
		{
			for ( i=0; i<runargs.datapopsize; i++)
			{
				if (datapop[i*NINPUTS+13] >= LOG10_FLUX_THRESHOLD)
					dataprob[i] = 1.0;
				else
					dataprob[i] = 0.0;
			}
		}
		else
		{
			fprintf(stderr, "Invalid method chosen: %d\n", runargs.method);
#ifdef PARALLEL
 			MPI_Finalize();
#endif
 			exit(-1);
		}

		// find list of detected GRBs
		detected(dataz, dataprob, runargs.datapopsize, 0.5, zdata, &ndetdata);
		printf("Simulated data population generated with %ld GRBs => %ld detected\n", runargs.datapopsize, ndetdata);
		
		// if printing a test population
		if (runargs.testpop)
		{
			FILE *fptr = fopen("population_test.txt","w");
			for ( i=0; i<runargs.datapopsize; i++)
			{
				for ( j=0; j<NINPUTS; j++ )
				{
					fprintf(fptr, "%lf ", datapop[i*NINPUTS+j]);
				}
				fprintf(fptr, "%lf\n", dataprob[i]);
			}
			fclose(fptr);
		}
		
		// free memory
		free(datapop);
		free(dataz);
		free(dataprob);

		// calculate new logL function at true values
		double logLnew = -1.0 * GRBRateIntegral(runargs.n0, runargs.n1, runargs.n2, runargs.z1);
		for (i = 0; i < ndetdata; i++)
		{
			logLnew += log(GRBRate(zdata[i], runargs.n0, runargs.n1, runargs.n2, runargs.z1));
		}
		printf("New logL = %lf at true values\n", logLnew);

		// save detected GRB redshifts
		char datasavefile[200];
		sprintf(datasavefile, "%sdetectedZdata.txt", outroot);
		FILE *datasave = fopen(datasavefile, "w");
		for (i = 0; i < ndetdata; i++)
		{
			fprintf(datasave, "%lf\n", zdata[i]);
		}
		fclose(datasave);
		printf("Detected GRB redshifts saved to %s\n", datasavefile);

		// save injected true values
		sprintf(datasavefile, "%sinjected_values.txt", outroot);
		datasave = fopen(datasavefile, "w");
		fprintf(datasave, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n", runargs.n0, runargs.n1, runargs.n2, runargs.z1,
				runargs.n0 * pow(1.0 + runargs.z1, runargs.n1), all_sky_rate);
		fclose(datasave);
		printf("Injected values saved to %s\n", datasavefile);
	}

	if (runargs.testpop)
	{
#ifdef PARALLEL
 		MPI_Finalize();
#endif
 		exit(0);
	}

	printf("Writing outputs to %s*\n", outroot);

	// set the MultiNest sampling parameters
	
	int mmodal = 0;					// do mode separation?
	
	int ceff = 0;					// run in constant efficiency mode?
	
	int nlive = runargs.nlive;				// number of live points
	
	double efr = 0.1;				// set the required efficiency
	
	double tol = 0.1;				// tol, defines the stopping criteria
	
	int ndims = 3;					// dimensionality (no. of free parameters)
	if (runargs.vary_z1) ndims++;
	
	int nPar = 9;					// total no. of parameters including free & derived parameters
	
	int nClsPar = 3;				// no. of parameters to do mode separation on
	if (runargs.vary_z1) nClsPar++;
	
	int updInt = 50;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	
	int maxModes = 1;				// expected max no. of modes (used only for memory allocation)
	
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndims; i++) pWrap[i] = 0;
	
	strcpy(root, outroot);			// root for output files
	strcpy(networkinputs, "net.inp");			// file with input parameters for network training
	
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	
	int fb = runargs.verbose;					// need feedback on standard output?
	
	resume = runargs.resume;					// resume from a previous job?
	
	int outfile = 1;				// write output files?
	
	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	
	logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass
	
	doBAMBI = 0;					// BAMBI?

	useNN = 0;
	
	// calling MultiNest

	nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, LogLike, dumper, bambi, context);

	// clean up allocated variables
	if (runargs.method == NEURALNET) delete GRBnn;
	free(zdata);
	free(sample);
	free(detcount);
	free(popcount);
	unload_splines();
	
#ifdef PARALLEL
 	MPI_Finalize();
#endif
}

/***********************************************************************************************************************/

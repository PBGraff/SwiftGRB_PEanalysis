#include "main.h"
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

NeuralNetwork *GRBnn;
double *zdata=NULL;
long int ndetdata=0, ndetpop;
double *population=NULL, *zpop=NULL, *ppop=NULL, *zdetpop=NULL;
float *sample=NULL;
//double n0_data, n1_data, n2_data;
double ksd, ksp, logpois, logpois0 = 0.0;
float prob;
//int popsize, dpopsize;
//bool nstar = false;
RunArgs runargs;

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
	// n1
	Cube[1] = CubeToFlatPrior(Cube[1], 1.60, 4.00);

	// n2
	Cube[2] = CubeToFlatPrior(Cube[2], -4.00, 0.00);

	// n0
	if (runargs.nstar) {
		double nstar = CubeToLogPrior(Cube[0], 1.00, 1000.0);
		Cube[0] = nstar * pow(1.0 + Z1DATA, -Cube[1]);
		Cube[7] = nstar;
	} else {
		if (runargs.flatn0) {
			Cube[0] = CubeToFlatPrior(Cube[0], 0.25, 2.00);
		} else {
			Cube[0] = CubeToLogPrior(Cube[0], 0.25, 2.00);
		}
	}
}

/******************************************** getallparams routine ****************************************************/

void getallparams(double *Cube, int &ndim, int &nPar, void *context)
{
	getphysparams(Cube,ndim,nPar,context);

	// z1
	Cube[3] = Z1DATA;

	// x
	Cube[4] = XDATA;

	// y
	Cube[5] = YDATA;

	// log10(L_star)
	Cube[6] = LOGLSTARDATA;
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
	//clock_t timer = clock();

	// Extract population parameters
	getallparams(Cube,ndim,npars,context);
	double n0,n1,n2,z1,x,y,logLstar;
	n0 = Cube[0];
	n1 = Cube[1];
	n2 = Cube[2];
	z1 = Cube[3];
	x = Cube[4];
	y = Cube[5];
	logLstar = Cube[6];

	int i, j;

	//printf("----------------------------------------\n");
	
	// calculate the population size
	runargs.popsize = GRBNumberIntegral(n0, n1, n2);
	//printf("%lf %lf %lf ==> %ld\n", n0, n1, n2, runargs.popsize);
	
	// allocate memory
	population = (double *) malloc(runargs.popsize * NINPUTS * sizeof(double));
	zpop = (double *) malloc(runargs.popsize * sizeof(double));
	ppop = (double *) malloc(runargs.popsize * sizeof(double));
	zdetpop = (double *) malloc(runargs.popsize * sizeof(double));

	// Make sample population
	long int seed = (long int) time(NULL);
	seed *= -34653;
	GeneratePopulation(population, runargs.popsize, n0, n1, n2, z1, x, y, logLstar, zpop, &seed);
	
	// Predict probability of detection for each GRB
	for ( i=0; i<runargs.popsize; i++)
	{
		for ( j=0; j<NINPUTS; j++ )
		{
			sample[j] = (float) population[i*NINPUTS+j];
		}
		GRBnn->forwardOne(1, &sample[0], &prob);
		ppop[i] = (double) prob;
		//printf("%f\n",prob);
	}

	// Find detected GRBs in population
	detected(zpop, ppop, runargs.popsize, 0.5, zdetpop, &ndetpop);

	// Calculate K-S test p-value
	kstwo(zpop-1, ndetpop, zdata-1, ndetdata, &ksd, &ksp);

	// Calculate Poisson probability for number count
	logpois = logPoisson((double) ndetpop, (double) ndetdata) - logpois0;
	//printf("%ld %ld ==> %lf\n", ndetpop, ndetdata, logpois);

	lnew = log(ksp) + logpois;
	//printf("%lf %lf ==> %lf\n", log(ksp), logpois, lnew);

	// free memory
	free(population);
	free(zpop);
	free(ppop);
	free(zdetpop);
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
	runargs.n0 = 0.84;
	runargs.n1 = 2.07;
	runargs.n2 = -0.7;
	runargs.popsize = 1000;
	runargs.datapopsize = 100;
	runargs.seed = 0;
	strcpy(runargs.datafile,"\0");
	runargs.nlive = 100;
	runargs.nstar = false;
	runargs.flatn0 = false;
	runargs.tobs = 1.0;

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
--resume     Resume the run from before (default=off)\n\
--nlive      Number of live points for BAMBI to use (default=100)\n\
\n\
Data Settings\n\
-------------------------------------------------------------------------------------------------\n\
--seed       Initial seed for generating simulated data, will use later defined n0, n1, and n2\n\
             (default=0 reads in data from file)\n\
--dpop       [deprecated] population size for simulated data (optional, default=100)\n\
--n0         n0 for generated simulated data (optional, default=0.84)\n\
--n1         n1 for generated simulated data (optional, default=2.07)\n\
--n2         n2 for generated simulated data (optional, default=-0.7)\n\
--file       data file to be read in (optional for when seed=0)\n\
\n\
Model Settings\n\
-------------------------------------------------------------------------------------------------\n\
--pop        [deprecated] population size for simulated models (optional, default=1000)\n\
--nstar      use GRB rate at peak instead of rate at z=0\n\
--flatn0     use flat prior (instead of log) on n0\n\
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

	// Read in saved neural network for GRB detection predictions
	GRBnn = new FeedForwardClassNetwork();
	GRBnn->read("support_data/Swift_NN_all_nhid-100-50_act330_network.txt");

	// read in values for the background possibilities of simulated GRBs
	read_background_values();

	// load the splines
	load_splines();

	// allocate memory
	sample = (float *) malloc(NINPUTS * sizeof(float));

	char outroot[100] = "";
	if ( runargs.seed == 0 )
	{
		sprintf(outroot, "chains/analysis_realdata_");
		
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
		sprintf(outroot, "chains/analysis_n0%d_n1%d_n2%d_seed%ld_", (int)round(fabs(runargs.n0*100)), 
				(int)round(fabs(runargs.n1*100)), (int)round(fabs(runargs.n2*100)), runargs.seed);
		
		// simulate data
		// calculate population size
		runargs.datapopsize = GRBNumberIntegral(runargs.n0, runargs.n1, runargs.n2);
		printf("Calculated population size of %ld\n", runargs.datapopsize);
		// allocate memory
		double *datapop=NULL, *dataz=NULL, *dataprob=NULL;
		datapop = (double *) malloc(runargs.datapopsize * NINPUTS * sizeof(double));
		dataz = (double *) malloc(runargs.datapopsize * sizeof(double));
		dataprob = (double *) malloc(runargs.datapopsize * sizeof(double));
		zdata = (double *) malloc(runargs.datapopsize * sizeof(double));
		float dprob;
		// simulate population
		GeneratePopulation(datapop, runargs.datapopsize, runargs.n0, runargs.n1, runargs.n2, Z1DATA, XDATA, YDATA, LOGLSTARDATA, dataz, &dataseed);
		//FILE *fptr = fopen("population_test.txt","w");
		for ( i=0; i<runargs.datapopsize; i++)
		{
			for ( j=0; j<NINPUTS; j++ )
			{
				sample[j] = (float) datapop[i*NINPUTS+j];
				//fprintf(fptr,"%f ",sample[j]);
			}
			GRBnn->forwardOne(1, &sample[0], &dprob);
			dataprob[i] = (double) dprob;
			//fprintf(fptr,"%f\n",dprob);
			//printf("%f\n",dprob);
		}
		//fclose(fptr);
		detected(dataz, dataprob, runargs.datapopsize, 0.5, zdata, &ndetdata);
		printf("Simulated data population generated with %d detected GRBs\n", (int) ndetdata);
		logpois0 = -0.5 * log(2.0 * M_PI * (double) ndetdata);
		// free memory
		free(datapop);
		free(dataz);
		free(dataprob);

		// simulate a similar population and evaluate the likelihood
		// calculate the population size
		runargs.popsize = runargs.datapopsize;
		// allocate memory
		population = (double *) malloc(runargs.popsize * NINPUTS * sizeof(double));
		zpop = (double *) malloc(runargs.popsize * sizeof(double));
		ppop = (double *) malloc(runargs.popsize * sizeof(double));
		zdetpop = (double *) malloc(runargs.popsize * sizeof(double));
		// simulate the population
		GeneratePopulation(population, runargs.popsize, runargs.n0, runargs.n1, runargs.n2, Z1DATA, XDATA, YDATA, LOGLSTARDATA, zpop, &dataseed);
		for ( i=0; i<runargs.popsize; i++)
		{
			for ( j=0; j<NINPUTS; j++ )
			{
				sample[j] = (float) population[i*NINPUTS+j];
			}
			GRBnn->forwardOne(1, &sample[0], &prob);
			ppop[i] = (double) prob;
		}
		detected(zpop, ppop, runargs.popsize, 0.5, zdetpop, &ndetpop);
		kstwo(zpop-1, ndetpop, zdata-1, ndetdata, &ksd, &ksp);
		logpois = logPoisson((double) ndetpop, (double) ndetdata) - logpois0;
		printf("Similar distribution has %d detected GRBs and logL = %lf\n", ndetpop, log(ksp)+logpois);
		// free memory
		free(population);
		free(zpop);
		free(ppop);
		free(zdetpop);
	}
	printf("Writing outputs to %s*\n", outroot);

	// set the MultiNest sampling parameters
	
	int mmodal = 0;					// do mode separation?
	
	int ceff = 0;					// run in constant efficiency mode?
	
	int nlive = runargs.nlive;				// number of live points
	
	double efr = 0.1;				// set the required efficiency
	
	double tol = 0.1;				// tol, defines the stopping criteria
	
	int ndims = 3;					// dimensionality (no. of free parameters)
	
	int nPar = 7;					// total no. of parameters including free & derived parameters
	if (runargs.nstar) nPar++;
	
	int nClsPar = 3;				// no. of parameters to do mode separation on
	
	int updInt = 50;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	
	int maxModes = 1;				// expected max no. of modes (used only for memory allocation)
	
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndims; i++) pWrap[i] = 0;
	
	strcpy(root, outroot);			// root for output files
	strcpy(networkinputs, "net.inp");			// file with input parameters for network training
	
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	
	int fb = 1;					// need feedback on standard output?
	
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
	delete GRBnn;
	free(zdata);
	free(sample);
	unload_splines();
	
#ifdef PARALLEL
 	MPI_Finalize();
#endif
}

/***********************************************************************************************************************/

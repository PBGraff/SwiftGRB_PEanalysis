#include "utils.h"

void detected(double *trigz, double *ptrig, long int ntrig, double pth, double *detz, long int *ndet)
{
	long int i;
	*ndet = 0;

	for ( i=0; i<ntrig; i++ )
	{
		if ( ptrig[i] > pth )
		{
			detz[*ndet] = trigz[i];
			(*ndet)++;
		}
	}
}

long int countlines(char filename[])
{
	FILE *ifp = fopen(filename,"r");
	if ( ifp == NULL )
	{
		fprintf(stderr,"No such file (%s) exists!\n");
		return 0;
	}

	long int lines = 0;
	
	char ch;
	while( !feof(ifp) )
	{
		ch = fgetc(ifp);
		if(ch == '\n') lines++;
	}
	
	fclose(ifp);
	
	return lines;
}

void read_options(int argc, char *argv[], RunArgs *args)
{
	int c;
	static int resflag=0, helpflag=0, nstarflag=0, ntotalflag=0, flatn0flag=0, testflag=0, zllflag=0, verboseflag=1;
	static int method=0, varyz1=0;
	
	while (1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"resume", no_argument, &resflag, 1},
			{"help", no_argument, &helpflag, 1},
			{"nstar", no_argument, &nstarflag, 1},
			{"ntotal", no_argument, &ntotalflag, 1},
			{"flatn0", no_argument, &flatn0flag, 1},
			{"test", no_argument, &testflag, 1},
			{"zeroLogLike", no_argument, &zllflag, 1},
			{"silent", no_argument, &verboseflag, 0},
			{"varyz1", no_argument, &varyz1, 1},
			/* These options don’t set a flag. We distinguish them by their indices. */
			{"n0", optional_argument, 0, 'n'},
			{"n1", optional_argument, 0, 'm'},
			{"n2", optional_argument, 0, 'l'},
			{"z1", optional_argument, 0, 'y'},
			{"seed", optional_argument, 0, 's'},
			{"pop", optional_argument, 0, 'p'},
			{"dpop", optional_argument, 0, 'd'},
			{"file", optional_argument, 0, 'f'},
			{"nlive", optional_argument, 0, 'e'},
			{"tobs", optional_argument, 0, 't'},
			{"bins", optional_argument, 0, 'b'},
			{"zpts", optional_argument, 0, 'z'},
			{"outfile", optional_argument, 0, 'o'},
			{"method", optional_argument, 0, 'a'},
			{0, 0, 0, 0}
		};

		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long(argc, argv, ":n:m:l:s:p:d:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1) break;

		switch (c)
		{
			case 'n':
				args->n0 = atof(optarg);
				break;
			case 'm':
				args->n1 = atof(optarg);
				break;
			case 'l':
				args->n2 = atof(optarg);
				break;
			case 'y':
				args->z1 = atof(optarg);
				break;
			case 's':
				args->seed = atol(optarg);
				break;
			case 'p':
				args->popsize = atoi(optarg);
				break;
			case 'd':
				args->datapopsize = atoi(optarg);
				break;
			case 'f':
				strcpy(args->datafile, optarg);
				break;
			case 'e':
				args->nlive = atoi(optarg);
				break;
			case 't':
				args->tobs = atoi(optarg);
				break;
			case 'b':
				args->nbins = atoi(optarg);
				break;
			case 'z':
				args->zpts = atoi(optarg);
				break;
				break;
			case 'o':
				strcpy(args->outfile,optarg);
				break;
			case 'a':
				method = atoi(optarg);
				break;
			case '?':
				break;
			default:
				break;
		}
	}

	args->resume = resflag;
	args->help = helpflag;

	if (nstarflag == 1) args->nstar = true;
	if (ntotalflag == 1) args->ntotal = true;
	if (flatn0flag == 1) args->flatn0 = true;
	if (testflag == 1) args->testpop = true;
	if (zllflag == 1) args->zeroLogLike = true;
	if (varyz1 == 1) args->vary_z1 = true;

	args->verbose = verboseflag;

	if (method>=0 && method<=3) args->method = method;

	if (args->z1 < 0.0) args->z1 = 0.0;
	else if (args->z1 > 10.0) args->z1 = 10.0;
}

double logPoisson(double k, double lambda)
{
	return k * log(lambda) - lambda - (k * log(k) - k + 0.5 * log(2.0 * M_PI * k));
}

double logPoisson2(double k, double lambda)
{
	return (k + 0.5) * log(lambda / k) + (k - lambda);
}

double logPoisson3(double k, double lambda)
{
	double offset = 0.0;
	if (lambda > 0)
	{
		offset = lambda * log(lambda) - lambda - logFactorial(lambda);
	}

	if (k == 0)
	{
		if (lambda == 0)
		{
			return 0.0;
		}
		else
		{
			return -lambda - offset;
		}
	}
	else
	{
		if (lambda == 0)
		{
			return -logFactorial(k) - offset;
		}
		else
		{
			return k * log(lambda) - lambda - logFactorial(k) - offset;
		}
	}
}

void bindetections(double *zdet, long int ndet, double zmin, double zmax, int nbins, int *count)
{
	int i;
	double dz = (zmax - zmin) / (nbins - 1);

	// initialize counts to zero
	for (i = 0; i < nbins; i++)
	{
		count[i] = 0;
	}

	int binid;
	for (i = 0; i < ndet; i++)
	{
		binid = (int) floor(zdet[i] / dz);
		count[binid]++;
	}
}

double logFactorial(double x)
{
	if (x <= 7)
	{
		return log(factorial(x));
	}
	else
	{
		return x * log(x) - x + 0.5 * log(2.0*M_PI*x);
	}
}

double factorial(double x)
{
	if (x == 0)
	{
		return 1.0;
	}
	else
	{
		return x * factorial(x - 1);
	}
}
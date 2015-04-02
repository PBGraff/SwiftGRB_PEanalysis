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

void read_options(int argc, char *argv[], double *n0, double *n1, double *n2, long int *seed, int *resume,
				  int *popsize, int *datpopsize, int *help, char *filename, int *nlive)
{
	int c;
	static int resflag=0, helpflag=0;
	
	while (1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"resume", no_argument, &resflag, 1},
			{"help", no_argument, &helpflag, 1},
			/* These options don’t set a flag. We distinguish them by their indices. */
			{"n0", optional_argument, 0, 'n'},
			{"n1", optional_argument, 0, 'm'},
			{"n2", optional_argument, 0, 'l'},
			{"seed", optional_argument, 0, 's'},
			{"pop", optional_argument, 0, 'p'},
			{"dpop", optional_argument, 0, 'd'},
			{"file", optional_argument, 0, 'f'},
			{"nlive", optional_argument, 0, 'e'},
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
				*n0 = atof(optarg);
				break;
			case 'm':
				*n1 = atof(optarg);
				break;
			case 'l':
				*n2 = atof(optarg);
				break;
			case 's':
				*seed = atol(optarg);
				break;
			case 'p':
				*popsize = atoi(optarg);
				break;
			case 'd':
				*datpopsize = atoi(optarg);
				break;
			case 'f':
				strcpy(filename, optarg);
				break;
			case 'e':
				*nlive = atoi(optarg);
			case '?':
				break;
			default:
				break;
		}
	}

	*resume = resflag;
	*help = helpflag;
}

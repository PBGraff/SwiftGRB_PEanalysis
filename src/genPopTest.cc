#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>

extern "C" {
	#include "utils.h"
	#include "mock_sample_functions.h"
}

#define POPSIZE		10
#define NINPUTS		15
#define NPOP		10

double population[POPSIZE*NINPUTS], zpop[POPSIZE];

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

int main(int argc, char *argv[])
{
	long int seed = (long int) time(NULL);
	double n0, n1, n2;
	double z1=3.6, x=-0.65, y=-3.0, log_lum_star=52.05;

	void read_background_values();

	double t[NPOP], tt[NPOP];

	int i;

	/*FILE *fp1 = fopen("splines_LuminosityDistance.txt","w");
	FILE *fp2 = fopen("splines_RedshiftRate.txt","w");
	double zz, zzld, zzR;
	for ( zz=0.0; zz<=10.0; zz+=1e-2 )
	{
		zzld = log10(lumi_distance(zz));
		fprintf(fp1, "%lf\t%lf\n", zz, zzld);
		fprintf(stderr, "%lf\t%lf\n", zz, zzld);
	}
	for ( zz=0.0; zz<=10.0; zz+=1e-3 )
	{
		zzR = Rate_dz_part(zz);
		fprintf(fp2, "%lf\t%lf\n", zz, zzR);
		fprintf(stderr, "%lf\t%lf\n", zz, zzR);
	}
	fclose(fp1);
	fclose(fp2);
	exit(-1);*/

	load_splines();

	for ( i=0; i<NPOP; i++ )
	{
		clock_t timer = clock();

		n0 = CubeToLogPrior(ran2d(&seed), 0.25, 2.0);
		n1 = CubeToFlatPrior(ran2d(&seed), 1.6, 4.0);
		n2 = CubeToFlatPrior(ran2d(&seed), -4.0, 0.0);

		GeneratePopulation(population, POPSIZE, n0, n1, n2, z1, x, y, log_lum_star, zpop, &seed);

		timer = clock() - timer;

		t[i] = (double) timer / CLOCKS_PER_SEC;
		tt[i] = t[i]*t[i];
	}

	unload_splines();

	double mean=0.0, stdev=0.0;
	for ( i=0 ; i<NPOP; i++ )
	{
		mean += t[i] / NPOP;
		stdev += tt[i] / NPOP;
	}
	stdev -= mean*mean;

	printf("%lf +/- %lf secs per population\n", mean, sqrt(stdev));

	return 0;
}

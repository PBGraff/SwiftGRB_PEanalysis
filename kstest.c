#include "kstest.h"

void kstwo(double *data1, double *data2, long int n1, long int n2, double *d, double *prob)
{
	long int j1=1, j2=1;
	double d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;

	sort(n1,data1);
	sort(n2,data2);

	en1 = n1;
	en2 = n2;
	*d = 0.0;

	while (j1 <= n1 && j2 <= n2)
	{
		d1 = data1[j1-1];
		d2 = data2[j2-1];

		if ( d1 <= d2 ) fn1 = j1++/en1;
		if ( d2 <= d1 ) fn2 = j2++/en2;

		dt = fabs(fn2-fn1);
		if ( dt > *d ) *d = dt;
	}

	en = sqrt(en1*en2/(en1+en2));
	*prob = probks((en + 0.12 + 0.11/en) * (*d));
}

double probks(double alam)
{
	int j;
	double a2,fac=2.0,sum=0.0,term,termbf=0.0;

	a2 = -2.0*alam*alam;
	for ( j=1; j<=100; j++ )
	{
		term = fac * exp(a2*j*j);
		sum += term;
		if ( fabs(term) <= EPSKS1*termbf || fabs(term) <= EPSKS2*sum ) return sum;
		fac = -fac;
		termbf = fabs(term);
	}

	return 1.0;
}

void sort(long int n, double *arr)
{
	// move the array back one to use the 1-indexed algorithm
	arr--;

	long int i,ir=n,j,k,l=1,*istack;
	int jstack=0;
	double a, temp;

	istack = (long int *)malloc(SORTSTACK*sizeof(long int));
	istack--;

	for ( ;; )
	{
		if ( ir-l <= SORTM )
		{
			for ( j=l+1; j<ir; j++ )
			{
				a = arr[j];
				for ( i=j-1; i>=l; i-- )
				{
					if ( arr[i] <= a ) break;
					arr[i+1] = arr[i];
				}
				arr[i+1] = a;
			}
			if ( jstack==0 ) break;
			ir = istack[jstack--];
			l = istack[jstack--];
		}
		else
		{
			k = (l+ir) >> 1;
			SORTSWAP(arr[k],arr[l+1]);
			if ( arr[l] > arr[ir] )
			{
				SORTSWAP(arr[l],arr[ir]);
			}
			if ( arr[l+1] > arr[ir] )
			{
				SORTSWAP(arr[l+1],arr[ir]);
			}
			if ( arr[l] > arr[l+1] )
			{
				SORTSWAP(arr[l],arr[l+1]);
			}
			i = l+1;
			j = ir;
			a = arr[l+1];
			for ( ;; )
			{
				do i++; while ( arr[i] < a );
				do j--; while ( arr[j] > a );
				if ( j < i ) break;
				SORTSWAP(arr[i],arr[j]);
			}
			arr[l+1] = arr[j];
			arr[j] = a;
			jstack += 2;
			if ( jstack > SORTSTACK )
			{
				fprintf(stderr,"SORTSTACK too small.");
				abort();
			}
			if ( ir-i+1 >= j-1 )
			{
				istack[jstack] = ir;
				istack[jstack-1] = i;
				ir = j-1;
			}
			else
			{
				istack[jstack] = j-1;
				istack[jstack-1] = l;
				l=i;
			}
		}
	}

	istack++;
	free(istack);

	// shift the array back into its original position
	arr++;
}

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

void SkyNetKStest(char filename[])
{
	long int i, numlines;

	//count lines in the file
	numlines = countlines(filename);

	// allocate memory
	double *z, *ptrue, *ppred, *z1, *z2;
	long int n1, n2;
	z = (double *)malloc(numlines * sizeof(double));
	z1 = (double *)malloc(numlines * sizeof(double));
	z2 = (double *)malloc(numlines * sizeof(double));
	ptrue = (double *)malloc(numlines * sizeof(double));
	ppred = (double *)malloc(numlines * sizeof(double));

	// read data files of predictions
	char line[500];
	double temp[16];
	FILE *ifp = fopen(filename,"r");
	for ( i=0; i<numlines; i++ )
	{
		fgets(line,500,ifp);
		sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			&temp[0], &z[i], &temp[1], &temp[2], &temp[3], &temp[4], &temp[5], &temp[6], &temp[7],
			&temp[8], &temp[9], &temp[10], &temp[11], &temp[12], &temp[13], &temp[14], &ptrue[i],
			&temp[15], &ppred[i]);
	}
	fclose(ifp);

	// find list of detected z's
	detected(z, ptrue, numlines, 0.5, z1, &n1);
	detected(z, ppred, numlines, 0.5, z2, &n2);

	// calculate K-S test probability
	double d, p;
	kstwo(z1, z2, n1, n2, &d, &p);

	// print results
	printf("N1=%d, N2=%d, KSprob=%lf, logL=%lf\n", n1, n2, p, log(p));

	// free memory
	free(z);
	free(z1);
	free(z2);
	free(ptrue);
	free(ppred);
}

void SciKitLearnKStest(char filename[])
{
	long int i, numlines;

	//count lines in the file
	numlines = countlines(filename);

	// allocate memory
	double *z, *ptrue, *ppred, *z1, *z2;
	long int n1, n2;
	z = (double *)malloc(numlines * sizeof(double));
	z1 = (double *)malloc(numlines * sizeof(double));
	z2 = (double *)malloc(numlines * sizeof(double));
	ptrue = (double *)malloc(numlines * sizeof(double));
	ppred = (double *)malloc(numlines * sizeof(double));

	// read data files of predictions
	char line[100];
	FILE *ifp = fopen(filename,"r");
	for ( i=0; i<numlines; i++ )
	{
		fgets(line,100,ifp);
		sscanf(line,"%lf %lf %lf\n", &z[i], &ptrue[i], &ppred[i]);
	}
	fclose(ifp);

	double d, p;

	// find list of detected z's
	detected(z, ptrue, numlines, 0.5, z1, &n1);
	detected(z, ppred, numlines, 0.5, z2, &n2);

	// calculate K-S test probability
	kstwo(z1, z2, n1, n2, &d, &p);
	
	// print results
	printf("N1=%d, N2=%d, KSprob=%lf, logL=%lf\n", n1, n2, p, log(p));

	for ( i=0; i<40004; i+=10001)
	{
		// find list of detected z's
		detected(&z[i], &ptrue[i], 10001, 0.5, z1, &n1);
		detected(&z[i], &ppred[i], 10001, 0.5, z2, &n2);

		// calculate K-S test probability
		kstwo(z1, z2, n1, n2, &d, &p);

		// print results
		printf("N1=%d, N2=%d, KSprob=%lf, logL=%lf\n", n1, n2, p, log(p));
	}
	// free memory
	free(z);
	free(z1);
	free(z2);
	free(ptrue);
	free(ppred);
}

#ifndef __KSTEST_H__
#define __KSTEST_H__ 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define EPSKS1		0.001
#define EPSKS2		1.0e-8
#define SORTM		7
#define SORTSTACK	500

#define SORTSWAP(a,b)	temp=(a);(a)=(b);(b)=temp;

// perform two-sample K-S test
void kstwo(double *data1, double *data2, long int n1, long int n2, double *d, double *prob);

// Find the probability of a K-S test D,Ne result
double probks(double alam);

// sort an array in place
void sort(long int n, double *arr);

// return an array of detected z from input array of z and detection prob, given a prob threshold
void detected(double *trigz, double *ptrig, long int ntrig, double pth, double *detz, long int *ndet);

// count lines in a file
long int countlines(char filename[]);

// wrapper function to perform end-to-end test given a set of predictions from SkyNet models
void SkyNetKStest(char filename[]);

// wrapper function to perform end-to-end test given a set of predictions from SciKitLearn models
void SciKitLearnKStest(char filename[]);

#endif

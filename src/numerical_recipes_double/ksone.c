#include <math.h>
#define NRANSI
#include "nrutil.h"

void ksone(double data[], unsigned long n, double (*func)(double), double *d,
	double *prob)
{
	double probks(double alam);
	void sort(unsigned long n, double arr[]);
	unsigned long j;
	double dt,en,ff,fn,fo=0.0;

	sort(n,data);
	en=n;
	*d=0.0;
	for (j=1;j<=n;j++) {
		fn=j/en;
		ff=(*func)(data[j]);
		dt=FMAX(fabs(fo-ff),fabs(fn-ff));
		if (dt > *d) *d=dt;
		fo=fn;
	}
	en=sqrt(en);
	*prob=probks((en+0.12+0.11/en)*(*d));
}
#undef NRANSI

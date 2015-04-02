#include <math.h>
#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 10

double qromb1(double (*func)(double), double a, double b)
{
	void polint1(double xa[], double ya[], int n, double x, double *y, double *dy);
	double trapzd1(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd1(func,a,b,j);
		if (j >= K) {
			polint1(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb1");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

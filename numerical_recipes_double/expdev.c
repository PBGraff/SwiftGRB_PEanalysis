#include <math.h>

double expdev(long *idum)
{
	double ran1(long *idum);
	double dum;

	do
		dum=ran1(idum);
	while (dum == 0.0);
	return -log(dum);
}

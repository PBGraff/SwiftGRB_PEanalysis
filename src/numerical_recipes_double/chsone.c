void chsone(double bins[], double ebins[], int nbins, int knstrn, double *df,
	double *chsq, double *prob)
{
	double gammq(double a, double x);
	void nrerror(char error_text[]);
	int j;
	double temp;

	*df=nbins-knstrn;
	*chsq=0.0;
	for (j=1;j<=nbins;j++) {
		if (ebins[j] <= 0.0) nrerror("Bad expected number in chsone");
		temp=bins[j]-ebins[j];
		*chsq += temp*temp/ebins[j];
	}
	*prob=gammq(0.5*(*df),0.5*(*chsq));
}

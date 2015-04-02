void chstwo(double bins1[], double bins2[], int nbins, int knstrn, double *df,
	double *chsq, double *prob)
{
	double gammq(double a, double x);
	int j;
	double temp;

	*df=nbins-knstrn;
	*chsq=0.0;
	for (j=1;j<=nbins;j++)
		if (bins1[j] == 0.0 && bins2[j] == 0.0)
			--(*df);
		else {
			temp=bins1[j]-bins2[j];
			*chsq += temp*temp/(bins1[j]+bins2[j]);
			//printf("%lf %lf\n",bins1[j],bins2[j]);
		}
	*prob=gammq(0.5*(*df),0.5*(*chsq));
}

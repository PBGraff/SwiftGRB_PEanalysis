//#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>

#define NSPLINELD	1001
#define NSPLINEEZ	1001
#define NSPLINERR	10001
#define NSPLINEDF	10001
gsl_interp_accel *accLD, *accRed, *accEz, *accDF;
gsl_spline *splineLD, *splineRed, *splineEz, *splineDF;

#define MPC_IN_M		3.08567758e22
#define GPC_IN_M		(1e3 * MPC_IN_M)
#define C_LIGHT			2.99792458e8
#define HUBBLE0			(1e5 * 0.71 / MPC_IN_M)
#define DH				(C_LIGHT / HUBBLE0 / GPC_IN_M)
#define DH3				(DH * DH * DH)

void load_splines()
{
	long int i;
	
	// Luminosity distribution

	accLD = gsl_interp_accel_alloc();
	splineLD = gsl_spline_alloc(gsl_interp_cspline, NSPLINELD);
	
	double z1[NSPLINELD], ld[NSPLINELD];
	FILE *fp1 = fopen("support_data/splines_LuminosityDistance.txt","r");
	for ( i=0; i<NSPLINELD; i++ )
	{
		fscanf(fp1, "%lf\t%lf\n", &z1[i], &ld[i]);
	}
	fclose(fp1);
	
	gsl_spline_init(splineLD, z1, ld, NSPLINELD);

	// Redshift distribution

	accRed = gsl_interp_accel_alloc();
	splineRed = gsl_spline_alloc(gsl_interp_cspline, NSPLINERR);

	double z2[NSPLINERR], red[NSPLINERR];
	FILE *fp2 = fopen("support_data/splines_RedshiftRate.txt","r");
	for ( i=0; i<NSPLINERR; i++ )
	{
		fscanf(fp2, "%lf\t%lf\n", &z2[i], &red[i]);
	}
	fclose(fp2);

	gsl_spline_init(splineRed, z2, red, NSPLINERR);

	// E(z) cosmology

	accEz = gsl_interp_accel_alloc();
	splineEz = gsl_spline_alloc(gsl_interp_cspline, NSPLINEEZ);

	double z3[NSPLINEEZ], Ez[NSPLINEEZ], temp1, temp2;
	FILE *fp3 = fopen("support_data/splines_Ez.txt","r");
	for ( i=0; i<NSPLINEEZ; i++ )
	{
		fscanf(fp3, "%lf %lf %lf %lf\n", &z3[i], &temp1, &temp2, &Ez[i]);
	}
	fclose(fp3);

	gsl_spline_init(splineEz, z3, Ez, NSPLINEEZ);

	// detection fraction (w.r.t z)

	accDF = gsl_interp_accel_alloc();
	splineDF = gsl_spline_alloc(gsl_interp_cspline, NSPLINEDF);

	double z4[NSPLINEDF], df[NSPLINEDF], temp3;
	FILE *fp4;
	if (runargs.method == NEURALNET)
	{
		fp4 = fopen("support_data/splines_detection_fraction_z.txt","r");
		for ( i=0; i<NSPLINEDF; i++ )
		{
			fscanf(fp4, "%lf %lf %lf\n", &z4[i], &df[i], &temp3);
		}
		fclose(fp4);
	}
	else if (runargs.method == RANDOMFOREST)
	{
		fp4 = fopen("support_data/splines_detection_fraction_z_RF.txt","r");
		for ( i=0; i<NSPLINEDF; i++ )
		{
			fscanf(fp4, "%lf %lf\n", &z4[i], &df[i]);
		}
		fclose(fp4);
	}
	else if (runargs.method == ADABOOST)
	{
		fp4 = fopen("support_data/splines_detection_fraction_z_AB.txt","r");
		for ( i=0; i<NSPLINEDF; i++ )
		{
			fscanf(fp4, "%lf %lf\n", &z4[i], &df[i]);
		}
		fclose(fp4);
	}
	else if (runargs.method == FLUXTHRESH)
	{
		fp4 = fopen("support_data/splines_detection_fraction_z_flux.txt","r");
		for ( i=0; i<NSPLINEDF; i++ )
		{
			fscanf(fp4, "%lf %lf\n", &z4[i], &df[i]);
		}
		fclose(fp4);
	}

	gsl_spline_init(splineDF, z4, df, NSPLINEDF);
}

void unload_splines()
{
	gsl_spline_free(splineLD);
	gsl_interp_accel_free(accLD);
	
	gsl_spline_free(splineRed);
	gsl_interp_accel_free(accRed);

	gsl_spline_free(splineEz);
	gsl_interp_accel_free(accEz);

	gsl_spline_free(splineDF);
	gsl_interp_accel_free(accDF);
}

static inline double MAX(double v1, double v2)
{
	return (((v1) > (v2)) ? (v1) : (v2));
}

double log_sum_exp(double x, double y)
{
	//convert to natural log
	x /= M_LOG10E;
	y /= M_LOG10E;
	
	double mx = MAX(x,y);
	double val = log(exp(x - mx) + exp(y - mx)) + mx;
	
	return val/M_LN10;
}

double log_subtract_exp(double x, double y)
{
	// x >= y must be true!
	if (y>x)
	{
		fprintf(stderr, "%lf < %lf !\n", x, y);
		exit(-1);
	}
	//convert to natural log
	x /= M_LOG10E;
	y /= M_LOG10E;
	double val = x + log1p(-exp(y-x));
	return val/M_LN10;
}

double lumdist_logpowint(double log_lum_min, double log_lum_max, double pwr)
{
	double val;
	if ( pwr > -1.0 )
	{
		val = -pwr * log_lum_star_global - log10(pwr+1.0) + log_subtract_exp((pwr+1.0)*log_lum_max, (pwr+1.0)*log_lum_min);
	}
	else
	{
		val = -pwr * log_lum_star_global - log10(-(pwr+1.0)) + log_subtract_exp((pwr+1.0)*log_lum_min, (pwr+1.0)*log_lum_max);
	}
	//fprintf(stderr, "%lf -> %lf : %lf\n", log_lum_min, log_lum_max, val);
	return val;
}

double lumdist_logpowinv(double prob, double log_lum_min, double log_lum_max, double pwr)
{
	double temp = lumdist_logpowint(log_lum_min, log_lum_max, pwr) + log10(prob);
	temp += pwr * log_lum_star_global;
	
	if ( pwr > -1.0 )
	{
		temp += log10(pwr+1.0);
		temp = log_sum_exp((pwr+1.0)*log_lum_min, temp);
	}
	else
	{
		temp += log10(-(pwr+1.0));
		temp = log_subtract_exp((pwr+1.0)*log_lum_min, temp);
	}

	return temp/(pwr+1.0);
}

double luminosity_distr_fast(double prob, double log_lum_min, double log_lum_max)
{
	//fprintf(stderr, "%lf, %lf, %lf, %lf, %lf, %lf\n", prob, log_lum_min, log_lum_max, log_lum_star_global, alpha1_global, beta1_global);

	double log_lum_int_lower = lumdist_logpowint(log_lum_min, log_lum_star_global, alpha1_global);
	double log_lum_int_upper = lumdist_logpowint(log_lum_star_global, log_lum_max, beta1_global);
	double log_lum_int_total = log_sum_exp(log_lum_int_lower, log_lum_int_upper);
	double lower_frac = pow(10.0, log_lum_int_lower - log_lum_int_total), upper_frac = 1.0 - lower_frac;

	//fprintf(stderr, "%lf, %lf, %lf, %lf, %lf\n", log_lum_int_lower, log_lum_int_upper, log_lum_int_total, lower_frac, upper_frac);

	if ( prob <= lower_frac )
	{
		return lumdist_logpowinv(prob/lower_frac, log_lum_min, log_lum_star_global, alpha1_global);
	}
	else
	{
		return lumdist_logpowinv((prob - lower_frac)/upper_frac, log_lum_star_global, log_lum_max, beta1_global);
	}
}

double lum2flux_integrand(double E, void *params)
{
	double *vals = (double *)params;
	double alpha = vals[0];
	double beta = vals[1];
	double E0 = vals[2];
	double Ebreak = vals[3];

	if ( E < Ebreak )
	{
		return pow(E/100.0,alpha) * exp(-E/E0) * E;
	}
	else
	{
		return pow((alpha-beta)*E0/100.0,alpha-beta) * exp(beta-alpha) * pow(E/100.0,beta) * E;
	}
}

double lum2flux_integral_numeric(double alpha, double beta, double Epeak, double Emin, double Emax)
{
	double E0 = Epeak / (2.0+alpha);
	double Ebreak = (alpha-beta)*E0;

	//gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	double result, error;
	unsigned long int neval;
	double vals[4] = {alpha, beta, E0, Ebreak};
	gsl_function F;
  	F.function = &lum2flux_integrand;
  	F.params = &vals[0];

  	gsl_set_error_handler_off();
  	//gsl_integration_qag(&F, Emin, Emax, 1e-7, 1e-6, 1000, 1, w, &result, &error);
  	gsl_integration_qng(&F, Emin, Emax, 1e-7, 1e-6, &result, &error, &neval);
  	//fprintf(stderr, "%d evaluations\t", neval);

  	//gsl_integration_workspace_free (w);

  	return result;
}

double Redshift_distribution_unnormalized(double z, double n1, double n2, double z1)
{
	if (z <= z1)
	{
		return pow(1.0 + z, n1);
	} else {
		return pow(1.0 + z1, n1 - n2) * pow(1.0 + z, n2);
	}
}

double Redshift_distribution_normalized(double z, double n0, double n1, double n2, double z1)
{
	return n0 * Redshift_distribution_unnormalized(z, n1, n2, z1);
}

double Redshift_rescaled(double z, void *params)
{
	double *pars = (double *) params;
	double n0 = pars[0];
	double n1 = pars[1];
	double n2 = pars[2];
	double z1 = pars[3];

	double Rprime = Redshift_distribution_normalized(z, n0, n1, n2, z1) / (1.0 + z);

	Rprime *= DH3 * gsl_spline_eval(splineEz, z, accEz);

	return Rprime;
}

double Redshift_rejection_sampler(long int *seed, double n0, double n1, double n2, double z1)
{
	double pars[4] = {n0, n1, n2, z1};
	
	double Rzmax = Redshift_rescaled(z1, (void *) &pars[0]);
	
	double z = 0.0, p, x;
	for (;;)
	{
		z = ran2d(seed) * 10.0;
		
		p = Redshift_rescaled(z, (void *) &pars[0]) / Rzmax;

		x = ran2d(seed);

		if (x < p) break;
	}

	return z;
}

double GRBNumberIntegral(double n0, double n1, double n2, double z1)
{
	double pars[4] = {n0, n1, n2, z1};

	double result, error;
	unsigned long int neval;
	gsl_function F;
  	F.function = &Redshift_rescaled;
  	F.params = (void *) &pars[0];

  	gsl_set_error_handler_off();
  	gsl_integration_qng(&F, ZMIN, ZMAX, 1e-7, 1e-6, &result, &error, &neval);

  	result *= 4.0 * M_PI;

  	return result;
}

double GRBRate(double z, double n0, double n1, double n2, double z1)
{
	double Rprime = Redshift_distribution_normalized(z, n0, n1, n2, z1) / (1.0 + z);

	Rprime *= DH3 * gsl_spline_eval(splineEz, z, accEz);

	Rprime *= 4.0 * M_PI * runargs.tobs * gsl_spline_eval(splineDF, z, accDF) / 6.0;

	return Rprime;
}

double GRBRateFunc(double z, void *params)
{
	double *pars = (double *) params;
	double n0 = pars[0];
	double n1 = pars[1];
	double n2 = pars[2];
	double z1 = pars[3];

	return GRBRate(z, n0, n1, n2, z1);
}

double GRBRateIntegral(double n0, double n1, double n2, double z1)
{
	double pars[4] = {n0, n1, n2, z1};

	double result, error;
	unsigned long int neval;
	gsl_function F;
  	F.function = &GRBRateFunc;
  	F.params = (void *) &pars[0];

  	gsl_set_error_handler_off();
  	gsl_integration_qng(&F, ZMIN, ZMAX, 1e-7, 1e-6, &result, &error, &neval);

  	return result;
}

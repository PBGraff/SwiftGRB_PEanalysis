//#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>

#define NSPLINELD	1001
#define NSPLINEEZ	1001
#define NSPLINERR	10001
gsl_interp_accel *accLD, *accRed, *accEz, *accEzI;
gsl_spline *splineLD, *splineRed, *splineEz, *splineEzI;

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
	accEzI = gsl_interp_accel_alloc();
	splineEzI = gsl_spline_alloc(gsl_interp_cspline, NSPLINEEZ);

	double z3[NSPLINEEZ], Ez[NSPLINEEZ], EzI[NSPLINEEZ];
	FILE *fp3 = fopen("support_data/splines_Ez.txt","r");
	for ( i=0; i<NSPLINEEZ; i++ )
	{
		fscanf(fp3, "%lf %lf %lf\n", &z3[i], &Ez[i], &EzI[i]);
	}
	fclose(fp3);

	gsl_spline_init(splineEz, z3, Ez, NSPLINEEZ);
	gsl_spline_init(splineEzI, z3, EzI, NSPLINEEZ);
}

void unload_splines()
{
	gsl_spline_free(splineLD);
	gsl_interp_accel_free(accLD);
	
	gsl_spline_free(splineRed);
	gsl_interp_accel_free(accRed);

	gsl_spline_free(splineEz);
	gsl_interp_accel_free(accEz);
	gsl_spline_free(splineEzI);
	gsl_interp_accel_free(accEzI);
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

/*double Redshift_denominator = -1.0;

double redshift_integrand(double z, void *params)
{
	return Redshift(z) * gsl_spline_eval(splineRed, z, accRed);
}

double redshift_cumprob_func(double z, void *params)
{
	double result, error;
	unsigned long int neval;
	gsl_function F;
  	F.function = &redshift_integrand;

  	gsl_integration_qng(&F, Z_i, z, 1e-7, 1e-6, &result, &error, &neval);

  	return result/Redshift_denominator - Prob_Z;
}

double redshift_distribution_fast()
{
	double z, z_lo, z_hi;

	if (Redshift_denominator == -1.0)
	{
		double error;
		unsigned long int neval;
		gsl_function F;
  		F.function = &redshift_integrand;

  		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

		//gsl_integration_qng(&F, Z_i, Z_f, 1e-5, 1e-4, &Redshift_denominator, &error, &neval);
		gsl_integration_qags(&F, Z_i, Z_f, 1e-5, 1e-4, 1000, w, &Redshift_denominator, &error);

		gsl_integration_workspace_free (w);
	}

	gsl_function F;
	F.function = redshift_cumprob_func;
	const gsl_root_fsolver_type * T = gsl_root_fsolver_bisection;
	gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, Z_i, Z_f);

	int status, iter=0, maxiter=100;
	do {
		iter++;
		status = gsl_root_fsolver_iterate (s);
		z = gsl_root_fsolver_root(s);
		z_lo = gsl_root_fsolver_x_lower(s);
		z_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(z_lo, z_hi, 0, 0.001);
	} while ( status == GSL_CONTINUE && iter < maxiter );

	gsl_root_fsolver_free(s);
}*/

/*
double lum2flux_integral_lower(double alpha, double E0, double Emin, double Emax)
{
	//double val = gsl_sf_gamma_inc_Q(alpha+2.0, Emin/E0) + gsl_sf_gamma_inc_P(alpha+2.0, Emax/E0) - 1.0;
	//val *= 100.0 * pow(E0/100.0, alpha+1.0) * gsl_sf_gamma(alpha+2.0);
	double val = gsl_sf_gamma_inc_Q(alpha+2.0, Emin/E0) - gsl_sf_gamma_inc_Q(alpha+2.0, Emax/E0);
	val = log(val) + log(100.0) + (alpha+1.0)*log(E0/100.0) + gsl_sf_lngamma(alpha+2.0);
	return val;
}

double lum2flux_integral_upper(double alpha, double beta, double E0, double Emin, double Emax)
{
	double val = (pow(Emax, beta+2.0) - pow(Emin, beta+2.0)) / (beta+2.0);
	//val *= pow((alpha-beta)*E0/100.0, alpha-beta) * exp(beta-alpha) * pow(100.0, -beta);
	val = log(val) + (alpha-beta)*log((alpha-beta)*E0/100.0) + (beta-alpha) - beta*log(100.0);
	return val;
}

double lum2flux_integral(double alpha, double beta, double Epeak, double Emin, double Emax)
{
	double val, val2;
	double E0 = Epeak / (2.0+alpha);
	double Ebreak = (alpha-beta)*E0;

	if (Emin < Ebreak)
	{
		if (Emax <= Ebreak)
		{
			val = lum2flux_integral_lower(alpha, E0, Emin, Emax);
		}
		else
		{
			val = lum2flux_integral_lower(alpha, E0, Emin, Ebreak);
			val2 = lum2flux_integral_upper(alpha, beta, E0, Ebreak, Emax);
			val = log_sum_exp(val, val2);
		}
	}
	else
	{
		val = lum2flux_integral_upper(alpha, beta, E0, Emin, Emax);
	}

	return exp(val);
}*/

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

double Redshift_distribution_unnormalized(double z)
{
	if (z <= z1_global)
	{
		return pow(1.0 + z, n1_global);
	} else {
		return pow(1.0 + z1_global, n1_global - n2_global) * pow(1.0 + z, n2_global);
	}
}

double Redshift_rejection_sampler(long int &seed)
{
	double Rzmax = Redshift_distribution_unnormalized(z1_global);
	double z = 0.0, p, x;
	for (;;)
	{
		z = ran2d(seed) * 10.0;
		
		p = Redshift_distribution_unnormalized(z) / Rzmax;

		x = ran2d(seed);

		if (x < p) break;
	}

	return z;
}



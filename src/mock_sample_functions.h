#ifndef __MOCK_H__
#define __MOCK_H__ 1

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
extern "C" {
	#include "nr.h"
	#include "nrutil.h"
}

#define Squ(x) ((x)*(x))
#define Cub(x) ((x)*(x)*(x))
#define C 2.99792458e+10 // cm/s
#define Plankh  6.6260755e-27 // erg*s
#define Boltzk 1.380658e-16 // erg/k
//#define MeV2erg 1.602e-6
#define Mpc2cm 3.09e+24
#define Gpc2cm 3.09e+27
#define H0 (71.0/(3.086e+19))  // 1/s
#define H0_unit 71.0 /* hubble cons in km/s/Mpc */
#define Pi 3.141592654
#define Omega_m 0.274
#define Omega_r 0.0
#define Omega_lambda 0.726

double Log_lum_i, Log_lum_f;
double Prob_lum;
double Z_i, Z_f;
double Prob_Z;
double Log_Epeak_max, Log_Epeak_min;
double Prob_logEpeak;
double Log_Liso_global, Z_global;
double Tau1_pulse1,Tau2_pulse1,Tau1_pulse2,Tau2_pulse2,Tau1_pulse3,Tau2_pulse3;
double *Lambdaa,*Qea,*Qea2,Lambda,Qe,Qeap1,Qeapn;
int N;

//parameters for GRB property distributions
double rate_GRB_0_global, z1_global, n1_global, n2_global;   //redshift distribution
double lum_star_global, log_lum_star_global, alpha1_global, beta1_global;  //luminosity distribution
double lum_step_global; //luminosity evolution
char Epeak_type[1000];  //Epeak distribution
char Lum_evo_type[1000]; //Luminosity evolution

double H(double z);
double rcomi(double z);
double rcom(double z);
double D_lum(double z);
double LF(double logL);
double lum_prob_func(double log_lum);
double lum_distribution(long *seed);
double Redshift(double z);
double Rate_dz(double z);
double Redshift_prob_func(double z);
double redshift_distribution(long *seed);
double alpha_distribution(long *seed);
double beta_distribution(long *seed);
double log_Epeak_func(double log_Epeak);
double logEpeak_prob_func(double log_Epeak);
double Epeak_distribution(long *seed);
int Angle_table(double num_ran);
double sfr(double z);
double lightcurve_pulse1(double t);
double lightcurve_pulse2(double t);
double lightcurve_pulse3(double t);
double lumi_distance (double z1);
double light_travel_time (double z1);
double ne_calc (double alpha, double beta, double e0, double ene, double *ne);
double ne_integral (double alpha, double beta, double e0, double emin, double emax, double *ne_int_value);
double lum2flux (double alpha, double beta, double epeak, double liso, double z, double band_start, double band_stop, double *output_flux);

#endif

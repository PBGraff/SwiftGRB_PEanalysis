#include "mock_sample_functions.h"
#include "myrand.h"

#include "poputils.c"

// functions for general cosmology
/*double H(double z) // *H0 s
{
        double val=sqrt(0.274*Cub(1.0+z)+0.726);
        
	      if(val < 0.0) {printf("H wrong!");}

        return val;
}*/

double rcomi(double z) // *(C/H0) cm
{
        //double val,y;
        //double H(double z);

        //y=sqrt(0.274*Cub(1.0+z)+0.726);
        //val=1.0/H(z);

        //return val;
        //return 1.0/H(z);

        return 1.0/sqrt(0.274*Cub(1.0+z)+0.726);
}

double rcom(double z) // *C/H0  cm
{
        /*double val,rcomi(double z);

        if(z == 0.0){val=0.0;}
        else{
                val=qromb(rcomi,0.0,z);
        }
        if(val < 0.0) {printf("z=%lf, rcom wrong!",z);}
        return val;*/
        
        if (z==0.0) return 0.0;
        
        double val = qromb(rcomi,0.0,z);
        if(val < 0.0) {printf("z=%lf, rcom wrong!",z);}
        return val;
}


double D_lum(double z)
{
        /*double rcom(double z);
        double d_lum_Mpc,val;
        val=(1.0+z)*C/H0*rcom(z);
        return val;*/
        return (1.0+z)*C/H0*rcom(z);
}


// functions for burst luminosity distribution

double LF(double logL)  //dN/dLogL from Wanderman et al. 2010
{
        double lum,lum_star,lum_factor,val;

	//lum_star_global = pow(10.0,52.58); //erg/s orignal one from Wanderman et la. 2010
	//alpha1_global = -0.18; //original one from Wanderman et al. 2010
        //beta1_global = -1.59; //original one from Wandermane et al. 2010
       
	lum = pow(10.0,logL);

	//luminosity evolution option
	if(strcmp(Lum_evo_type,"lum_evo_log")==0){	
		//lum_factor = (int)Z_global;
		//lum_factor = Z_global*lum_step_global; //try1
		//lum_factor = Z_global*Z_global*lum_step_global; //try2
		lum_factor = lum_step_global*log10(Z_global); //try3
		lum_star = lum_star_global*pow(10.0,lum_factor);
	}
	if(strcmp(Lum_evo_type,"lum_evo_linear")==0){
                //lum_factor = (int)Z_global;
                lum_factor = Z_global*lum_step_global; //try1
                //lum_factor = Z_global*Z_global*lum_step_global; //try2
		lum_star = lum_star_global*pow(10.0,lum_factor);
        }
	if(strcmp(Lum_evo_type,"no_lum_evo")==0){
		lum_star = lum_star_global;
	}

	//luminosity function 	 
	if(lum<lum_star){
		val=pow((lum/lum_star),alpha1_global);
	}
	else{
		val=pow((lum/lum_star),beta1_global);
	}
	val = val*lum*log(10.0);
	
	return val;

}
	
double lum_prob_func(double log_lum)
{
	double val;
	double LF_integrate, LF_deno;
	
	LF_integrate = qromb(LF,Log_lum_i,log_lum);
        LF_deno = qromb(LF, Log_lum_i, Log_lum_f);
        val = LF_integrate/LF_deno - Prob_lum;

	return val;
}

double lum_distribution(long *seed)
{
	/*double Log_lum;
	
	Log_lum_i = 30.0; //erg/s
        Log_lum_f = 70.0; //erg/s
        Prob_lum = ran2d(seed);
        zbrac(lum_prob_func,&Log_lum_i,&Log_lum_f);
        Log_lum=rtbis(lum_prob_func,Log_lum_i,Log_lum_f,1.0e-3);

	//Log_lum = (70.0-30.0)*ran2d(seed)+30.0;

	return Log_lum;*/

  Log_lum_i = 30.0;
  Log_lum_f = 70.0;

  Prob_lum = ran2d(seed);
  //printf("%lf\t", Prob_lum);

  double lumval = luminosity_distr_fast(Prob_lum, Log_lum_i, Log_lum_f);

  /*zbrac(lum_prob_func,&Log_lum_i,&Log_lum_f);
  printf("%g/%g\t", Log_lum_i, Log_lum_f);
  double Log_lum=rtbis(lum_prob_func,Log_lum_i,Log_lum_f,1.0e-9);

  printf("%g =?= %g\t(%g)\n", Log_lum, lumval, (Log_lum-lumval)/Log_lum);*/

  return lumval;

}

// functions for burst redshift distribution

double Redshift(double z)  //dN/dz from Wanderman et al. 2010
{
        double val;


        //rate_GRB_0_global = 1.3; // Gpc^{-3} yr^{-1}      
        
	//parameters from Table 1, Wanderman et al. 2010
	// 1/3 bins
	//z1_global = 3.45;
        //n1_global = 1.74;
        //n2_global = -1.47;
	
	//1/2 bins
	//z1_global = 3.11; //original one from Wandermane et al. 2010
	//n1_global = 2.07;
	//n2_global = -1.36;
	
        if(z<=z1_global){
                val=pow((1.0+z),n1_global);
        }
        else{
                val=pow((1+z1_global),(n1_global-n2_global))*pow((1+z),n2_global);
        }

	val=val*rate_GRB_0_global;

	return val;

/*	//use the shape of star-formation rate
	double a,c,b,d;
        double rho_0,alpha,beta,gamma,z1,z2,eta,y1,y2,y3;
        double log_rho,norm;

	//sfr from Horiuchi, Beacom, Dwek 09 [M_sun/yr/Mpc^3]
        rho_0 = 0.0178; alpha = 3.4; beta = -0.3; gamma = -3.5;
        z1 = 1.0; z2 = 4.0;
        eta = -10.0;
        b = pow((1.0+z1),(1.0-alpha/beta));
        c = pow((1.0+z1),((beta-alpha)/gamma))*pow((1.0+z2),(1.0-beta/gamma));
        y1=pow((1.0+z),alpha*eta);
        y2=pow(((1.0+z)/b),beta*eta);
        y3=pow(((1.0+z)/c),gamma*eta);
        val = rho_0*pow((y1+y2+y3),1.0/eta);

        //calculate snr from sfr 
        val=0.007*val;   // number/yr/Mpc^3

	norm = rate_GRB_0_global/0.0001246;  //renormalized to GRB_rate at z=0
        val=norm*val; 
	return val;
*/

/*	//sfr from Hopkins and Beacom 06, piecewise fit
        double log_rho, z1, z2, norm;
	z1 = 1.04;
        z2 = 4.48;
        if(z <= z1){
                log_rho = 3.28*log10(1.0+z)-1.82;
	}
        if(z > z1 && z <= z2){
                log_rho = -0.26*log10(1.0+z)-0.724;
	}
        if(z > z2){
                log_rho = -8.0*log10(1.0+z)+4.99;
	}
        val = pow(10.0,log_rho);

	norm = rate_GRB_0_global/0.0151356124844;  //renormalized to GRB_rate at z=0
	val=norm*val; 
        return val;
*/

}

double Rate_dz(double z) // number/yr/sr/redshift
{
        double rate_dzdomega, rate_dzdomega2;
        rate_dzdomega=Cub(C/H0)*Squ(rcom(z))*rcomi(z)/(1.0+z)*Redshift(z)/Cub(Gpc2cm);
        //rate_dzdomega2 = gsl_spline_eval(splineRed, z, accRed) * Redshift(z);
        //fprintf(stderr, "%g =?= %g\t(%g)\n", rate_dzdomega, rate_dzdomega2, (rate_dzdomega - rate_dzdomega2)/rate_dzdomega);
        //if ( fabs(rate_dzdomega-rate_dzdomega2)>1e-3 ) fprintf(stderr, "WRONG!\t%g\t%g\t%g\n", z, rate_dzdomega, rate_dzdomega2);
        return rate_dzdomega;
}

double Rate_dz_part(double z) // number/yr/sr/redshift
{
        double rate_dzdomega;
        rate_dzdomega=Cub(C/H0)*Squ(rcom(z))*rcomi(z)/(1.0+z)/Cub(Gpc2cm);

        return rate_dzdomega;
}

double Redshift_deno=-1.0;
double Redshift_prob_func(double z)
{
        double val;
        double Redshift_integrate;

        Redshift_integrate = qromb(Rate_dz,Z_i,z);
        //Redshift_deno = qromb(Rate_dz, Z_i, Z_f);
        val = Redshift_integrate/Redshift_deno - Prob_Z;

        return val;
}

double redshift_distribution(long *seed)
{
	double z;
	Z_i = 1.0e-3;
  Z_f = 10.0;
  // this needs to be re-done each time n0, n1, or n2 change
  //if (Redshift_deno==-1.0) 
  Redshift_deno = qromb(Rate_dz, Z_i, Z_f);
  Prob_Z = ran2d(seed);
  zbrac(Redshift_prob_func,&Z_i,&Z_f);
  z=rtbis(Redshift_prob_func,Z_i,Z_f,1.0e-3);

	//flat distribution
	//z = (10.0-1.0e-3)*ran2d(seed)+1.0e-3;

  //double z2 = redshift_distribution_fast();
  //fprintf(stderr, "%g =?= %g\t(%g)\n", z, z2, (z-z2)/z);
	
	return z;
}

//find alpha; Gaussian distribution based on Sakamoto et al. 2009
double alpha_distribution(long *seed)  
{
	double alpha_avg,alpha_sigma,val;
	
	alpha_avg = -0.87;
        alpha_sigma = 0.33;
        val = alpha_avg + alpha_sigma*gasdev(seed);

	
	//flat distribution
	//val = (2.0-(-2.0))*ran2d(seed)+(-2.0);

	return val;
}

//find beta; Gaussian distribution based on Sakamoto et al. 2009
double beta_distribution(long *seed)  
{
	double beta_avg,beta_sigma,val;

	beta_avg = -2.36;
        beta_sigma = 0.31;
        val = beta_avg + beta_sigma*gasdev(seed);

	//flat distribution
	//val = (0.0-(-10.0))*ran2d(seed)+(-10.0);
	
	return val;
}

//find E_peak
double log_Epeak_func(double log_Epeak)
{
        double val;

        if(log_Epeak >= 1.0){
                val = 1.0;
        }
        else{
                val = log_Epeak;
        }

        return val;

}

double logEpeak_prob_func(double log_Epeak)
{
        double val;
        double logEpeak_integrate, logEpeak_deno;

        logEpeak_integrate = qromb(log_Epeak_func,Log_Epeak_max,log_Epeak);
        logEpeak_deno = qromb(log_Epeak_func, Log_Epeak_max, Log_Epeak_min);
        val = logEpeak_integrate/logEpeak_deno - Prob_logEpeak;

        return val;
}

double Epeak_distribution(long *seed)
{
	double E_peak_min,E_peak_max,log_E_peak,E_peak;
	double log_Liso_52, Liso_52;
	int type_check;

	type_check = 0;

	E_peak_min = 1.0; //keV
        Log_Epeak_min = log10(E_peak_min);
        E_peak_max = 10.0e+3; //keV  
	Log_Epeak_max = log10(E_peak_max);
        
	//Flat Epeak distribution in linear space
	if(strcmp(Epeak_type,"flat_linear")==0){
		E_peak = ran2d(seed)*(E_peak_max - E_peak_min) + E_peak_min;
		E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
		type_check = 1;
	}


	//flat Epeak distribution in log space
        if(strcmp(Epeak_type,"flat_log")==0){
		log_E_peak = ran2d(seed)*(Log_Epeak_max - Log_Epeak_min) + Log_Epeak_min;
		E_peak = pow(10.0,log_E_peak);
		E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
		type_check = 1;
	}

	//Gaussian Epeak distribution
	if(strcmp(Epeak_type,"Gaussian")==0){
		log_E_peak = log10(300.0) + 1.0*gasdev(seed);
		E_peak = pow(10.0,log_E_peak);
		E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
		type_check = 1;
	}

	//Epeak_Liso relation from Yunetoku et al. 2004
        if(strcmp(Epeak_type,"Yonetoku")==0){
		log_Liso_52 = Log_Liso_global - 52.0;
	        Liso_52 = pow(10.0,log_Liso_52);
	        E_peak = 1.0/sqrt(2.34e-5)*sqrt(Liso_52);
		E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
		type_check = 1;
	}

	//Epeak_Liso relation from a modified Yunetoku relation
        if(strcmp(Epeak_type,"Yonetoku_mod20")==0){
                log_Liso_52 = Log_Liso_global - 52.0;
                Liso_52 = pow(10.0,log_Liso_52);
                E_peak = 1.0/sqrt(2.34e-5)*sqrt(Liso_52);
		E_peak = E_peak*2.0;
                E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
                type_check = 1;
        }

	//Epeak_Liso relation from a modified Yunetoku relation
        if(strcmp(Epeak_type,"Yonetoku_mod18")==0){
                log_Liso_52 = Log_Liso_global - 52.0;
                Liso_52 = pow(10.0,log_Liso_52);
                E_peak = 1.0/sqrt(2.34e-5)*sqrt(Liso_52);
                E_peak = E_peak*1.8;
                E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
                type_check = 1;
        }

	//Epeak_Liso relation from a modified Yunetoku relation
        if(strcmp(Epeak_type,"Yonetoku_mod17")==0){
                log_Liso_52 = Log_Liso_global - 52.0;
                Liso_52 = pow(10.0,log_Liso_52);
                E_peak = 1.0/sqrt(2.34e-5)*sqrt(Liso_52);
                E_peak = E_peak*1.7;
                E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
                type_check = 1;
        }

	//Epeak_Liso relation from a modified Yunetoku relation
        if(strcmp(Epeak_type,"Yonetoku_mod15")==0){
                log_Liso_52 = Log_Liso_global - 52.0;
                Liso_52 = pow(10.0,log_Liso_52);
                E_peak = 1.0/sqrt(2.34e-5)*sqrt(Liso_52);
                E_peak = E_peak*1.5;
                E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
                type_check = 1;
        }

	//Epeak_Liso relation from a modified Yunetoku relation
        if(strcmp(Epeak_type,"Yonetoku_mod11")==0){
                log_Liso_52 = Log_Liso_global - 52.0;
                Liso_52 = pow(10.0,log_Liso_52);
                E_peak = 1.0/sqrt(2.34e-5)*sqrt(Liso_52);
                E_peak = E_peak*1.1;
                E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
                type_check = 1;
        }
	//use other function as specified in log_Epeak_func
	if(strcmp(Epeak_type,"specified")==0){
		Prob_logEpeak = ran2d(seed);
       		zbrac(logEpeak_prob_func,&Log_Epeak_min, &Log_Epeak_max);
       		log_E_peak=rtbis(logEpeak_prob_func,Log_Epeak_min,Log_Epeak_max,1.0e-3);
		E_peak = pow(10.0,log_E_peak);
		E_peak = E_peak/(1.0+Z_global); //move Epeak to observor's frame, because Epeak in sim_lc_v2 is in observor's frame.
		type_check = 1;
	}


	if(type_check == 0){
		printf("%s is an invalid Epeak_type\n", Epeak_type);
	}	

	return E_peak;
}

int Angle_table(double num_ran)
{
	//double prob[31]={0.005042,0.025210,0.030252,0.023529,0.010084,0.003361,0.030252,0.048739,0.107563,0.060504,0.015126,0.001681,0.000000,0.026891,0.084034,0.110924,0.095798,0.011765,0.008403,0.003361,0.018487,0.048739,0.075630,0.050420,0.020168,0.001681,0.010084,0.013445,0.035294,0.018487,0.005042};
	//double total_prob;
	int table[31]={1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,30,31,32,33};
	int i, i_prob, val;

	i_prob = (int)num_ran;
	val = table[i_prob];
	return val;	
	
}

double snr(double z) 
{
        double val,a,c,b,d;
        double rho_0,alpha,beta,gamma,z1,z2,eta,y1,y2,y3;
        double log_rho;

        //sfr from Horiuchi, Beacom, Dwek 09 [M_sun/yr/Mpc^3]
        rho_0 = 0.0178; alpha = 3.4; beta = -0.3; gamma = -3.5;
        z1 = 1.0; z2 = 4.0;
        eta = -10.0;
        b = pow((1.0+z1),(1.0-alpha/beta));
        c = pow((1.0+z1),((beta-alpha)/gamma))*pow((1.0+z2),(1.0-beta/gamma));
        y1=pow((1.0+z),alpha*eta);
        y2=pow(((1.0+z)/b),beta*eta);
        y3=pow(((1.0+z)/c),gamma*eta);
        val = rho_0*pow((y1+y2+y3),1.0/eta);

	//calculate snr from sfr 
	val=0.007*val;   // number/yr/Mpc^3
	return val;
}

double lightcurve_pulse1(double t) // Norris et al. (2005) lightcurve, without normalization factor (added later)
{
	double val;
	double mu,lambda,lum,y;
	
	mu = sqrt(Tau1_pulse1/Tau2_pulse1);
	lambda = exp(2.0*mu);
	y = Tau1_pulse1/t + t/Tau2_pulse1;
	val = lambda*exp(-1.0*y);

	return val;

}

double lightcurve_pulse2(double t) // Norris et al. (2005) lightcurve, without normalization factor (added later)
{
        double val;
        double mu,lambda,lum,y;

        mu = sqrt(Tau1_pulse2/Tau2_pulse2);
        lambda = exp(2.0*mu);
        y = Tau1_pulse2/t + t/Tau2_pulse2;
        val = lambda*exp(-1.0*y);

        return val;

}

double lightcurve_pulse3(double t) // Norris et al. (2005) lightcurve, without normalization factor (added later)
{
        double val;
        double mu,lambda,lum,y;

        mu = sqrt(Tau1_pulse3/Tau2_pulse3);
        lambda = exp(2.0*mu);
        y = Tau1_pulse3/t + t/Tau2_pulse3;
        val = lambda*exp(-1.0*y);

        return val;

}

/// Below are functions needed for routine lum2flux (convert luminosity to flux), written by Taka

double lumi_distance (double z1)
{
  double z, dz, dm, z_integral;
  double hubble_cons_cgs;
  double dl;

  hubble_cons_cgs = H0_unit * 1e5 / (1e6 * 3.0857e18);

  z_integral = 0.0;
  dz=1.0e-7;

  for(z=0.0; z <= z1; z+=dz){
    z_integral += 1.0 / sqrt((1+z)*(1+z)*(1+Omega_m*z)-z*(2+z)*Omega_lambda) * dz;
  }

  dm = 2.9979e10 / hubble_cons_cgs * z_integral;

  dl = (1+z1) * dm;

  return dl;
}

double light_travel_time (double z1)
{
  double ltt, z, z_integral, dz;
  double hubble_cons_cgs;

  hubble_cons_cgs = H0_unit * 1e5 / (1e6 * 3.0857e18);

  z_integral = 0.0;
  dz=0.00001;

  for(z=0.0; z <= z1; z+=dz){
    z_integral += 1.0 / sqrt((1.0+z)*(1.0+z)*(1.0+Omega_m*z)-z*(2.0+z)*
                             Omega_lambda) / (1.0+z) * dz;
    //    printf("%lf %lf\n", z, z_integral);
  }

  //  printf("hubble = %e\n", 1.0/hubble_cons_cgs);

  ltt = z_integral / hubble_cons_cgs;

  return ltt;
}

double ne_calc (double alpha, double beta, double e0, double ene, double *ne)
{
  double norm_ene;
  double ne_calc;
  norm_ene = 100.0;
  ne_calc = 0.0;

  if (ene <= (alpha - beta) * e0){
    ne_calc = pow(ene/norm_ene, alpha) * exp(-ene/e0);
  }
  else {
    ne_calc = pow((alpha - beta) * e0 / norm_ene, alpha - beta) *
      exp(beta - alpha) * pow(ene / norm_ene, beta);
  }

  *ne = ne_calc;
}

double ne_integral (double alpha, double beta, double e0,
                    double emin, double emax, double *ne_int_value)
{
  double ne_integral_value;
  double ene, dene;
  double ne;

  ne_integral_value = 0.0;

  dene=0.001;

  for(ene = emin+dene/2.0; ene <= emax+dene/2.0; ene+=dene) {
    ne_calc(alpha, beta, e0, ene, &ne);
    /* printf("%lf %lf\n", ene, ne); */
    /* printf("ne = %e\n", ne); */
    ne_integral_value += ene * ne * dene;
  }

  *ne_int_value = ne_integral_value;
}

double lum2flux (double alpha, double beta, double epeak, double liso, double z, double band_start, double band_stop, double *output_flux)
{
  double A, e0, fluence;
  double emin, emax;
  double ne_int_value;
  double fluence_kev, flux, flux_obs;
  double duration;
  double ld, epeak_src;
  double Log_lum_peak, t_bin,time;
  int i;

  //double ld2 = lumi_distance(z);
  ld = pow(10.0, gsl_spline_eval(splineLD, z, accLD));
  //fprintf(stderr, "%g =?= %g\t(%g)\n", log10(ld2), log10(ld), (log10(ld2)-log10(ld))/log10(ld2));

  flux = liso / (4.0 * Pi * ld * ld);

  //  emin    = atof(argv[5]);
  //  emax    = atof(argv[6]);

  // epeak = epeak_src / (1.0+z);

  e0 = epeak / (2.0+alpha);

  emin = 1.0 / (1.+z);
  emax = 10000.0/ (1.+z);

  /*ne_integral(alpha, beta, e0, emin, emax, &ne_int_value);

  A = flux / (ne_int_value * 1.6022e-9);

  //  flux = A * ne_int_value * 1.6022e-9;
  // fluence = flux * duration;
  //  fluence = flux;
  //   printf("%lf %lf %e\n", emin, emax, fluence);

  ne_integral(alpha, beta, e0, band_start, band_stop, &ne_int_value);

  flux_obs = A * ne_int_value * 1.6022e-9;*/

  //*output_flux = flux_obs;

  /*double test = lum2flux_integral(alpha, beta, epeak, emin, emax);
  A = flux / (test * 1.6022e-9);
  test = lum2flux_integral(alpha, beta, epeak, band_start, band_stop);
  test *= A * 1.6022e-9;
  //printf("%e =?= %e\t(%e)\n", flux_obs, test, (flux_obs-test)/flux_obs);
  *output_flux = test;*/

  double test = lum2flux_integral_numeric(alpha, beta, epeak, emin, emax);
  A = flux / (test * 1.6022e-9);
  test = lum2flux_integral_numeric(alpha, beta, epeak, band_start, band_stop);
  test *= A * 1.6022e-9;
  //printf("%e =?= %e\t(%e)\n", flux_obs, test, (flux_obs-test)/flux_obs);
  *output_flux = test;

  return 0.0;

}


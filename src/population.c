double grid_thetas[31] = {50.694, 40.675, 35.017749, 40.675039, 50.662786, 56.990242, 46.628525, 31.393231,
						  19.267549, 31.393231, 46.628525, 56.990242, 56.291844, 44.994283, 26.557574,
						  0.075701, 26.557574, 44.994283, 56.291844, 56.990242, 46.628525, 31.393231,
						  19.267549, 31.393231, 46.628525, 56.990242, 50.662786, 40.675039, 35.017749,
						  40.675039, 50.662786};

void angle_to_grid(int angle, double *r, double *phi)
{
	double x = (angle%7-3.0)/2.0;
    double y = ((angle/7)-2.0)/3.0;
    *r = sqrt(x*x+y*y);
    *phi = atan2(y,x);
}

// obtain the background values
#define NBACKGROUND		371
double bkgd_logvals[4*NBACKGROUND];
void read_background_values()
{
	char bkgdfilename[] = "support_data/background_count_table.txt";
	char line[100], bgd_name[50];
	int i = 0;
	
	FILE *bkgdfile = fopen(bkgdfilename,"r");
	
	fgets(line, 100, bkgdfile);
	
	for ( i=0; i<NBACKGROUND; i++ )
	{
		fgets(line, 100, bkgdfile);
		sscanf(line,"%s %lf %lf %lf %lf\n", bgd_name, &bkgd_logvals[4*i], &bkgd_logvals[4*i+1],
			   &bkgd_logvals[4*i+2], &bkgd_logvals[4*i+3]);
	}
	
	for ( i=0; i<4*NBACKGROUND; i++ )
	{
		bkgd_logvals[i] = log10(bkgd_logvals[i]);
	}
	
	fclose(bkgdfile);
}

void GeneratePopulation(double population[], long int population_size, double n0, double n1, double n2, double z1,
						double x, double y, double log_lum_star, double zpop[], long int *seed)
{
	/*
		Input parameters:
		0  log_L 					log-luminosity
		1  z 						redshift
		2  r 						radius of grid id from center
		3  phi 						aziumthal angle of grid id
		4  bin_size_emit 			source time bin size
		5  alpha 					Band func param
		6  beta 					Band func param
		7  E_peak 					energy peak (log this)
		8  bgd_15-25keV 			bkg in band (log this)
		9  bgd_15-50keV 			bkg in band (log this)
		10 bgd25-100keV 			bkg in band (log this)
		11 bgd50-350keV 			bkg in band (log this)
		12 theta 					incoming angle
		13 flux 					flux of burst (log this)
		14 ndet 					number of detectors active
	*/
	
	/* UNUSED variables
	double Eiso, lum_peak, lum_t, log_E_min, log_E_max, log_E_peak, alpha_avg, alpha_sigma, beta_avg, beta_sigma;
	double E_min, E_max, bgd_ran, lc_ran, lambdadata, qedata, time_start, time_end, N_t_bin_cal, Qe_max;
	double E_peak_t, test, *time, *Qe_array;
	int N_t_bin, bgd_ran_i,size, lc_ran_i, N_total;
	*/
	//long int seed = (long int) time(NULL);
	double z, Log_lum, t_bin_obs, grid_r, grid_phi, lum, num_ran, alpha, beta, E_peak, flux_15_150, E0, t_bin, ndet_ran;
	int startid, ipop, i, ndet, angle;

	// set global variables
	rate_GRB_0_global = n0;
	n1_global = n1;
	n2_global = n2;
	z1_global = z1;
	alpha1_global = x;
	beta1_global = y;
	lum_star_global = pow(10.0, log_lum_star);
	log_lum_star_global = log_lum_star;
	lum_step_global = 0.0;

	Redshift_deno=-1.0;

	sprintf(Epeak_type, "Yonetoku_mod18");
	sprintf(Lum_evo_type, "no_lum_evo");

	//time bin in obserbed frame
	t_bin_obs = 1600.0*1.0e-3; //s

	//fprintf(stderr, "Done ");
	for ( ipop=0; ipop<population_size; ipop++ )
	{
		//find z
		do {
			if (runargs.flatzpop) {
				z = ran2d(seed) * (runargs.zbin_max - runargs.zbin_min) + runargs.zbin_min;
			} else {
				//z = redshift_distribution(seed);
				z = Redshift_rejection_sampler(seed, n0, n1, n2, z1);
			}
		} while ( std::isnan(z)==1 );

		//passing z to all the subroutines in mock_sample_functions.c
		Z_global = z;

		//find grid-id
		num_ran = 31.0*ran2d(seed);
		angle = Angle_table(num_ran);
		
		i=11;
		do {
			//if can't match requirements after 10 tries, find new luminosity
            if(i>10)
            {
                //find luminosity
                //passsing luminosity to all subroutines in mock_sample_functions.c
                Log_lum = lum_distribution(seed);
                Log_Liso_global = Log_lum; //passing luminosity to subroutines
                lum = pow(10.0, Log_lum);
                i=1;
            }
			
			//find alpha in the observed frame; Gaussian distribution based on Sakamoto et al. 2009
			//make sure alpha is greater than -2, otherwise E0 = E_peak/(2+alpha) is negative
            do {
                alpha = alpha_distribution(seed);
            } while (alpha <= -2.0);

			//find beta in the observed frame; Gaussian distribution based on Sakamoto et al. 2009
			//get new beta
            beta = beta_distribution(seed);

			//find E_peak in the observed frame
			//get new Epeak
            E_peak = Epeak_distribution(seed);
			E0 = E_peak/(2.0+alpha);
                
			lum2flux(alpha,beta,E_peak,lum,z,15.0,150.0,&flux_15_150);
			
			i++;

			//below are some boundary conditions to make sure: 
			//(1) E0 is in range(0.01,10000), as set in bat_sim_lc_evo.sh
			//(2) E_peak >= 0.1, as set in calc_liso_norm100kev_wmap_evo.c, which is used in bat_sim_lc_evo.sh
			//(3) burst does not have nan flux or flux greater than 1.0e+10
		} while ( E_peak < 0.1 || E0 >= 1.0e+4 || E0 <=1.0e-2 || flux_15_150 > 1.0e+10 || flux_15_150 == 0.0 || 
				  std::isnan(flux_15_150)==1 );

		int bkgd_ran = (int) (NBACKGROUND * ran2d(seed));

		t_bin = t_bin_obs/(1.0+z); //This is the observational bin size redshifted to the emitted frame

		//find ndet
        ndet_ran = 10000.0*ran2d(seed)+20000.0;
		ndet = (int)ndet_ran;

		angle_to_grid(angle, &grid_r, &grid_phi);

		startid = ipop * NINPUTS;

		population[startid + 0]  = Log_lum;
		population[startid + 1]  = z;
		population[startid + 2]  = grid_r;
		population[startid + 3]  = grid_phi;
		population[startid + 4]  = t_bin;
		population[startid + 5]  = alpha;
		population[startid + 6]  = beta;
		population[startid + 7]  = log10(E_peak);
		population[startid + 8]  = bkgd_logvals[4*bkgd_ran];
		population[startid + 9]  = bkgd_logvals[4*bkgd_ran+1];
		population[startid + 10] = bkgd_logvals[4*bkgd_ran+2];
		population[startid + 11] = bkgd_logvals[4*bkgd_ran+3];
		population[startid + 12] = grid_thetas[(int) num_ran];
		population[startid + 13] = log10(flux_15_150);
		population[startid + 14] = 26884.0;

		zpop[ipop] = z;

		/*if ( (ipop + 1) % 100 == 0 )
		{
			fprintf(stderr, "%d / %d done\n", ipop+1, population_size);
		}*/
	}
	//fprintf(stderr, "\n");
}

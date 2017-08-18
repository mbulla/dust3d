 /*
 
 Effect of dust on reddening
 This is a 3D version of Goobar 2008 and Amanullah & Goobar 2011 

 Mattia Bulla, 2016-11-15
 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "mylib.h"
#include "draine.h"

int main(int argc, char *argv[]) {

	long seed;
	seed = time(NULL);
	rng = gsl_rng_alloc(gsl_rng_ran3);
	gsl_rng_set (rng, seed);

	const float B_eff = 0.4405, V_eff = 0.547; // Effective wavelength (in microns)
	const double flux2mag = 1.085736; 	
	const int nmax=50000;             // Maximum number of interactions per photon
	
	char *dusttype = NULL, *is_dusttype = NULL;
	char str;
	long int nph;
	int iwave,iphot,B_i,V_i,is_B_i,is_V_i,i,k,istep,inside_inner,inside_outer,is;
	unsigned int dustresol=1,nwavel,nextphoton=0;	
	double axis_ratio,rinn,rout,minlambda,maxlambda;
	double al_B,al_V,a_B,a_V,is_al_B,is_al_V,is_a_B,is_a_V,factor,is_factor,scatconst,tau,s_tau,s_forw,s_back,gran;
	double g_dust,costh,th,phi,path,tottheta,thetafirst,pos_angles[2],dir_angles[2];
	double x[3],axis_inn[3],axis_out[3],ndir[3],ndir_out[3],ndir_back[3],xinn_DLOS[3];
	double ebv,ebv_cs,ebv_is,ebv_tot,ebv_unsc,ebv_cs_err,ebv_tot_err,ebv_unsc_err;
	double path_inside;
	double tmin=0,tmax=1e3;//e3;
	int dobin=0,thetabin,totbin,bin_index,break_inn;
	double step_costh,dOmega_source,dOmega_collect,norm_factor;//z,D
	unsigned short *nscatter=NULL,save_pkt_info=0,savescatter=0,saveangles=0,savefirstangle=0;
	scatext *draine, *is_draine;
	char filename[128];
	FILE *ccstest, *cistest, *scstest, *packet_info;

	int dopolar =0;
	polarisation *polar;
  	float pvalue;
	double I,Q,U,avgI,avgQ,avgU;

	clock_t begin, end, ttmp;
	double time_spent;
	begin = clock();

	// -------------  Flags ------------- 
	save_pkt_info = 1;	
	if (save_pkt_info) {
		saveangles = 1;
		savescatter = 1;
		savefirstangle = 1;
	}
	dopolar = 0;
	dobin = 0;
	break_inn = 0;  // If = 1, kill pkts that enter the inner region (for is, back-scattered packets are killed)
	is = 0; // = 1 if you want to switch IS dust on

	// ------------- Input parameters ------------- 
	nph = 5e7;
	ebv = 1.5;	
	if (is == 1) ebv_is = 0.3;
	else ebv_is = 0.0 ;

  	// Dust type1
  	dusttype = "MW3";
	if (is == 1) is_dusttype = "MW3";

	// See notes 18/08/17
	dOmega_source = 4 * PI ; 
	dOmega_collect = 4 * PI ;	
	
    // Geometry
    rout = 1000.0;
    rinn = 0.95;

    axis_ratio = 1;

	if (rinn < 0.0 || rinn > 1.0) {
		fprintf(stderr,"Error: -R must be given a value between 0.0 and 1.0!\n");
		return 1;
	} 
	else rinn *= rout;

	axis_inn[0] = rinn ;
	axis_inn[1] = rinn * axis_ratio ;
	axis_inn[2] = rinn * axis_ratio ;

	axis_out[0] = rout ;
	axis_out[1] = rout * axis_ratio ;
	axis_out[2] = rout * axis_ratio ;

	// Observer viewing angle

	if (dobin) {

		// For the moment this is assuming axial symmetry
		thetabin = 7;
		totbin = thetabin ;
		step_costh = 2. / thetabin;
	}

	else {

		thetabin = 1;
		totbin = 1 ;
		step_costh = 2;	
	}
	
    // Wavelength range (0.17-2.5 for SN2011fe template)
    minlambda = 0.25;//0.17;
    maxlambda = 0.85;//2.5; 
    if (minlambda < 0.0001 || minlambda > 9.9999) {
    	fprintf(stderr,"Warning: wavelengths < 0.0001 or > 9.9999 mu not supported!\n");
    	minlambda = 0.0001;
	}
  	if (maxlambda > 9.9999 || maxlambda < 0.0001) {
    	fprintf(stderr,"Warning: wavelengths < 0.0001 or > 9.9999 mu not supported!\n");
    	maxlambda = 9.9999;
  	}
  	if (minlambda>B_eff) {
    	fprintf(stderr,"Warning: B_eff < wavelengths < V_eff for factor to make sense!\n");
    	minlambda = B_eff;
  	}
  	if (maxlambda<V_eff) {
    	fprintf(stderr,"Warning: B_eff < wavelengths < V_eff for factor to make sense!\n");
    	maxlambda = V_eff;
  	}

	// Normalization factor (see notes 18/08/17)
	norm_factor =  dOmega_collect / dOmega_source ;

	// Initialize the scattering, extinction, albedo table
	draine = NULL;
	draine = init_draine_extscat(dusttype,minlambda,maxlambda,dustresol,&nwavel);
	if (draine == NULL) {
		fprintf(stderr,"Error: failed to load Draine's extinction/scattering table!\n");
		exit(1);
	}

	if (is == 1) {

		is_draine = NULL;
		is_draine = init_draine_extscat(is_dusttype,minlambda,maxlambda,dustresol,&nwavel);
	
		if (is_draine == NULL) {
			fprintf(stderr,"Error: failed to load Draine's extinction/scattering table!\n");
			exit(1);
		}

		find_draine_lambda_idx(is_draine,nwavel,B_eff,&is_B_i);
		find_draine_lambda_idx(is_draine,nwavel,V_eff,&is_V_i);
	}

	// Initialize polarisation table
	//
	polar = NULL;
	if (dopolar) {
		polar = init_draine_pol(dusttype);
		if (polar == NULL) {
			fprintf(stderr,"Error: failed to load the polarisation table!\n");
			exit(1);
		}
	}

	// Packet info file
	snprintf(filename,sizeof filename,"work/06-sample/runs/packets_%s_rinn%g_ratio%g_ebv%g_dOmsource%gpi_dOmcoll%gpi_Nph%.0e_BVonly.out"
		,dusttype,rinn/rout,axis_ratio,ebv,dOmega_source/PI,dOmega_collect/PI,(double) nph);   

	//snprintf(filename,sizeof filename,"test.out");   

    packet_info=fopen(filename,"a+");
    fseek(packet_info, 0, SEEK_END);
    if (ftell(packet_info) != 0) {

    	printf("\n%s already exists!!\n\n",filename);
    	printf("Should I overwrite this [Y/N]? ");
	    scanf("%s",&str);
		if (strcmp(&str,"Y") != 0) return -1;
		packet_info=fopen(filename,"w");
	}

	printf("\nOutput file: %s \n\n",filename);
	
	fprintf(packet_info,"%s \t\t # Dust type\n",dusttype);
	fprintf(packet_info,"%ld \t\t # Number of packets\n",nph);
	fprintf(packet_info,"%d \t\t # Number of wavelengths\n",nwavel);
	fprintf(packet_info,"%g %g \t # Min and max wavelengths\n",minlambda,maxlambda);
	fprintf(packet_info,"%g \t\t # Input E(B-V)\n",ebv);
	fprintf(packet_info,"%g \t\t # Binning factor \n",norm_factor);
	//fprintf(packet_info,"%g \t\t # D \n",D);
	//fprintf(packet_info,"%g \t\t # z \n",z);
	fprintf(packet_info,"%g \t\t # dOmega_collect \n",dOmega_collect);


	// Find which indices correspond to B and V filters
	find_draine_lambda_idx(draine,nwavel,B_eff,&B_i);
	find_draine_lambda_idx(draine,nwavel,V_eff,&V_i);

	// Things to keep track of for each wavelength
	float *fractime = NULL, *stokesI = NULL, *stokesQ = NULL, *stokesU = NULL,
		*totangle = NULL, *firstangle = NULL;
	double *scatdust, *absdust, *avt, *tau_is, *a_x_unsc, *a_x_cs, *a_x_is, *a_x_tot, *a_x_away, *a_x_abs, *a_x_ins;
	double *a_x_unsc_avg, *a_x_cs_avg, *a_x_is_avg, *a_x_tot_avg, *a_x_away_avg, *a_x_abs_avg, *a_x_ins_avg;
	unsigned int *sum_cs, *sum_cs_avg, *sum_unsc, *sum_unsc_avg, *sum_abs, *sum_away, *sum_ins;
	unsigned int *sum_abs_avg, *sum_away_avg, *sum_ins_avg;

	sum_cs_avg   = (unsigned int*) malloc(sizeof(unsigned int) * nwavel );
	sum_cs       = (unsigned int*) malloc(sizeof(unsigned int) * nwavel * totbin );
	sum_unsc_avg  = (unsigned int*) malloc(sizeof(unsigned int) * nwavel );
	sum_unsc      = (unsigned int*) malloc(sizeof(unsigned int) * nwavel * totbin );
	sum_away_avg   = (unsigned int*) malloc(sizeof(unsigned int) * nwavel );
	sum_away       = (unsigned int*) malloc(sizeof(unsigned int) * nwavel * totbin );
	sum_abs_avg  = (unsigned int*) malloc(sizeof(unsigned int) * nwavel );
	sum_abs      = (unsigned int*) malloc(sizeof(unsigned int) * nwavel * totbin );
	sum_ins_avg  = (unsigned int*) malloc(sizeof(unsigned int) * nwavel );
	sum_ins      = (unsigned int*) malloc(sizeof(unsigned int) * nwavel * totbin );

	avt          = (double*) malloc(sizeof(double) * nwavel);
	scatdust     = (double*) malloc(sizeof(double) * nwavel);
	absdust      = (double*) malloc(sizeof(double) * nwavel);
	tau_is       = (double*) malloc(sizeof(double) * nwavel);
	
	a_x_cs     = (double*) malloc(sizeof(double) * nwavel * totbin);
	a_x_cs_avg = (double*) malloc(sizeof(double) * nwavel);
	a_x_is     = (double*) malloc(sizeof(double) * nwavel * totbin);
	a_x_is_avg = (double*) malloc(sizeof(double) * nwavel);
	a_x_tot     = (double*) malloc(sizeof(double) * nwavel * totbin);
	a_x_tot_avg = (double*) malloc(sizeof(double) * nwavel);
	a_x_unsc     = (double*) malloc(sizeof(double) * nwavel * totbin);
	a_x_unsc_avg = (double*) malloc(sizeof(double) * nwavel);
	a_x_away     = (double*) malloc(sizeof(double) * nwavel * totbin);
	a_x_away_avg = (double*) malloc(sizeof(double) * nwavel);
	a_x_abs     = (double*) malloc(sizeof(double) * nwavel * totbin);
	a_x_abs_avg = (double*) malloc(sizeof(double) * nwavel);
	a_x_ins     = (double*) malloc(sizeof(double) * nwavel * totbin);
	a_x_ins_avg = (double*) malloc(sizeof(double) * nwavel);

	fractime     = (float*) malloc(sizeof(float) * nph);
	if (saveangles) totangle = (float*) malloc(sizeof(float) * nph);
	if (savescatter) nscatter = (unsigned short*) malloc(sizeof(unsigned short) * nph);
	if (savefirstangle) firstangle = (float*) malloc(sizeof(float)* nph);
	if (dopolar) {
		stokesI = (float*) malloc(sizeof(float)* nph);
		stokesQ = (float*) malloc(sizeof(float)* nph);
		stokesU = (float*) malloc(sizeof(float)* nph);
	}	
	

	// Circumstellar dust

	al_B = draine[B_i].albedo;
	al_V = draine[V_i].albedo;
	a_B  = draine[B_i].abs_k;  // This is alpha_B in the notes
	a_V  = draine[V_i].abs_k;  // This is alpha_V in the notes

	factor = (double) a_B/(1.- al_B) - a_V/(1. - al_V);

	// This is the number density (actually density times M in the notes)
	if (axis_inn[0]>axis_inn[2])
		scatconst = ebv/flux2mag/(axis_out[2]-axis_inn[2])/factor; 
	else
		scatconst = ebv/flux2mag/(axis_out[0]-axis_inn[0])/factor; 

	if (al_B == 1 && al_V == 1) 
		printf("Warning: no absorption! Equations are meaningless for albedo = 1 \n"); 

	// Interstellar dust

	if (is ==1) {

		is_al_B = is_draine[is_B_i].albedo;
		is_al_V = is_draine[is_V_i].albedo;
		is_a_B  = is_draine[is_B_i].abs_k;  // This is alpha_B in the notes
		is_a_V  = is_draine[is_V_i].abs_k;  // This is alpha_V in the notes

		is_factor = (double) is_a_B/(1.- is_al_B) - is_a_V/(1. - is_al_V);
	}


	// ------------------------------------------------
	// --------------- START PROPAGATION --------------
	// ------------------------------------------------

	for (iwave=0;iwave<nwavel;iwave++) {   // loop over wavelengths
	
		// Initialize
		
		for (i=0;i<totbin;i++) {
			sum_cs[i * nwavel + iwave]    = 0;
			sum_abs[i * nwavel + iwave]    = 0;
			sum_away[i * nwavel + iwave]    = 0;
			sum_unsc[i * nwavel + iwave]   = 0;
			sum_ins[i * nwavel + iwave]   = 0;
		}

		// CS opacity

		//printf("%g %g \n",draine[iwave].lambda,ebv / flux2mag / factor *  draine[iwave].abs_k / ( 1 - draine[iwave].albedo ));

		// IS opacity
		if (is == 1) 
			tau_is[iwave] = ebv_is / flux2mag / is_factor *  is_draine[iwave].abs_k / ( 1 - is_draine[iwave].albedo ) ;		
		else 
			tau_is[iwave] = 0.0 ;

		sum_cs_avg[iwave]      = 0.0 ;
		sum_away_avg[iwave]      = 0 ;
		sum_unsc_avg[iwave]    = 0.0 ;
		sum_abs_avg[iwave]     = 0 ;
		sum_ins_avg[iwave]     = 0 ;
		avt[iwave]      	   = 0.;
	
		avgI = 0;
		avgQ = 0;
		avgU = 0;


		g_dust = (double) draine[iwave].cost;
		
	    // Mean free path (lambda_s and lambda_a)
		scatdust[iwave] = ((double) 1./draine[iwave].abs_k*(1./draine[iwave].albedo-1.))/scatconst;
		absdust[iwave]  = ((double) 1./draine[iwave].abs_k)/scatconst;
		
		//g_dust = 0.6 ;
		//draine[iwave].albedo = 0.6 ;

		for (iphot=0;iphot<nph;iphot++) {   // loop over photons

			// Generate an unpolarised packet
			I = 1;
			Q = 0;
			U = 0;

			path = 0.;
			nextphoton = 0;
			fractime[iphot] = 0.0; 
			if (saveangles) tottheta = 0.0;
			if (savefirstangle) thetafirst = 0.0;
			pos_angles[0] = 0.0;
			pos_angles[1] = 0.0;
			dir_angles[0] = 0.0;
			dir_angles[1] = 0.0;

			// Generate packet
			source(&dOmega_source,axis_inn,x,ndir,xinn_DLOS);

			//theta0 = acos( norm_plus_dot(x,xinn_DLOS) );

			// Increment path (we have moved packet from 0 to inner surface)
			path += sqrt( pow(x[0],2.) + pow(x[1],2.) + pow(x[2],2.) ) ;

			// Propagate photon until it is absorbed or leaves the dust sphere
			istep = 0;
			
			while (istep < nmax && nextphoton == 0) {	

				// Generate random optical depth and calculate path length 
				tau = - log( 1 - gsl_rng_uniform(rng) );

				s_tau = tau_to_s(x,xinn_DLOS,ndir,&rinn,&rout,&tau,&scatdust[iwave],&absdust[iwave],&dOmega_source) ;
				//s_tau = tau / ( 1. / absdust[iwave] + 1. / scatdust[iwave] );

				// Move packet to the new position
				x[0] = x[0] + s_tau * ndir[0] ;
				x[1] = x[1] + s_tau * ndir[1] ;
				x[2] = x[2] + s_tau * ndir[2] ;

				inside_inner = check(x,axis_inn);
				inside_outer = check(x,axis_out);

				// Packet leaves the sphere
				if (inside_inner == 0 && inside_outer == 0) {

				    nextphoton = 1 ;

					// Move packet back to the previous position
					x[0] = x[0] - s_tau * ndir[0] ;
					x[1] = x[1] - s_tau * ndir[1] ;
					x[2] = x[2] - s_tau * ndir[2] ;

					// You already went back. Now increment path (s+d in the notes)
					path += finalpath(&istep,x,ndir,axis_out);
					// Packet left sphere without interacting at all
				    if (istep==0) {

						escape(xinn_DLOS,x,pos_angles) ;
						escape(xinn_DLOS,ndir,dir_angles) ;

						if ( pos_angles[0] <= dOmega_source / 4. && pos_angles[1] <= dOmega_source / 4. 
							&& dir_angles[0] <= dOmega_collect / 4. && dir_angles[1] <= dOmega_collect / 4. ) { 

							if ( path/rout>=tmin && path/rout<=tmax ) {
					    		sum_unsc_avg[iwave] += 1 ;
					    		sum_cs_avg[iwave] += 1 ;
								
								// Save Stokes parameters
								if (dopolar) {

									stokesI[iphot] = I;
									stokesQ[iphot] = Q;
									stokesU[iphot] = U;

									avgI += I;
									avgQ += Q;
									avgU += U;
								}

							}

							bin_index = (1 - ndir[2]) / step_costh ;
							if (bin_index <= thetabin) {
								sum_unsc[bin_index * nwavel + iwave] += 1 ;
								sum_cs[bin_index * nwavel + iwave] += 1 ;
							}
							else  
								printf("Packet escaped but was not collected in an angular bin %d %d \n",bin_index,thetabin);
							
							// When packet leaves the sphere, thetafirst is updated only if no interaction occurred
							if (savefirstangle == 1) thetafirst = 0.0 ; 

							
						}

				    	// To ensure that packets that are not selected have fractime = -2 (see below)
				    	else {

				    		sum_away_avg[iwave] += 1 ;
				    		path = - rout * 2;
				    	}

				    }

				    else {

				    	// -------- Selecting packets ------
				
						escape(xinn_DLOS,x,pos_angles) ;
						escape(xinn_DLOS,ndir,dir_angles) ;

						if ( pos_angles[0] <= dOmega_source / 4. && pos_angles[1] <= dOmega_source / 4. 
							&& dir_angles[0] <= dOmega_collect / 4. && dir_angles[1] <= dOmega_collect / 4. ) { 

				    		if ( path/rout>=tmin && path/rout<=tmax ) {
 								
 								sum_cs_avg[iwave] += 1 ;

	 							// Save Stokes parameters
								if (dopolar) {

									stokesI[iphot] = I;
									stokesQ[iphot] = Q;
									stokesU[iphot] = U;

									avgI += I;
									avgQ += Q;
									avgU += U;
								}
							}

							bin_index = (1 - ndir[2]) / step_costh ;
							if (bin_index <= thetabin) sum_cs[bin_index * nwavel + iwave] += 1 ;
							else printf("Packet escaped but was not collected in an angular bin \n");
							
				    	}

				    	// To ensure that packets that are not selected have fractime = -2 (see below)
				    	else {

				    		sum_away_avg[iwave] += 1 ;
				    		path = - rout * 2;
				    	}

				    }	    
				}

				//Packet entered the inner region --> move to inner surface and then to interaction
				else if (inside_inner == 1 && inside_outer == 1) {

					if (break_inn) {

						nextphoton = 1 ;

						sum_ins_avg[iwave] += 1 ;
						sum_ins[iwave] += 1 ;

						if (saveangles) tottheta = 0.0;
						if (savefirstangle == 1) thetafirst = 0.0 ; 

						I = 1;
						Q = 0;
						U = 0;

						// To ensure that packets that passed thorugh the inner region have fractime = -3 (see below)
						path = -rout * 3;

						break ;
					}

					// Path inside (see notes at 02/12/16)
					ndir_back[0] = -ndir[0];
					ndir_back[1] = -ndir[1];
					ndir_back[2] = -ndir[2];

					// s_back + s_forw = 2q in Rahman's note for spherical case
					s_forw = dist_to_ellipsoid(axis_inn,x,ndir);
					s_back = dist_to_ellipsoid(axis_inn,x,ndir_back);

					// Increment path to inner surface
					path += s_tau + s_forw;

					//if (s_forw + s_back > 2*rinn) 
					//	fprintf(stderr,"Path travelled inside the inner region is too long %g %g %g\n",s_forw,s_back,rinn);

					// Move to inner surface (s_forw) and then to interaction (s_back)
					x[0] = x[0] + (s_forw + s_back) * ndir[0] ;
					x[1] = x[1] + (s_forw + s_back) * ndir[1] ;
					x[2] = x[2] + (s_forw + s_back) * ndir[2] ;

					inside_outer = check(x,axis_out);

					// Packet leaves the sphere
					if (inside_outer == 0) {

						nextphoton = 1 ;
						
						// Move packet back to the previous position (to inner surface)
						x[0] = x[0] - s_back * ndir[0] ;
						x[1] = x[1] - s_back * ndir[1] ;
						x[2] = x[2] - s_back * ndir[2] ;

						if (fabs(sqrt(pow(x[0]/axis_inn[0],2.) + pow(x[1]/axis_inn[1],2.) + pow(x[2]/axis_inn[2],2.)) - 1.)>1e-4) 
							fprintf(stderr,"Packet should be on inner surface and is not \n");


				    	// -------- Selecting packets ------

						escape(xinn_DLOS,x,pos_angles) ;
						escape(xinn_DLOS,ndir,dir_angles) ;

						if ( pos_angles[0] <= dOmega_source / 4. && pos_angles[1] <= dOmega_source / 4. 
							&& dir_angles[0] <= dOmega_collect / 4. && dir_angles[1] <= dOmega_collect / 4. ) { 

				    		if ( path/rout>=tmin && path/rout<=tmax ) {
 								
 								sum_cs_avg[iwave] += 1 ;

	 							// Save Stokes parameters
								if (dopolar) {

									stokesI[iphot] = I;
									stokesQ[iphot] = Q;
									stokesU[iphot] = U;

									avgI += I;
									avgQ += Q;
									avgU += U;
								}
							}

							bin_index = (1 - ndir[2]) / step_costh ;
							if (bin_index <= thetabin) sum_cs[bin_index * nwavel + iwave] += 1 ;
							else printf("Packet escaped but was not collected in an angular bin \n");		

							// You already went back. Now increment path (s+d in the notes)
							path += finalpath(&istep,x,ndir,axis_out);
				    	}

				    	// To ensure that packets that are not selected have fractime = -2 (see below)
				    	else {

				    		sum_away_avg[iwave] += 1 ;
				    		path = - rout * 2;
				    	}


						// This should never happened (at istep=0, packets is on inner surface and will not cross 
						// the inner region unless it scatters at least once: istep > 0)
						if (istep==0) fprintf(stderr,"This should never happen! \n") ;
					  
					}
					
					// Packet interacts with dust
					else {

						gran = gsl_rng_uniform(rng);
						
						// Absorption
						if (gran > draine[iwave].albedo) {

							nextphoton = 1 ;

							sum_abs_avg[iwave] += 1 ;
							sum_abs[iwave] += 1 ;

							if (saveangles) tottheta = 0.0;
							if (savefirstangle == 1) thetafirst = 0.0 ; 

							I = 1;
							Q = 0;
							U = 0;

							// To ensure that absorbed packets have fractime = -1 (see below)
							path = -rout;
							
						}
						
						// Dust scattering
						else {

							istep++;							
							nextphoton = 0 ;

					   		/* Costheta. As in Code+1995, we do not sample this angle
							   from the polarised scattering matrix and using a rejection
							   technique (not so efficient). Instead we sample the angle
							   from the unpolarised scattering matrix */ 
							costh = cthenyey_greenstein(g_dust); 
					   		th = acos(costh);
							phi = 2 * PI * gsl_rng_uniform(rng);  // phi

					   		// Increment theta
					   		if (saveangles) tottheta += fabs(th);
					   		// First angle if istep==1
					   		if (istep==1 && savefirstangle) thetafirst = th;

					   		// New direction
					   		newdir(&th,&phi,ndir,ndir_out);

					   		// Increment path from rinn to interaction point
					   		path += s_back ;

					   		// Update Stokes parameters
							if (dopolar) {
	     						
	     						pvalue = draine_pvalue(polar,draine[iwave].lambda,th*180./PI);

								if (pvalue == -1) 
									fprintf(stderr,"Warning: failed to calculate P for %2f microns and %1f deg!\n",draine[iwave].lambda,th*180./PI);
	      
					   			update_stokes(ndir,ndir_out,&costh,&g_dust,&pvalue,&I,&Q,&U);				   			
					   		}

					   		// Update new direction
					   		ndir[0] = ndir_out[0];
					   		ndir[1] = ndir_out[1];
					   		ndir[2] = ndir_out[2];

						}	

					}				
				}

				// Packet interacts in the dust region
				else {

					// If this is not the first step...
					// First check that you did not arrive here passing
					// through the inner region (see notes 02/12/16)
					// through_inner: walk backwards and check if you went inside inner region
					
					path_inside = through_inner(axis_inn,&s_tau,x,ndir) ;
					
					if (istep!=0) {
						
						if (path_inside > 0) {

							// Need to move packet of path_inside
							x[0] = x[0] + path_inside * ndir[0] ;
							x[1] = x[1] + path_inside * ndir[1] ;
							x[2] = x[2] + path_inside * ndir[2] ;

							if (break_inn) {
								
								nextphoton = 1 ;
						
								sum_ins_avg[iwave] += 1 ;
								sum_ins[iwave] += 1 ;

								if (saveangles) tottheta = 0.0;
								if (savefirstangle == 1) thetafirst = 0.0 ;

								I = 1;
								Q = 0;
								U = 0;

								// To ensure that packets that passed thorugh the inner region have fractime = -3 (see below)
								path = -rout * 3;
								
								break ;
							}
						
						}
					}

					inside_outer = check(x,axis_out);			

					// Packet leaves the sphere
					if (inside_outer == 0) {

					    path += s_tau - path_inside;
						
						nextphoton = 1 ;

					    // Move packet back to the previous position
						x[0] = x[0] - path_inside * ndir[0] ;
						x[1] = x[1] - path_inside * ndir[1] ;
						x[2] = x[2] - path_inside * ndir[2] ;


				    	// -------- Selecting packets ------

						escape(xinn_DLOS,x,pos_angles) ;
						escape(xinn_DLOS,ndir,dir_angles) ;

						if ( pos_angles[0] <= dOmega_source / 4. && pos_angles[1] <= dOmega_source / 4. 
							&& dir_angles[0] <= dOmega_collect / 4. && dir_angles[1] <= dOmega_collect / 4. ) {  

				    		if ( path/rout>=tmin && path/rout<=tmax ) {
 								
 								sum_cs_avg[iwave] += 1 ;

	 							// Save Stokes parameters
								if (dopolar) {

									stokesI[iphot] = I;
									stokesQ[iphot] = Q;
									stokesU[iphot] = U;

									avgI += I;
									avgQ += Q;
									avgU += U;
								}
							}


							bin_index = (1 - ndir[2]) / step_costh ;
							if (bin_index <= thetabin) sum_cs[bin_index * nwavel + iwave] += 1 ;
							else printf("Packet escaped but was not collected in an angular bin \n");

							// You already went back. Now increment path (s+d in the notes)
							path += finalpath(&istep,x,ndir,axis_out);
		    	
				    	}

				    	// To ensure that packets that are not selected have fractime = -2 (see below)
				    	else {

				    		sum_away_avg[iwave] += 1 ;
				    		path = - rout * 2;
				    	}

					}

					// Packet interacts in the dust region
					else { 

						gran = gsl_rng_uniform(rng);
						
						// Absorption
						if (gran > draine[iwave].albedo) {

							nextphoton = 1 ;

							sum_abs_avg[iwave] += 1 ;
							sum_abs[iwave] += 1 ;

							if (saveangles) tottheta = 0.0;
							if (savefirstangle) thetafirst = 0.0 ; 

							// To ensure that absorbed packets have fractime = -1 (see below)
							path = -rout;
						
						}
						
						// Dust scattering
						else {

					   		istep++;
							nextphoton = 0 ;

							costh = cthenyey_greenstein(g_dust); 
					   		th = acos(costh);
							phi = 2 * PI * gsl_rng_uniform(rng);  // phi

					   		// Increment theta
					   		if (saveangles) tottheta += fabs(th);
					   		// First angle if istep==1
					   		if (istep==1 && savefirstangle) thetafirst = th;

					   		// New direction
					   		newdir(&th,&phi,ndir,ndir_out);

					   		// Increment path
					   		path += s_tau ;

					   		// Update Stokes parameters
							if (dopolar) {
	     						
	     						pvalue = draine_pvalue(polar,draine[iwave].lambda,th*180./PI);

								if (pvalue == -1) 
									fprintf(stderr,"Warning: failed to calculate P for %2f microns and %1f deg!\n",draine[iwave].lambda,th*180./PI);
	      
					   			update_stokes(ndir,ndir_out,&costh,&g_dust,&pvalue,&I,&Q,&U);				   			
					   		}

					   		// Update new direction
					   		ndir[0] = ndir_out[0];
					   		ndir[1] = ndir_out[1];
					   		ndir[2] = ndir_out[2];

						}	
					}
				} 

				if (istep==nmax) 
					fprintf(stderr,"Warning: you have reached the maximum number of steps \n");
			
			} 

			// Travel time for the escaping packet
			fractime[iphot] = (float) path / rout;
			avt[iwave] += fractime[iphot];
			//if (fractime[iphot]<0 && fractime[iphot]>-2) printf("%g \n",fractime[iphot]);
			//printf("%g \n",fractime[iphot]);

			// Total scattering angle for the escaping packet
			if (saveangles) totangle[iphot] = (float) tottheta;

			// Total number of scatterings for the escaping packet
			if (savescatter) nscatter[iphot] = (unsigned short) istep;

			// Total number of scatterings for the escaping packet
			if (savefirstangle) firstangle[iphot] = (float) thetafirst;

			// Save packet infos: delay times, totangles and nscatt
			if (save_pkt_info){  // && fractime[iphot]>=1) {

				if (fractime[iphot]>0) {

					if (is == 1) {
						//fprintf(packet_info,"%g %g %d %g %g %g %g %g \n",
						//	draine[iwave].lambda,fractime[iphot],nscatter[iphot],tau_is[iwave],pos_angles[0],pos_angles[1],dir_angles[0],dir_angles[1]);
						fprintf(packet_info,"%g %g %d %g \n",
							draine[iwave].lambda,fractime[iphot],nscatter[iphot],tau_is[iwave]);
					}

					else {
						//fprintf(packet_info,"%g %g %d %g %g %g %g \n",
						//	draine[iwave].lambda,fractime[iphot],nscatter[iphot],pos_angles[0],pos_angles[1],dir_angles[0],dir_angles[1]);
						fprintf(packet_info,"%g %g %d \n",
							draine[iwave].lambda,fractime[iphot],nscatter[iphot]);
					}

				}
			}

		}  // -> new photon

	// Calculate average time for each wavelength
	avt[iwave] = avt[iwave]/sum_cs_avg[iwave];

	if (dopolar == 1) {

		avgQ = avgQ / avgI ;
		avgU = avgU / avgI ;
		avgI = avgI / avgI ;

		printf("%g %g %g \n",draine[iwave].lambda,avgQ*100,avgU*100);

	}


	ttmp = clock();
	time_spent = (double)(ttmp - begin) / CLOCKS_PER_SEC;
	printf("%d out of %d wavelenghts (%g) have been processed. Time: %3.2f \n",iwave+1,nwavel,draine[iwave].lambda,time_spent);

	} // -> new wavelength

				   		
	// Calculate E(B-V)

	ebv_cs =   -2.5*log10(sum_cs_avg[B_i]/((double) nph * norm_factor )) + 
				2.5*log10(sum_cs_avg[V_i]/((double) nph * norm_factor ));
	ebv_cs_err = (2.5/log(10.)) * sqrt( 1./sum_cs_avg[B_i] + 1./sum_cs_avg[V_i] );

	ebv_tot =   -2.5*log10(sum_cs_avg[B_i] * exp(-tau_is[B_i]) /((double) nph * norm_factor )) + 
				2.5*log10(sum_cs_avg[V_i] * exp(-tau_is[V_i]) /((double) nph * norm_factor ));
	ebv_tot_err = (2.5/log(10.)) * sqrt( 1./(sum_cs_avg[B_i] * exp(-tau_is[B_i]) ) + 1./ (sum_cs_avg[V_i] * exp(-tau_is[B_i]) ) );
	
	ebv_unsc =   -2.5*log10(sum_unsc_avg[B_i]/((double) nph * norm_factor )) + 
				2.5*log10(sum_unsc_avg[V_i]/((double) nph * norm_factor ));
	ebv_unsc_err = (2.5/log(10.)) * sqrt( 1./sum_unsc_avg[B_i] + 1./sum_unsc_avg[V_i] );
	
	// ------------------------------------------------
	// ------------------ OUTPUT FILES ----------------
	// ------------------------------------------------


	int print_flag = 0 ;
	
	if (print_flag == 1) {

		scstest = fopen("test_scattering_LMC_CS.out","w");
		ccstest = fopen("test_color_LMC_CS.out","w");     // Output file for CS dust
		cistest = fopen("test_color_LMC_DLOS.out","w");     // Output for IS sanity

		fprintf(scstest,"EBV lambda tau_s tau_a albedo costhe time\n");
	}
	

	for (k=nwavel-1; k>=0; k--) {

		a_x_cs_avg[k] = -2.5*log10(sum_cs_avg[k]/((double) nph * norm_factor  ));
		a_x_tot_avg[k] = -2.5*log10(sum_cs_avg[k] * exp(-tau_is[k]) /((double) nph * norm_factor  ));
		a_x_unsc_avg[k] = -2.5*log10(sum_unsc_avg[k]/( (double) nph * norm_factor ));
		a_x_is_avg[k] = -2.5*log10( exp(-tau_is[k]) );
		a_x_away_avg[k] = -2.5*log10(sum_away_avg[k]/((double) nph ));
		a_x_abs_avg[k] = -2.5*log10(sum_abs_avg[k]/((double) nph ));
		a_x_ins_avg[k] = -2.5*log10(sum_ins_avg[k]/((double) nph ));
		
		if (print_flag == 1) fprintf(scstest,"%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
				ebv_cs,draine[k].lambda,rout/scatdust[k],rout/absdust[k],
				draine[k].albedo,draine[k].cost,avt[k]);

	}


	printf("\nInitial E(B-V)_CS = %.2f E(B-V)_IS = %.2f\n\n",ebv,ebv_is);

	printf("DLOS: E(B-V)_DLOS = %.2f (%.2f), RV_DLOS = %.3f\n\n",ebv_unsc,ebv_unsc_err,a_x_unsc_avg[V_i]/ebv_unsc);

	printf("CS: E(B-V)_CS = %.2f (%.2f), RV_CS = %.2f \n",ebv_cs,ebv_cs_err,a_x_cs_avg[V_i]/ebv_cs);
	printf("IS: E(B-V)_IS = %.2f, RV_IS = %.2f \n",ebv_is,a_x_is_avg[V_i]/ebv_is);
	printf("CS+IS: E(B-V)_TOT = %.2f (%.2f), RV_TOT = %.2f \n",ebv_tot,ebv_tot_err,a_x_tot_avg[V_i]/ebv_tot);
	
	if (print_flag == 1) {

		fprintf(ccstest,"# Summary of colors \n");
		fprintf(ccstest,"# E(B-V) = %6.3f\n# R_V = %6.3f\n#\n",ebv_cs,a_x_cs_avg[V_i]/ebv_cs);
		fprintf(ccstest,"# lambda (microns) :\n# A_X (mag):\n");

		fprintf(cistest,"# Summary of colors \n");
		fprintf(cistest,"# E(B-V) = %6.3f\n# R_V = %6.3f\n#\n",ebv_unsc,a_x_unsc_avg[V_i]/ebv_unsc);
		fprintf(cistest,"# lambda (microns) :\n# A_X (mag):\n");


		for (k=nwavel-1; k>=0; k--) {

			fprintf(ccstest,"%6.4f %6.3f %6.3f %6.3f %6.3f\n",draine[k].lambda,a_x_cs_avg[k],a_x_away_avg[k],a_x_abs_avg[k],a_x_ins_avg[k]);
			fprintf(cistest,"%6.4f %6.3f\n",draine[k].lambda,a_x_unsc_avg[k]);

		}
		


		for (i=0;i<totbin;i++) {
			
			for (k=nwavel-1; k>=0; k--) {
				
				a_x_cs[i*nwavel + k] = -2.5*log10(sum_cs[i*nwavel + k]/((double) nph * step_costh/2. ));
				a_x_tot[i*nwavel + k] = -2.5*log10(sum_cs[i*nwavel + k] * exp(-tau_is[k])/((double) nph * step_costh/2. ));
				a_x_unsc[i*nwavel + k] = -2.5*log10(sum_unsc[i*nwavel + k]/((double) nph * step_costh/2. ));
				a_x_is[i*nwavel + k] = -2.5*log10( exp(-tau_is[k]) );
				a_x_away[i*nwavel + k] = -2.5*log10(sum_away[i*nwavel + k]/((double) nph * step_costh/2. ));
				a_x_abs[i*nwavel + k] = -2.5*log10(sum_abs[i*nwavel + k]/((double) nph * step_costh/2. ));
				a_x_ins[i*nwavel + k] = -2.5*log10(sum_abs[i*nwavel + k]/((double) nph * step_costh/2. ));
				
				fprintf(ccstest,"%6.4f %6.3f %6.3f %6.3f %6.3f\n",draine[k].lambda,a_x_cs[i*nwavel + k],a_x_away[i*nwavel + k],a_x_abs[i*nwavel + k],a_x_ins[i*nwavel + k]);
				fprintf(cistest,"%6.4f %6.3f\n",draine[k].lambda,a_x_unsc[i*nwavel + k]);

				//printf("%d %d\n",sum_cs_avg[i*nwavel + k],sum_cs[i*nwavel + k]);
				//if (i==0) printf("%d %g \n",sum_cs[i*nwavel+k],(double) nph * step_costh/2.);
				//if (k==0) printf("%d %g \n",sum_cs[i*nwavel+k],(double) nph * step_costh/2.);

			}



		} 
		
		fprintf(ccstest,"\n");
		fprintf(cistest,"\n");
		
		fclose(ccstest);
		fclose(cistest);
		fclose(scstest);

	}

	//fprintf(packet_info,"%g \t\t # CS E(B-V) \n",ebv_cs);
	//fprintf(packet_info,"%g \t\t # CS RV \n",a_x_cs_avg[V_i]/ebv_cs);
	//fprintf(packet_info,"%g \t\t # IS E(B-V) \n",ebv_unsc);
	//fprintf(packet_info,"%g \t\t # IS RV \n",a_x_unsc_avg[V_i]/ebv_unsc);

	fclose(packet_info);

	// ------------------------------------------------
	// -------- FREE OUT THE MEMORY ALLOCATED ---------
	// ------------------------------------------------

	free(draine);
  	free(scatdust);
	free(absdust);
	free(avt);
	free(tau_is);
	free(a_x_unsc);
	free(a_x_cs);
	free(a_x_is);
	free(a_x_tot);
	free(a_x_away);
	free(a_x_abs);
	free(a_x_ins);
	free(a_x_unsc_avg);
	free(a_x_cs_avg);
	free(a_x_is_avg);
	free(a_x_tot_avg);
	free(a_x_away_avg);
	free(a_x_abs_avg);
	free(a_x_ins_avg);
	free(sum_cs);
	free(sum_unsc);
	free(sum_away);
	free(sum_abs);
	free(sum_ins);
	free(sum_cs_avg);
	free(sum_unsc_avg);
	free(sum_away_avg);
	free(sum_abs_avg);
	free(sum_ins_avg);
	free(fractime);
	if (saveangles) free(totangle);
	if (savescatter) free(nscatter);
	if (savefirstangle) free(firstangle);
	if (polar != NULL) free_draine_polar(polar);
	if (dopolar) {
		free(stokesI);
		free(stokesQ);
		free(stokesU);
	}


    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nTime: %3.3f \n\n",time_spent);


    return 0;
    
}








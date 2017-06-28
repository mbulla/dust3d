#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "draine.h"


// Read the extinction and scattering table

scatext* init_draine_extscat(char *dusttype,float minlambda, float maxlambda,
			unsigned int resolution, unsigned int *nwavel) {
 
	char *fname = NULL, *drainedir = NULL;

	drainedir = getenv("DRAINEPATH");
	if (drainedir == NULL) drainedir = ".";

	// Scatter and Extinction parameters
	
	const char *lmcextscatf = "draine/kext_albedo_WD_LMCavg_20";
	const char *smcextscatf = "draine/kext_albedo_WD_SMCbar_0";
	const char *mw3extscatf = "draine/kext_albedo_WD_MW_3.1_60_D03.all";
	const char *amcextscatf = "bladh/opt_amc.dat";
	const char *forextscatf = "bladh/opt_for.dat";
	const char *mixextscatf = "bladh/opt_mix.dat";

	FILE *f;
	unsigned int veclen = 25, i = 0;
	
	const unsigned int buf = 128;
	char line[buf];
	
	float twavel, talbedo, tcost, tC_ext, tabs_K, tcos2t;

	// Variables if you want to select only U, B, V, R and I
	double eff_wave[5],diff[5];
	int j;

	scatext *scatextlist = NULL;
	scatextlist = (scatext*) malloc(sizeof(scatext)*veclen);

	// Determine dust type
	//
	fname = (char*) calloc(strlen(drainedir)+255,sizeof(char));
	strcat(fname,drainedir);
	strcat(fname,"/");
	if (strcmp(dusttype,"LMC") == 0)
		strcat(fname,lmcextscatf);
	else if (strcmp(dusttype,"SMC") == 0)
		strcat(fname,smcextscatf);
	else if (strcmp(dusttype,"MW3") == 0)
		strcat(fname,mw3extscatf);
	else if (strcmp(dusttype,"AMC") == 0)
		strcat(fname,amcextscatf);
	else if (strcmp(dusttype, "FOR") ==0)
		strcat(fname, forextscatf);
	else if (strcmp(dusttype,"MIX") == 0)
		strcat(fname,mixextscatf);
	else {
		printf("Warning: unknown dust type!\n");
		fname = NULL;
	}

	f = NULL;
	if (fname != NULL) f = fopen(fname, "rt");

	if (f != NULL) {
		*nwavel = 0;
		i = 0;
		while (fgets(line,buf,f) != NULL) {
			if (line[0] != '#') {
				if (strcmp(dusttype,"LMC") == 0)
					sscanf(line,"%11e %6f  %6f %9e %9e", 
						&twavel, &talbedo, &tcost, &tC_ext, &tabs_K);
				// sscanf(line,"%11e %6f  %6f %9e %9e %7f", 
				// &twavel, &talbedo, &tcost, &tC_ext, &tabs_K, &tcos2t);
				  //sscanf(line,"%9e %6e  %6e %9e %9e %9e",  //should work but doesn't
				//	 &twavel, &talbedo, &tcost, &tC_ext, &tabs_K);
				else if (strcmp(dusttype,"SMC") == 0)
				    sscanf(line,"%11e %6f  %6f %9e %9e %7f", 
				    	&twavel, &talbedo, &tcost, &tC_ext, &tabs_K, &tcos2t);
				  //sscanf(line,"%9e %6f  %6f %9e %9e %7f", 
					// &twavel, &talbedo, &tcost, &tC_ext, &tabs_K);
				else if (strcmp(dusttype,"MW3") == 0)
					sscanf(line,"%11e %6f  %6f %9e %9e %7f", 
						&twavel, &talbedo, &tcost, &tC_ext, &tabs_K, &tcos2t);
				else if (strcmp(dusttype,"AMC") == 0)
					sscanf(line,"  %12e  %12e  %12e  %12e  %12e", 
						&twavel, &talbedo, &tcost, &tC_ext, &tabs_K);
				else if (strcmp(dusttype,"FOR") == 0)
					sscanf(line,"  %12e  %12e  %12e  %12e  %12e", 
						&twavel, &talbedo, &tcost, &tC_ext, &tabs_K);
				else if (strcmp(dusttype,"MIX") == 0)
					sscanf(line,"  %12e  %12e  %12e  %12e  %12e", 
						&twavel, &talbedo, &tcost, &tC_ext, &tabs_K);
				else {
					printf("Warning: unknown dust type!\n");
					fname = NULL;
				}
	
				//if (twavel >= minlambda && twavel <= maxlambda && (i++ % resolution) == 0) {

				//	if ((*nwavel) == veclen) {  // Increase size of vector
				//		veclen *= 2;
				//		scatextlist = realloc(scatextlist,sizeof(scatext)*veclen);
				//	}
	
				//	printf("%d %g \n",*nwavel,);
	
				//	if ((*nwavel) == veclen) {  // Increase size of vector
				//		veclen *= 2;
				//		scatextlist = realloc(scatextlist,sizeof(scatext)*veclen);
				//	}
	
				//	scatextlist[*nwavel].lambda = twavel;
				//	scatextlist[*nwavel].albedo = talbedo;
				//	scatextlist[*nwavel].cost   = tcost;
				//	scatextlist[*nwavel].c_ext  = tC_ext;
				//	scatextlist[*nwavel].abs_k  = tabs_K;
				//	scatextlist[*nwavel].cos2t  = tcos2t;
	
				//	(*nwavel)++;
				//}



				// Select only U, B, V, R and I
				
				*nwavel = 5;

				eff_wave[0] = 0.365;// 0.442;//0.365;
				eff_wave[1] = 0.442;// 0.547;//0.442;
				eff_wave[2] = 0.546;// 1.22; //0.546;
				eff_wave[3] = 0.658;// 1.63; //0.658;
				eff_wave[4] = 0.806;// 2.19; //0.806;
						
				if (twavel >= minlambda && twavel <= maxlambda && (i++ % resolution) == 0) {

					if (i==1) {

						diff[0] = 1e10;
						diff[1] = 1e10;
						diff[2] = 1e10;
						diff[3] = 1e10;
						diff[4] = 1e10;
					}

					if ((*nwavel) == veclen) {  // Increase size of vector
						veclen *= 2;
						scatextlist = realloc(scatextlist,sizeof(scatext)*veclen);
					}

					for (j=0;j<(*nwavel);j++) {

						if ( fabs(twavel-eff_wave[j]) < diff[j] ) {

							scatextlist[j].lambda = twavel;
							scatextlist[j].albedo = talbedo;
							scatextlist[j].cost   = tcost;
							scatextlist[j].c_ext  = tC_ext;
							scatextlist[j].abs_k  = tabs_K;
							scatextlist[j].cos2t  = tcos2t;

							diff[j] = fabs(twavel-eff_wave[j]) ;

						}

					}
						
				}

	    	}
	
	  	}

	fclose(f);
	
	} 

	else {
		fprintf(stderr,"Warning: failed to open file extinction file!\n");
		free(scatextlist);
		scatextlist = NULL;
	}

	return scatextlist;
}


// Find the index for the wavelength closest to the given wavelength

void find_draine_lambda_idx(scatext *list, unsigned int len, float lambda, int *idx) {
	int i;
	float dlambda = 1e9;

	for (i=0; i<len; i++) {
		if (fabs(list[i].lambda - lambda) < dlambda) {
			dlambda = fabs(list[i].lambda - lambda);
			*idx = i;
		}
	}
}


// Read polarisation table
//
polarisation* init_draine_pol(char *dusttype) {
	char *tablef = NULL, *anglesf = NULL, *lambdaf = NULL;
	char *drainedir;
	const unsigned int buf = 2048;
	char line[buf];
	int nlambda, nangles;
	unsigned int maxnlambda, maxnangles;
	int col, row, i;
	float *angles = NULL, *lambda = NULL;
	float **table = NULL;
	float min, max;

	FILE *f;

	char *current;

	polarisation *p = NULL;

	// Polarisation table
	//
	const char *lmcpolpf = "draine/callscat_init_p_LMC_avg.dat";
	const char *smcpolpf = "draine/callscat_init_p_SMC_bar.dat";
	const char *mw3polpf = "draine/callscat_init_p_MW_3.1.dat";
	
	
	// Wavelength vector for polarisation table
	//
	const char *lmcpollf = "draine/callscat_init_p_LMC_avg_lambda.dat";
	const char *smcpollf = "draine/callscat_init_p_SMC_bar_lambda.dat";
	const char *mw3pollf = "draine/callscat_init_p_MW_3.1_lambda.dat";
	
	
	// Angle vector for polarisation table
	
	const char *lmcpolaf = "draine/callscat_init_p_LMC_avg_angle.dat";
	const char *smcpolaf = "draine/callscat_init_p_SMC_bar_angle.dat";
	const char *mw3polaf = "draine/callscat_init_p_MW_3.1_angle.dat";
	

	
	polarisation* failed(char *tablef, char *lambdaf, char *anglesf, polarisation* p,
			float *angles, float *lambda, float **table) {
		if (tablef  != NULL) free(tablef);
		if (lambdaf != NULL) free(lambdaf);
		if (anglesf != NULL) free(anglesf);
		if (p       != NULL) free(p);
		if (angles  != NULL) free(angles);
		if (lambda  != NULL) free(lambda);
		if (table   != NULL) free(table);
		return NULL;
	}

	// Read angle or lambda vector
	//
	void readvector(const char *filename, float *vec, int *vlen, unsigned int *maxvlen) {
		FILE *f;
		const unsigned int buf = 256;
		char line[buf];
		float val;
		
		*vlen = 0;

		f = fopen(filename, "rt");
		if ( f != NULL ) {    
			while (fgets(line,buf,f) != NULL) {
				if( line[0] != '#' ) {            // Comment
	  				sscanf(line,"%f", &val);
					/* This doesn't work since we are overwriting the address of vec
					   if ((*vlen) == (*maxvlen))      // Allocate more memory if needed
					   vec = realloc(vec,sizeof(float) * 2 * (*maxvlen));
					*/
				  	if ((*vlen) == (*maxvlen)) {
						printf("Error: oops this is where we should have dynamic mem allocation!\n");
						exit(1);
					}
					vec[(*vlen)++] = val;
				}
			}
			fclose(f);
		} 

		else *vlen = -1;               // Failed to read file
	    
	}
	  
	drainedir = getenv("DRAINEPATH");
	if (drainedir == NULL) drainedir = ".";

	tablef  = (char*) calloc(strlen(drainedir)+255,sizeof(char));
	strcat(tablef,drainedir);
	strcat(tablef,"/");

	lambdaf = (char*) calloc(strlen(drainedir)+255,sizeof(char));
	strcat(lambdaf,drainedir);
	strcat(lambdaf,"/");

	anglesf = (char*) calloc(strlen(drainedir)+255,sizeof(char));
	strcat(anglesf,drainedir);
	strcat(anglesf,"/");
	
	// Determine dust type
	//
	if (strcmp(dusttype,"LMC") == 0) {
		strcat(tablef,lmcpolpf);
		strcat(lambdaf,lmcpollf);
		strcat(anglesf,lmcpolaf);
	} else if (strcmp(dusttype,"SMC") == 0) {
		strcat(tablef,smcpolpf);
		strcat(lambdaf,smcpollf);
		strcat(anglesf,smcpolaf);
	} else if (strcmp(dusttype,"MW3") == 0) {
		strcat(tablef,mw3polpf);
		strcat(lambdaf,mw3pollf);
		strcat(anglesf,mw3polaf);
	} 
	else if (strcmp(dusttype,"AMC") == 0) {
	  
	}
	else {
		fprintf(stderr,"Error: unknown dusttype '%s'!\n",dusttype);
		return failed(tablef,lambdaf,anglesf,p,angles,lambda,table);
	}

	p = (polarisation*) malloc(sizeof(polarisation));
	
	// Allocate memory
	//
	maxnlambda = 256;
	maxnangles = 256;
	lambda = (float*) malloc(sizeof(float)*maxnlambda);
	angles = (float*) malloc(sizeof(float)*maxnangles);
	
	// Read lambda and angle vectors
	//
	readvector(lambdaf,lambda,&nlambda,&maxnlambda);
	if (nlambda < 0) {               // Failed to read vectors
		fprintf(stderr,"Error: failed to read lambda vector '%s'\n!",lambdaf);
		return failed(tablef,lambdaf,anglesf,p,angles,lambda,table);
		exit(1);
	}
	readvector(anglesf,angles,&nangles,&maxnangles); 
	if (nlambda < 0) {               // Failed to read vectors
		fprintf(stderr,"Error: failed to read angles vector '%s'\n!",anglesf);
		return failed(tablef,lambdaf,anglesf,p,angles,lambda,table);
		exit(1);
	}

	
	p->lambda  = lambda;
	p->nlambda = nlambda;
	p->angles  = angles;
	p->nangles = nangles;
	
	// Find the min and max values of the vectors
	//
	min = lambda[0]; max = lambda[0];
	for (i=0; i<nlambda; i++) {
		if (min > lambda[i]) min = lambda[i];
		if (max < lambda[i]) max = lambda[i];
	}
	p->minlambda = min;
	p->maxlambda = max;
	
	min = angles[0]; max = angles[0];
	for (i=0; i<nangles; i++) {
		if (min > angles[i]) min = angles[i];
		if (max < angles[i]) max = angles[i];
	}
	p->minangles = min;
	p->maxangles = max;

	// Read the table
	//
	row = col = 0;
	f = fopen(tablef, "rt");
	if (f == NULL) {
		fprintf(stderr,"Error: failed to read polarisation table '%s'\n!",tablef);
		return failed(tablef,lambdaf,anglesf,p,angles,lambda,table);
	}

	// Allocate memory for table
	//
	table = (float **) malloc(sizeof(float *) * nlambda);
	for (i = 0; i < nlambda; i++) table[i] = (float *) malloc(sizeof(float) * nangles);
	
	while(fgets(line,buf,f) != NULL && row < nangles ) {
		if (line[0] != '#') {
			current = strtok(line," ");
			col = 0;
			while (current != NULL && col < nlambda) {
				if (current[0] != ' ' && current[0] != '\n') table[col++][row] = (float) atof(current);

				current = strtok(NULL, " ");
	    	}
	    	row++;
		}
	} // done reading the table
	
	p->table = table;
	fclose(f);

	return p;
}


// Free up the memory used by the polarisation structure
//
void free_draine_polar(polarisation *p) {
	unsigned int i;
	for (i = 0; i < p->nlambda; i++) free(p->table[i]);
	
	free(p->table);
	free(p->lambda);
	free(p->angles);
}


// Calculate the p-value through bilinear interpolation for
// a given wavlength (microns) and scattering angle (degrees);
//
float draine_pvalue(polarisation *p, float l, float a) {
	float pvalue = -1.0;
	float l1, l2, a1, a2;
	float f11, f12, f21, f22;
	float d;
	unsigned int n;
	unsigned int i,j;

	if (l >= p->minlambda && l <= p->maxlambda &&
		  a >= p->minangles && a <= p->maxangles ) {

		// Find indices for grid that surrounds the given point
		//
		n = 1;
		while (p->lambda[n] < l && n < p->nlambda-1) n++;
		i = n-1;
		l1 = p->lambda[i];
		l2 = p->lambda[i+1];

		n = 1;
		while (p->angles[n] < a && n < p->nangles-1) n++;
		j = n-1;
		a1 = p->angles[j];
		a2 = p->angles[j+1];

		// Bilinear interpolation
		//
		f11 = p->table[i][j];
		f12 = p->table[i][j+1];
		f21 = p->table[i+1][j];
		f22 = p->table[i+1][j+1];
		
		d = (l2-l1) * (a2-a1);
		
		pvalue =
			(l2-l)*(a2-a) * f11/d +
			(l-l1)*(a2-a) * f21/d +
			(l2-l)*(a-a1) * f12/d +
			(l-l1)*(a-a1) * f22/d;
	}

	return pvalue;
}



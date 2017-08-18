#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "mylib.h"
#include "draine.h"

// Distance to a spherical surface
double dist_to_sphere(double *rsph, double *x0, double *ndir) {

	double r0_sq,rsph_sq,x0_dot_ndir,s;

	r0_sq = pow(x0[0],2.) + pow(x0[1],2.) + pow(x0[2],2.);
	rsph_sq = pow(*rsph,2.);
	x0_dot_ndir = dot(x0,ndir);

	s = -x0_dot_ndir + sqrt( pow(x0_dot_ndir,2.) + rsph_sq - r0_sq );	

	return s;

}

// Distance to a ellipsoid surface
double dist_to_ellipsoid(double *axis, double *x0, double *ndir) {

	double xell[3],nell[3],A,B,C,s;

	xell[0] = x0[0] / axis[0] ;
	xell[1] = x0[1] / axis[1] ;
	xell[2] = x0[2] / axis[2] ;

	nell[0] = ndir[0] / axis[0] ;
	nell[1] = ndir[1] / axis[1] ;
	nell[2] = ndir[2] / axis[2] ;

	A = pow(nell[0],2.) + pow(nell[1],2.) + pow(nell[2],2.) ;
	B = 2 * dot(xell,nell) ;
	C = pow(xell[0],2.) + pow(xell[1],2.) + pow(xell[2],2.) - 1 ;

	s = ( - B + sqrt( pow(B,2.) - 4 * A * C ) ) / 2. / A ;	

	if (s<0) 
		fprintf(stderr,"Error: negative distance.");
	if ( (axis[0]>axis[2] && s>2*axis[0]) || (axis[0]<axis[2] && s>2*axis[2]) )
		fprintf(stderr,"Error: distance to the ellipsoidal surface is too large \n");

	return s;

}


// Source emission
void source(double *dOmega_source, double *axis, double *x, double *ndir, double *xinn_DLOS) {

	double u,v,s,rfinal,nobs[3],dphi,dcosth;
	//double K,sint2;

	// -------------- POSITION --------------

	/* Create at the center */
	
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;	

    xinn_DLOS[0] = 0;
    xinn_DLOS[1] = 0;
    xinn_DLOS[2] = 0;

	// -------------- DIRECTION  (see notes 23-25/01/17 and 18/08/17) -------------
	
	dphi = ( *dOmega_source / 8. ) ;
	dcosth = sin( dphi );

	// Random direction
    u = 2 * dcosth * gsl_rng_uniform(rng) - dcosth ;
	if (gsl_rng_uniform(rng)<0.5) v = 2 * dphi * gsl_rng_uniform(rng) - dphi ;
	else v = PI + 2 * dphi * gsl_rng_uniform(rng) - dphi ;

	//u = 0;
	//v = PI/2.;

	ndir[0] = sqrt(1-pow(u,2.)) * cos(v);
	ndir[1] = sqrt(1-pow(u,2.)) * sin(v);
	ndir[2] = u;

	nobs[0] = 1 ;
	nobs[1] = 0 ;
	nobs[2] = 0 ;

	// ------- MOVE PACKET TO INNER SURFACE ---------
	
    if (axis[0]==0 && axis[1]==0 && axis[2]==0) s = 0 ;

    else {

	    s = dist_to_ellipsoid(axis,x,ndir); 

	    x[0] += s * ndir[0];
	    x[1] += s * ndir[1];
	    x[2] += s * ndir[2];

	    s = dist_to_ellipsoid(axis,xinn_DLOS,nobs); 

	    xinn_DLOS[0] += s * nobs[0];
	    xinn_DLOS[1] += s * nobs[1];
	    xinn_DLOS[2] += s * nobs[2];

	    /* Check you are at rinn */
	    rfinal = sqrt(pow(x[0]/axis[0],2.) + pow(x[1]/axis[1],2.) + pow(x[2]/axis[2],2.));
		if (fabs(rfinal - 1)>1e-4) fprintf(stderr,"Error in source function \n");
  	} 	

}

/* Calculate last step before leaving the sphere (s and d in Ariel/Rahman notes) */
/* See my notes (18/11/16 and 11/01/17) for calculations */
double finalpath(int *istep, double *x, double *ndir, double *axis_out) {

	double s,d;
	double xbound[3];
	double costh,totpath;
	double origin[3],sdir,sbound;

	// Calculate s to the boundary
	s = dist_to_ellipsoid(axis_out,x,ndir);

	if (s<0) printf("Warning: s<0 \n");

	// Calculate d

	if (*istep==0) d = 0 ; // Packets escaping after creation -> we know already that s=rout-rinn and d=0;

	else {

		origin[0] = 0. ;
		origin[1] = 0. ;
		origin[2] = 0. ;

		sdir = dist_to_ellipsoid(axis_out,origin,ndir) ;

		xbound[0] = x[0] + s * ndir[0];
		xbound[1] = x[1] + s * ndir[1];
		xbound[2] = x[2] + s * ndir[2];

		sbound = sqrt( dot(xbound,xbound) );

		costh = norm_plus_dot(xbound,ndir);

		d = sdir - costh * sbound ;	

		//if (fabs(d) > 2*axis_out[0]) { 

		//	printf("Error: d/r = %.2f > 1 should not happen!\n",d/axis_out[0]);
		//	printf("%g %g %g \n",xbound[0],xbound[1],xbound[2]);	
		//	printf("%g %g %g %g %g %g \n",s,sdir,costh,sbound,d,axis_out[0]);	
		//}

	}

	totpath = s + d;

	return totpath;

}



// Calculate the continuum optical depth s taking into account different regions with different densities

double tau_to_s(double *x, double *xinn_DLOS, double *ndir, double *rinn, double *rout, double *tau, double *scatdust, double *absdust, double *dOmega_source) {

	double s=0,eps;
	double tautmp=0,xt[3];
	double angles[2];

	double dens_ratio = 1;
	double ang_dens = (*dOmega_source) ;

	xt[0]=x[0];   xt[1]=x[1];   xt[2]=x[2];

	eps = ( *rout - *rinn ) / 50. ;

	while( tautmp < (*tau) ) {

		xt[0] += eps * ndir[0];
		xt[1] += eps * ndir[1];
		xt[2] += eps * ndir[2];

		escape(xinn_DLOS,x,angles);

		if ( angles[0] <= ang_dens && angles[1] <= ang_dens ) tautmp += eps * ( 1. / (*absdust) + 1. / (*scatdust) ) ;
		else tautmp += dens_ratio * eps * ( 1. / (*absdust) + 1. / (*scatdust) ) ;

		if ( tautmp < (*tau) ) s += eps ;

	}

	escape(xinn_DLOS,x,angles);

	if ( angles[0] <= ang_dens && angles[1] <= ang_dens ) {

		// Remove last contribution - that brought tautmp above tau
		tautmp = tautmp - eps * ( 1. / (*absdust) + 1. / (*scatdust) ) ;
		// Add last path length to reach the event
		s = s + ( (*tau) - tautmp ) / ( 1. / (*absdust) + 1. / (*scatdust) ) ;
	}

	else {
	
		tautmp = tautmp - dens_ratio * eps * ( 1. / (*absdust) + 1. / (*scatdust) ) ;
		s = s + ( (*tau) - tautmp ) / ( 1. / (*absdust) + 1. / (*scatdust) ) / dens_ratio ;
    }
 


	return s;

}

// See notes 23-25/01/17
void escape(double *xinn_DLOS, double *vec, double *angles) {

	double zaxis[3],costh_obs,phi_obs,costh_esc,phi_esc;
	double dtheta,dphi;

	zaxis[0] = 0;
	zaxis[1] = 0;
	zaxis[2] = 1;

	// ------- Angle between escaping position and DLOS -------

	/* d(costh) selection */

	costh_obs = norm_plus_dot(zaxis,xinn_DLOS);
	costh_esc = norm_plus_dot(zaxis,vec);

	/* d(phi) selection */

	phi_obs = atan2(xinn_DLOS[1],xinn_DLOS[0]);
	phi_esc = atan2(vec[1],vec[0]);

	/* Angular distance to the DLOS (see notes 03/02/17) */
	// dcosth = sin(theta), but acos is also OK because I am taking the difference

	dtheta = fabs( asin(costh_obs) - asin(costh_esc) );   
	dphi = fabs(phi_obs - phi_esc) ;

	angles[0] = dtheta ;
	angles[1] = dphi ;

	//double angle = 45*PI/180.;
	//double dcosang = sin(angle);
	//int esc;
//
	//if ( (fabs(costh_obs - costh_esc) < dcosang) && ( fabs(phi_obs - phi_esc) < angle ) ) esc = 1 ;
	//else esc = 0;
//
//if (esc==1 && final_ang>angle ) printf("%g %g  %g %g \n",dtheta,dphi,final_ang,angle);


}



int check(double *x, double *axis) {
    
    double result;
    
    if (axis[0]==0 && axis[1]==0 && axis[2]==0) return 0 ;
    
    else {
    
    	result = pow(x[0]/axis[0],2.) + pow(x[1]/axis[1],2.) + pow(x[2]/axis[2],2.) ;
	   
	    if (result>=1.0) return 0;

	    else return 1; 

    }

    
}


double cthenyey_greenstein(double g) {
	
	double r,gsq,costh;
	
	r = gsl_rng_uniform(rng);
	
	gsq = pow(g,2.);
	
	costh = ( 1 + gsq - pow( (1 - gsq) / (1 - g + 2 * g * r) , 2. ) ) / 2. / g ;
	
	if (fabs(costh) > 1.0) {
		
		printf("Warning !! costh=%f\n",costh);
		
		if (costh > 1.) costh = 1.;
		else if (costh < -1.) costh = -1.;
	
	}
	
	return costh;

}

// This was checked to give the same result as cthenyey_greenstein
double cthenyey_greenstein_rahman(double tau) {  
  double r, rt, ct;
  r  = gsl_rng_uniform(rng);
  rt = r*tau;
  ct = (1. - tau + rt)*r*pow((1. + tau),2)/pow((1. - tau + 2.*rt),2);
  ct = 2.*ct-1.;
  if (fabs(ct) > 1.0) {
    printf("Warning !! ct=%f\n",ct);
    if (ct > 1.) {
      ct = 1.;
    } else {
      if (ct < -1.) ct = -1.;
    }
  }
  return ct;
}


/* ----- New direction after scattering (Kalos & Whitlock 2008) ------------ */
void newdir(double *tsc, double *phisc, double *n_in, double *n_out) {
    
    double a,b,c;
    
    
    if( fabs(n_in[2]) < 0.99999 ) {
        
        a = sin((*tsc))/sqrt(1.-pow(n_in[2],2.)) * ( n_in[1] * sin(*phisc) - n_in[0] * n_in[2] * cos(*phisc) ) + n_in[0] * cos(*tsc) ;
        b = sin(*tsc)/sqrt(1-pow(n_in[2],2.)) * ( - n_in[0] * sin(*phisc) - n_in[1] * n_in[2] * cos(*phisc) ) + n_in[1] * cos(*tsc) ;
        c = sin(*tsc) * cos(*phisc) * sqrt(1-pow(n_in[2],2.))  +  n_in[2] * cos(*tsc) ;
        
    }
    
    else {
        
        a = sin(*tsc) * cos(*phisc) ;
        b = sin(*tsc) * sin(*phisc) ;
        
        if( n_in[2]>0 ) c =  cos(*tsc) ;
        else c = - cos(*tsc) ;
        
    }
    
    n_out[0]=a;
    n_out[1]=b;
    n_out[2]=c;
    
}


/* ----- Update Stokes parameters (similar to toycode in Bulla+2015 but for dust scattering) ------- */

void update_stokes(double *n_in, double *n_out, double *costsc, double *g_dust, float *pvalue, double *I, double *Q, double *U) {

	double ref1[3],ref2[3];
	double i1,i2,cos2i1,sin2i1,cos2i2,sin2i2;
	double P1,P2,P3;
	double gtmp,gtmp_sq;
	double Itmp,Qtmp,Utmp;

	/* 1) Compute meridian frame axes ref1 and ref2 */
	meridian(n_in,ref1,ref2);

	/* 2) Rotate Stokes parameters in the scattering plane --> calculate i1 angle
		This is the i1 angle of Bulla+2015, obtained computing the angle between the
		reference axes ref1 and ref2 in the meridian frame and the corresponding axes
		ref1_sc and ref2_sc in the scattering plane. It is the supplementary angle of the
		scatt angle phisc chosen in the rejection technique */
	i1 = rot_angle(n_in,n_out,ref1,ref2);
	cos2i1 = cos(2 * i1) ;
	sin2i1 = sin(2 * i1) ;
	stokes_rotation_counterclock(I,Q,U,&cos2i1,&sin2i1);

	/* 3) Scattering matrix */

	gtmp = (*g_dust);
	gtmp_sq = pow(gtmp,2.);

	P1 = (1- gtmp_sq) / pow(1 + gtmp_sq - 2 * gtmp * (*costsc),3/2.);
	P2 = - (*pvalue) * P1 * (1 - pow(*costsc,2.) ) / (1 + pow(*costsc,2.) ) ;
	P3 = P1 * 2 * (*costsc) / (1 + pow(*costsc,2.) ) ;

	Itmp = 3./4. * ( P1 * (*I) + P2 * (*Q) );
	Qtmp = 3./4. * ( P2 * (*I) + P1 * (*Q) );
	Utmp = 3./4. * P3 * (*U) ;

	Qtmp = Qtmp / Itmp;
	Utmp = Utmp / Itmp;
	Itmp = Itmp / Itmp;

	/* 4) Compute meridian frame axes ref1 and ref2 relative to outgoing direction*/
	meridian(n_out,ref1,ref2);

	/* 5) Rotate Stokes parameters back to the meridian frame. 
		This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
		reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
		meridian frame. NB: we need to add PI to transform THETA to i2 */
	i2 = PI + rot_angle(n_out,n_in,ref1,ref2);
	cos2i2 = cos(2 * i2) ;
	sin2i2 = sin(2 * i2) ;
	stokes_rotation_clock(&Itmp,&Qtmp,&Utmp,&cos2i2,&sin2i2);


	*I = Itmp;
	*Q = Qtmp;
	*U = Utmp;

}


/* ----------------------- Routine to compute the meridian frame axes ref1 and ref2 ----------------------------------------*/

void meridian(double *n, double *ref1, double *ref2){
    
    
    // for ref_1 use (from triple product rule)
        
    ref1[0] = -1. * n[0] * n[2]/ sqrt(n[0]*n[0] + n[1]*n[1]);
    ref1[1] = -1. * n[1] * n[2]/ sqrt(n[0]*n[0] + n[1]*n[1]);
    ref1[2] = (1 - (n[2] * n[2]))/ sqrt(n[0]*n[0] + n[1]*n[1]);
    
    // for ref_2 use vector product of n_cmf with ref1
        
    ref2[0] = n[2] * ref1[1] - n[1] * ref1[2];
    ref2[1] = n[0] * ref1[2] - n[2] * ref1[0];
    ref2[2] = n[1] * ref1[0] - n[0] * ref1[1];
    
    
}


/* --------------------------------- Rotation angle from the scattering plane ---------------------------------------------*/
/* -------- We need to rotate Stokes Parameters to (or from) the scattering plane from (or to) the meridian frame -------- */
/* ------------------------------- such that Q=1 is in the scattering plane and along ref1 -------------------------------- */

double rot_angle(double *n1, double *n2, double *ref1, double *ref2) {
    
    double ref1_sc[3],cos_stokes_rot_1, cos_stokes_rot_2, i=0;
    
    // ref1_sc is the ref1 axis in the scattering plane ref1 = n1 x ( n1 x n2 )
    ref1_sc[0] = n1[0] * norm_plus_dot(n1,n2) - n2[0];
    ref1_sc[1] = n1[1] * norm_plus_dot(n1,n2) - n2[1];
    ref1_sc[2] = n1[2] * norm_plus_dot(n1,n2) - n2[2];

    cos_stokes_rot_1 = norm_plus_dot(ref1_sc,ref1);
    cos_stokes_rot_2 = norm_plus_dot(ref1_sc,ref2);
    
    if (cos_stokes_rot_1<-1) cos_stokes_rot_1=-1;
    if (cos_stokes_rot_1>1) cos_stokes_rot_1=1;
    
    if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 > 0)) i = acos(cos_stokes_rot_1);
    if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 < 0)) i = 2 * acos(-1.) - acos(cos_stokes_rot_1);
    if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 < 0)) i = acos(-1.) + acos(fabs(cos_stokes_rot_1));
    if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 > 0)) i = acos(-1.) - acos(fabs(cos_stokes_rot_1));
    if (cos_stokes_rot_1 == 0) i = acos(-1.)/2.;
    if (cos_stokes_rot_2 == 0) i = 0.0 ;
    
    if (i!=i ) printf("Warning NaN: %3.6f \t %3.6f \t %3.6f \n",cos_stokes_rot_1,cos_stokes_rot_2,acos(cos_stokes_rot_1));
    
    
    return i;
    
}

/* --------------------------- Stokes Parameter after rotation ------------------------------ */

void stokes_rotation_counterclock(double *I, double *Q, double *U, double *a, double *b) {
    
    double It,Qt,Ut;
    
    It = *I;
    Qt = (*Q) * (*a) - (*U) * (*b);
    Ut = (*Q) * (*b) + (*U) * (*a);

	*I = It;    
    *Q = Qt;
    *U = Ut;
    
}

void stokes_rotation_clock(double *I, double *Q, double *U, double *a, double *b) {
    
    double It,Qt,Ut;
    
    It = *I;
    Qt = (*Q) * (*a) + (*U) * (*b);
    Ut = -(*Q) * (*b) + (*U) * (*a);
    
	*I = It;    
    *Q = Qt;
    *U = Ut;
    
}


double through_inner(double *axis_inn, double *s_tau, double *x, double *ndir) {

	double s_step, xtmp[3];
	int Nsteps,i,inside_inner;
	double path_inside;

	Nsteps = 100;

	path_inside = 0;

	xtmp[0] = x[0];
	xtmp[1] = x[1];
	xtmp[2] = x[2];

	s_step = (*s_tau)/Nsteps;
	
	for (i=0;i<Nsteps;i++) {

		xtmp[0] = xtmp[0] - s_step * ndir[0];
		xtmp[1] = xtmp[1] - s_step * ndir[1];
		xtmp[2] = xtmp[2] - s_step * ndir[2];

		inside_inner = check(xtmp,axis_inn);

		if (inside_inner == 1) path_inside += s_step;

	}

	return path_inside;

}



double norm_plus_dot(double *x, double *y){
    
    double nx,ny,result;

    nx = sqrt( (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]) );
    ny = sqrt( (y[0] * y[0]) + (y[1] * y[1]) + (y[2] * y[2]) );

    result = ( (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]) ) / nx / ny;

    if (result>1.) result = 1.;
    else if (result<-1.) result = -1.;
    
    return result;
    
}


double dot(double *x, double *y){
    
    double result;

    result = (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]) ;
    
    return result;
    
}

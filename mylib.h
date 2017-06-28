const gsl_rng *rng;
#define PI 3.14159265
#define CLIGHT 2.99792458e10
#define ARCSEC_TO_RAD (1 / 3600. * PI / 180.)
#define PC_TO_CM  (3.26 * 365 * 24 * 3600 * CLIGHT)

double dist_to_sphere(double *rsph, double *x, double *ndir);
double dist_to_ellipsoid(double *axis, double *x, double *ndir);
void source(double *sn_dust_dphi, double *axis, double *x, double *ndir, double *xinn_DLOS);
double finalpath(int *istep, double *x, double *ndir, double *axis_out);
double tau_to_s(double *x, double *xfirst, double *ndir, double *rinn, double *rout, double *tau, double *scatdust, double *absdust, double *sn_dust_dphi);
void escape(double *xinn_DLOS, double *vec, double *angles);
int check(double *x, double *axis);
double cthenyey_greenstein(double g);
double cthenyey_greenstein_rahman(double tau);
void newdir(double *tsc, double *phisc, double *n_in, double *n_out);
void update_stokes(double *n_in, double *n_out, double *costsc, double *g_dust, float *pvalue, double *I, double *Q, double *U);
void meridian(double *n, double *ref1, double *ref2);
double rot_angle(double *n1, double *n2, double *ref1, double *ref2);
void stokes_rotation_counterclock(double *I, double *Q, double *U, double *a, double *b);
void stokes_rotation_clock(double *I, double *Q, double *U, double *a, double *b);
double through_inner(double *axis_inn, double *s_tau, double *x, double *ndir);
double norm_plus_dot(double *x, double *y);
double dot(double *x, double *y);
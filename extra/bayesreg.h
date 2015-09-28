/* *
   File:        bayesreg.h
   Author:      Andrew W. Moore
   Created:     Wed Apr 23 15:43:43 EDT 1997
   Description: Header for Bayesian Regression Analysis

   Copyright 1997, Carnegie Mellon University
*/

#ifndef BAYESREG_H
#define BAYESREG_H

#include "./utils/amdmex.h"

/* A regtype denotes one of four kinds of regression.

     LINEAR_REGTYPE is simple old linear regression: Find c and b to fit

          y = c + B^T x


     CIRCLE_REGTYPE has one extra coefficient: mean-of-squares of x's:
         Find c and b and a to fit

          y = c + B^T x + 0.5 * x^T A x where A = Diag(a,a,a ... a)

     ELLIPSE_REGTYPE has a coefficient for each pure square term but
       no coefficient for crossterms: Find c b and vector a to fit

          y = c + b^T x + 0.5 * x^T A x where A = Diag(a)

     QUADRATIC_REGTYPE is full quadratic. Find c, b and symmetric A to fit

          y = c + b^T x + 0.5 * x^T A x
*/

typedef struct brdat
{
  dyv *x;
  dyv *w;
  dyv *y;
  double cx;
  double rx;
  double n0;
  double ss0;
  double sigsqd;
  int regtype;
  dyv *b0;
  dyv *b0v;
  double xlo;
  double xhi;
} brdat;

#define LINEAR_REGTYPE    0
#define CIRCLE_REGTYPE    1
#define ELLIPSE_REGTYPE   2
#define QUADRATIC_REGTYPE 3
#define NUM_REGTYPES      4

/* A regtype denotes one of four kinds of regression.

     LINEAR_REGTYPE is simple old linear regression: Find c and b to fit

          y = c + B^T x


     CIRCLE_REGTYPE has one extra coefficient: mean-of-squares of x's:
         Find c and b and a to fit

          y = c + B^T x + 0.5 * x^T A x where A = Diag(a,a,a ... a)

     ELLIPSE_REGTYPE has a coefficient for each pure square term but
       no coefficient for crossterms: Find c b and vector a to fit

          y = c + b^T x + 0.5 * x^T A x where A = Diag(a)

     QUADRATIC_REGTYPE is full quadratic. Find c, b and symmetric A to fit

          y = c + b^T x + 0.5 * x^T A x
*/


/* my_error if no such regtype */
char *regtype_to_name(int regtype);

int num_regtypes();

/* Return -1 if no such name */
int name_to_regtype(char *regtype_name);

int regtype_to_num_terms(int num_inputs,int regtype);

int regtype_from_args(char *key,int argc,char *argv[],int default_regtype);

/***** Main Bayesian Regression Code ****/

/* rows may be NULL denoting "use all datapoints" */
/* Sets xlo and xhi to define a hyperrectangle surrounding all the
   specified points. Also makes sure that this hyperrectangle is
   small but contains nicely rounded numbers.

   The set of points is the set of rows in x which are mentioned in
   "rows". If rows is NULL, that means "use all points". */
void make_sensible_vector_limits(dym *x,ivec *rows,dyv **xlo,dyv **xhi);


void bayesian_regression(double n,
                         dym *ztwz,
                         dyv *ztwy,
                         double ytwy,
                         dyv *beta_0,
                         dyv *beta_0_variance,
                         double n0,
                         double sigma_squared0,
                         dyv **r_beta_hat,
                         dym **r_vbeta,
                         double *r_nu,
                         double *r_sigma_n_sqd);

double sample_from_gamma(double alpha,double beta);
double sample_from_chi_squared(double nu);
double sample_from_scaled_inverse_chi_squared(double nu,double s_squared);
dyv *mk_vector_of_unit_normal_samples(int size);
dyv *mk_sample_from_normal(dyv *mu,dym *covariance);
dyv *mk_sample_beta(dyv *beta_hat,dym *vbeta,double nu,double sigma_n_sqd);

/* Marginal on specified coefficient is a t-distribution
   with mean *r_mu, variance *r_sigma_squared and dof *r_nu */
void coefficient_marginal(dyv *beta_hat,dym *vbeta,
                          double nu,double sigma_n_sqd,
                          int coeff_index,
                          double *r_nu,double *r_mu,double *r_sigma_squared);
/* Marginal on prediction (at scaled polynomial point z); is a t-distribution
   with mean *r_mu, variance *r_sigma_squared and dof *r_nu */
/*              mean = poi -> beta_hat . z
            variance = z^T poi->vbeta z sigma_n_sqd
                 dof = nu
*/
void predict_marginal(dyv *beta_hat,dym *vbeta,double nu,double sigma_n_sqd,
                      dyv *z,
                      double *r_nu,double *r_mu,double *r_sigma_squared);

double scaled_inverse_chi_squared_cdf_inv(double prob,double nu,
                                          double s_squared);

typedef struct bprior
{
  bool non_informative;
  bool known_sdev;
  double noise_frac_prior;
  double nu0;
  double constant_coeff_prior_sdev;
  double linear_coeff_prior_sdev;
  double quadratic_coeff_prior_sdev;
} bprior;

typedef struct poster
{
  dym *vbeta;
  dyv *beta_hat;
  double nu;   /* = n0 + n */
  double sigma_n_sqd;
} poster;

#define MAX_POSTERS 50

typedef struct reginfo
{
  int regtype;
  dyv *mid;
  dyv *range;
  dyv *ymid;
  dyv *yrange;
} reginfo;

typedef struct posters
{
  int num_posters;
  reginfo *ri;
  poster *ps[MAX_POSTERS];
} posters;

/* In the following , i,j index over elements of input points (x's);
                   they also index over z vectors (vectors of basis function evaluations);

              k indexes over datapoints

              l indexes over outputs.


   Define v[l] as a version of y[l] scaled thusly:   v[l] = ( y[l] - ymid[l] ) / yrange[l]
                                               so    y[l] = ymid[l] + yrange[l] * v[l]

   A z vector corresponding to an x input in a quadratic or linear regression is
   defined shortly (assume x has components x[0] ... x[M-1]);

   Define u[i] = (x[i] - mid[i])/range[i] for i = 0 .. M-1

     z0 = 1
     z1 = u[0]
           .
     z[j] = u[j-1]   for j = 1 ... M
           .
           .
     z[M+1] = u[0] * u[0]
     z[M+2] = u[0] * u[1]
           .
           .
     z[2M]  = u[0] * u[M-1]
     z[2M+1] = u[1] * u[1]
     z[2M+2] = u[1] * u[2]
            .
     z[3M-1] = u[1] * u[M-1]
     z[3M] = u[2] * u[2]
        .
        .
     z[M + M(M+1);/2 - 1] = u[M-2] * u[M-1]
     z[M + M(M+1)/2] = u[M-1] * u[M-1]

   where the quadratic terms exist if and only if regtype!=LINEAR_REGTYPE.

     This means we can interpret the beta vector of coefficients from
     a regression (which will be the same length as z) as


   v = b[0] + ( b[1] )T( u[0]   ) + 0.5 ( u[0]   )T( 2b[M+1] b[M+2]  ... b[2M]         ) ( u[0]   )
              ( b[2] ) ( u[1]   )       ( u[1]   ) ( b[M+2]  b[2M+1] ... b[3M-1]       ) ( u[1]   )
              (    . ) (  .     )       (   .    ) (  .            .       .           ) (  .     )
              (    . ) (  .     )       (   .    ) (  .               .    .           ) (  .     )
              ( b[M] ) ( u[M-1] )       ( u[M-1] ) { b[2M]   b[3M-1]     b[M+M(M+1)/2] ) ( u[M-1] )

   which we write as

     y = CP + BP^T u + 0.5 u^T AP u

   Remember that u = R ( x - m ), v = ( y - ymid ) / yrange

     where R = Diag( 1/range[0] , 1/range[1] ... 1/range[M-1] )
     and   m = ( mid[0]   )
               ( mid[1]   )
               (  .       )
               (  .       )
               ( mid[M-1] )

   So we may write

     y = ym + yr * [ CP + BP^T R ( x - m ) + 0.5 (R(X - m))^T AP R(X - m) ]

   or y = C + B^T x + 0.5 x^T A x where...

     C = ym + yr ( CP - BP^T (Rm) + 0.5 (Rm)^T AP (Rm) )

     B = yr R ( BP - AP(Rm) )

     A = yr R AP R
*/

dyv *mk_u_from_x(reginfo *ri,dyv *x);

dyv *mk_z_from_u(reginfo *ri,dyv *u);
dyv *mk_z_from_u_and_regtype(dyv *u,int regtype);

dyv *mk_z_from_usin(reginfo *ri,dyv *x);

double v_from_y(reginfo *ri,double y,int out_index);
double y_from_v(reginfo *ri,double v,int out_index);

dyv *mk_v_from_usout(reginfo *ri,dyv *y);

void print_cba(char *cname,char *bname,char *aname,
               double c,dyv *b,dym *a);

int reginfo_num_inputs(reginfo *ri);

int reginfo_num_terms(reginfo *ri);

/*** bpriors are used to make posters ***/

bprior *mk_non_informative_bprior();
bprior *mk_known_sdev_non_informative_coeffs_prior(double noise_frac_prior);

bprior *mk_informative_bprior(double nu0,double noise_frac_prior,
                              double constant_coeff_prior_sdev,
                              double linear_coeff_prior_sdev,
                              double quadratic_coeff_prior_sdev);

void free_bprior(bprior *bp);

void fprintf_bprior(FILE *s,char *m1,bprior *bp,char *m2);

bprior *mk_copy_bprior(bprior *bp);

/*** Posters contain the regression information obtained  ****/

reginfo *posters_reginfo(posters *po);

int reginfo_regtype(reginfo *ri);

int posters_num_inputs(posters *po);

void make_cba_from_beta(reginfo *ri,dyv *beta,int out_index,
                        double *c,dyv **b,dym **a);

/* rows may be NULL denotinf "use all datapoints" */
/* Uses only the points that are in "rows" and that are
   not outside the hyper-rectangle
   defined by xlo and xhi.
   xlo and xhi may be NULL denoting "infinitely large hyperrect" */
dym *mk_ztz_from_usins(reginfo *ri,dym *usins,ivec *rows,dyv *xlo,dyv *xhi);

/* rows may be NULL denotinf "use all datapoints" */
/* Uses only the points that are in "rows" and that are
   not outside the hyper-rectangle
   defined by xlo and xhi.
   xlo and xhi may be NULL denoting "infinitely large hyperrect" */
void make_reg_matrices_from_usins_and_usouts(reginfo *ri,dym *usins,
                                             dym *usouts,int out_index,ivec *rows,
                                             dyv *xlo,dyv *xhi,
                                             dyv **ztv,double *vtv);

/*
   The regression is trying to fit

     y[l] = c[l] + b[l]^T x + 0.5 x^T A[l] x

   for each output (where 0 <= l < dym_cols(y))
   (A[l] is always symmetric)

    NOTE THE ORDR IN THE ARGUMENT LIST: xlo ylo xhi yhi. ROOM FOR CONFUSION

   regtype is as defined at top of file

Call with xlo[i] = the minimum value of the i'th component of any input
                   point we're going to see.

Call with ylo[i] = the minimum value of the i'th component of any output
                   point we're going to see.

Call with xhi[i] = the maximum value of the i'th component of any input
                   point we're going to see.

Call with yhi[i] = the maximum value of the i'th component of any output
                   point we're going to see.

And let REGTYPE be one of LINEAR_REGTYPE, CIRCLES_REGTYPE,
ELLIPSES_REGTYPE or QUADRATIC_REGTYPE.
*/
reginfo *mk_reginfo(dyv *xlo,dyv *ylo,
                    dyv *xhi,dyv *yhi,
                    int regtype);


posters *mk_posters_from_reginfo_and_matrices(dym *ztz,dym *ztv,dym *vtv,
                                              reginfo *ri,bprior *bp);

/* rows may be NULL denotinf "use all datapoints" */
posters *mk_posters_from_reginfo(dym *x,dym *y,ivec *rows,reginfo *ri,bprior *bp);

/* This function does a bayesian linear or quadratic regression on
   all or a portion of the data in x and y.

   The regression is trying to fit

     y[l] = c[l] + b[l]^T x + 0.5 x^T A[l] x

   for each output (where 0 <= l < dym_cols(y))
   (A[l] is always symmetric)

   The regression only pays attention to datapoints that have
   their inputs weakly within the hypercube defined by xlo and xhi.

   The regression internally rescales all the data so that the computations
   should be well-conditioned even if xlo xhi brackets a small region far
   from the origin.

   ylo and yhi serve a slightly different role. They too are used
   for scaling the data to well-conditioned range.

    NOTE THE ORDR IN THE ARGUMENT LIST: xlo ylo xhi yhi. ROOM FOR CONFUSION

   regtype is as defined at top of file

   nu0 is the n0 prior noise strength parameter used in bayesreg.tex
   and Gelman et al.

   noise_frac_prior indirectly specifies the prior estimated noise
   value. The prior standard deviation of the noise is (roughly)

         noise_frac_prior * (yhi[l] - ylo[l])/2

   so if your third output column (l = 2) varied between
   5 and 8 and you expected the noise to be around 0.8, you might
   set   ylo[2] == 5, yhi[2] == 8, noise_frac = 8/15

   How unsure are you about the coefficients? Well first you should
   know that the coefficients will be with respect to input terms
   that have been scaled to the range [-1,+1].

     Good values are

                    2.0, /+ constant_coeff_prior_sdev +/
                    5.0, /+linear_coeff_prior_sdev +/
                    5.0 /+quadratic_coeff_prior_sdev +/
*/
/* rows may be NULL denotinf "use all datapoints" */
posters *mk_posters(dym *x,dym *y,ivec *rows,
                    dyv *xlo,dyv *ylo,
                    dyv *xhi,dyv *yhi,
                    int regtype,bprior *bp);

poster *posters_ref(posters *pos,int out_index);


/* Draw the expected (maximum aposteriaiai); values of the
   coefficients.
     po --- created with mk_posters
     out_index --- the output we are predicting (0 <= out_index < dym_cols(y););
                    where y was the matrix of outputs used in making the
                    posters.
		    *c --- get filled with the sample's c
                    *b -- a dyv is AM_MALLOCKED and stored in here
                    *a -- a dym is AM_MALLOCKED and stored in here.

  Our model y[out_index] = c + b^T x + x^t A ^x + noise
     has a posterior distribution on c,b,A. Here we give the
     most likely value from that distribution.

  NB... if mk_posters was created with regtype==LINEAR_REGDEG, then
  *a is filled with NULL.
*/
void make_posters_mean_cba(posters *po,int out_index,double *c,dyv **b,dym **a);

/* Draw a random sample quadratic form from the posteriors of a fit.
     po --- created with mk_posters
     out_index --- the output we are predicting (0 <= out_index < dym_cols(y););
                    where y was the matrix of outputs used in making the
                    posters.
		    *c --- get filled with the sample's c
                    *b -- a dyv is AM_MALLOCKED and stored in here
                    *a -- a dym is AM_MALLOCKED and stored in here.

  Our model y[out_index] = c + b^T x + x^t A ^x + noise
     has a posterior distribution on c,b,A. Here we draw from
     that distribution.

  NB... if mk_posters was created with regtype==LINEAR_REGDEG, then
  *a is filled with NULL.
*/
void make_posters_sample_cba(posters *po,int out_index,double *c,dyv **b,dym **a);


/* Returns the estimated noise. Looks at the posterior distribution
   on the noise, and returns three values.

    The low limit of the 95% confidence interval on the noise
    The median value
    The high limit.

   These values are all in original variables (y-values); coordinate system.

   NOTE THEY ARE STANDARD DEVIATIONS NOT VARIANCES.
*/
void posters_noise_sdev(posters *po,int out_index,double *lo95,double *median,double *hi95);

void free_reginfo(reginfo *ri);

reginfo *mk_copy_reginfo(reginfo *ri);

/* Free p and all its subcomponents */
void free_posters(posters *p);

void fprintf_reginfo(FILE *s,char *m1,reginfo *p,char *m2);

void fprintf_posters(FILE *s,char *m1,posters *p,char *m2);

dyv *reginfo_middle(reginfo *ri);
dyv *reginfo_range(reginfo *ri);
dyv *reginfo_ymiddle(reginfo *ri);
dyv *reginfo_yrange(reginfo *ri);

double reginfo_dsqd(reginfo *ri,dyv *x1,dyv *x2);

int reginfo_num_outputs(reginfo *ri);

dyv *posters_middle(posters *po);

int posters_num_outputs(posters *po);

int posters_regtype(posters *po);

void explain_posters(posters *po, int out_index);

void explain_posters_detailed(posters *po, int out_index);

/* rows may be NULL denotinf "use all datapoints" */
void make_sensible_vector_limits(dym *x,ivec *rows,dyv **xlo,dyv **xhi);

/* Same function as mk_reginfo() above except includes all
   datapoints...not only those inside a hypercube */
/* rows may be NULL denotinf "use all datapoints" */
reginfo *mk_reginfo_from_all_datapoints(dym *x,dym *y,ivec *rows,int regtype);

/* The prediction of the out_index'th output attribute has a Prob Density Function
   the is a t-distribution. This function computes and returns the
   dof, mean, variance of this t-distribution */
void posters_predict_marginal(posters *po,int out_index,dyv *in,
                              double *r_dof,double *r_mean,double *r_variance);

/* Predicts the value of the given input. Also gives the
   upper and lower XXX-%ile bounds on this prediction
   where XXX is supplied by the called (e.g. 0.975 for 95%-ile bounds);. */
double posters_predict(posters *po,int out_index,dyv *in,
                       double confidence_level, /* eg 0.975 for 95% bounds */
                       double *r_low_bound,double *r_high_bound);

/* You have terms matrix Z and
   outputs vector V.

     Z_ki = i'th term of kth datapoint
     V_k  = kth output

  Suppose you guess a set of regression coefficients beta, so your
  model is predicted output = beta^T terms of input

  Then the sum of squared residuals = (V - Z beta)^T (V - Z beta)
                                    = V^T V - 2 beta^T Z^T V + beta^T Z^T Z beta

  THIS IS TRUE NO MATTER HOW YOU GENERATED BETA. IT DOES NOT DEPEND ON
  BETA BEING THE LEAST SQUARES SOLUTION.
  */
double compute_sum_sqd_residuals(dym *ztz,dyv *ztv,double vtv,dyv *beta);

void diag_times_dyv(dyv *a,dyv *b,dyv *x);
dyv *mk_diag_times_dyv(dyv *a,dyv *b);

/* What is the probability we'd have discovered such a high-magnitude
   sample mean if the true mean was really 0.0? */
double t_significance_level(double dof,double mean,double variance);

/* What is the probability we'd have discovered such a high-magnitude
   fitted coefficient if the true coefficient was really 0.0? */
double posters_coefficient_significance_level(posters *po,
					      int coeff_index,int out_index);

/* What is the probability we'd have discovered such a high-magnitude
   fitted slope if the true slope was really 0.0? */
double posters_slope_significance_level(posters *po,
					int in_index,int out_index);

/* rows may be NULL meaning "use all rows" */
double posters_sum_sqd_residuals(posters *po,dym *x,dym *y,
				 ivec *rows,int out_index);

/* rows may be NULL meaning "use all rows" */
dyv *mk_dyv_from_subset_of_dym_rows(dym *x,ivec *rows,int col);

/* rows may be NULL meaning "use all rows" */
double posters_out_variance(dym *y,ivec *rows,int out_index);

/* rows may be NULL meaning "use all rows" */
double posters_fraction_variance_explained(posters *po,dym *x,dym *y,
					   ivec *rows,int out_index);

#endif /* #ifndef BAYESREG_H */

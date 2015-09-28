/* *
   File:        bayesreg.c
   Author:      Andrew W. Moore
   Created:     Wed Apr 23 15:43:43 EDT 1997
   Description: Bayesian Regression Analysis

   Copyright 1997, Carnegie Mellon University
*/

#include "bayesreg.h"
#include "./utils/stats.h"
/*#include "ongr.h"*/

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
char *regtype_to_name(int regtype)
{
  char *result = NULL;
  switch ( regtype )
  {
    case LINEAR_REGTYPE: result = "linear"; break;
    case CIRCLE_REGTYPE: result = "circle"; break;
    case ELLIPSE_REGTYPE: result = "ellipse"; break;
    case QUADRATIC_REGTYPE: result = "quadratic"; break;
    default: my_error("bad regtype");
  }
  return(result);
}

/* Return -1 if no such name */
int name_to_regtype(char *regtype_name)
{
  int result = -1;
  if ( eq_string(regtype_name,regtype_to_name(LINEAR_REGTYPE)) )
    result = LINEAR_REGTYPE;
  else if ( eq_string(regtype_name,regtype_to_name(CIRCLE_REGTYPE)) )
    result = CIRCLE_REGTYPE;
  else if ( eq_string(regtype_name,regtype_to_name(ELLIPSE_REGTYPE)) )
    result = ELLIPSE_REGTYPE;
  else if ( eq_string(regtype_name,regtype_to_name(QUADRATIC_REGTYPE)) )
    result = QUADRATIC_REGTYPE;
  return(result);
}

int num_regtypes()
{
  return NUM_REGTYPES;
}

int regtype_to_num_terms(int num_inputs,int regtype)
{
  int result = -1;
  switch ( regtype )
  {
    case LINEAR_REGTYPE: result = 1 + num_inputs; break;
    case CIRCLE_REGTYPE: result = 2 + num_inputs; break;
    case ELLIPSE_REGTYPE: result = 1 + 2 * num_inputs; break;
    case QUADRATIC_REGTYPE: result = (1+num_inputs) * (2+num_inputs) / 2;
                            break;
    default: my_error("bad regtype");
  }
  return(result);
}

int regtype_from_args(char *key,int argc,char *argv[],int default_regtype)
{
  char *s = string_from_args(key,argc,argv,regtype_to_name(default_regtype));
  int regtype = name_to_regtype(s);
  if ( regtype < 0 )
  {
    fprintf(stderr,"Unknown regtype: %s\n",s);
    my_error("regtype_from_args");
  }
  return(regtype);
}

/***** Main Bayesian Regression Code ****/

#define IFVERB(f) if (Verbosity >= (f))

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
                         double *r_sigma_n_sqd)
{
  int m = dym_cols(ztwz);
  dym *covw = mk_copy_dym(ztwz);
  int i;
  double s_squared;

  my_error("Not doing bayesian regression!");

  IFVERB(15) {
    printf("n = %g\n",n);
    fprintf_dym(stdout,"ztwz",ztwz,"\n");
    fprintf_dyv(stdout,"ztwy",ztwy,"\n");
    printf("ytwy = %g\n",ytwy);
  }

#ifndef AMFAST
  if ( dym_rows(ztwz) != m ||
       dyv_size(ztwy) != m ||
       dyv_size(beta_0) != m ||
       dyv_size(beta_0_variance) != m )
    my_error("Bad matrix/vector shape in call to bayesian regression");
#endif

  for ( i = 0 ; i < m ; i++ )
    dym_increment(covw,i,i,1.0 / dyv_ref(beta_0_variance,i));

  IFVERB(15) fprintf_dym(stdout,"covw",covw,"\n");

  *r_vbeta = mk_invert_spd_cholesky(covw);

  if ( *r_vbeta == NULL )
  {
    fprintf_dym(stdout,"covw",covw,"\n");
    my_error("Covariance plus prior should be symmetric positive definite");
  }

  IFVERB(15) fprintf_dym(stdout,"*r_vbeta",*r_vbeta,"\n");

  *r_beta_hat = mk_dyv(m);
  for ( i = 0 ; i < m ; i++ )
    dyv_set(*r_beta_hat,i, dyv_ref(beta_0,i) / dyv_ref(beta_0_variance,i) +
                           dyv_ref(ztwy,i));
  dym_times_dyv(*r_vbeta,*r_beta_hat,*r_beta_hat);

  IFVERB(15) fprintf_dyv(stdout,"*r_beta_hat",*r_beta_hat,"\n");

  s_squared = ytwy - 2 * dyv_scalar_product(*r_beta_hat,ztwy);

  for ( i = 0 ; i < m ; i++ )
  {
    double beta_hat_i = dyv_ref(*r_beta_hat,i);
    int j;

    s_squared += beta_hat_i * beta_hat_i * dym_ref(ztwz,i,i);

    for ( j = 0 ; j < i ; j++ )
      s_squared += 2 * beta_hat_i * dyv_ref(*r_beta_hat,j) * dym_ref(ztwz,i,j);
  }

  for ( i = 0 ; i < m ; i++ )
  {
    double dbi = dyv_ref(beta_0,i) - dyv_ref(*r_beta_hat,i);
    s_squared += dbi * dbi / dyv_ref(beta_0_variance,i);
  }

  s_squared /= real_max(1e-10,n);

  *r_nu = n0 + n;
  *r_sigma_n_sqd = (n0 * sigma_squared0 + n * s_squared) / (n0 + n);
  free_dym(covw);
}

brdat *mk_default_brdat(int regtype)
{
  brdat *br = AM_MALLOC(brdat);
  int num_inputs = 1;
  int num_terms = regtype_to_num_terms(num_inputs,regtype);
  br->x = mk_dyv(0);
  br->w = mk_dyv(0);
  br->y = mk_dyv(0);
  br->cx = 0.0;
  br->rx = 1.0;
  br->n0 = 1.0;
  br->ss0 = 1.0;
  br->regtype = regtype;
  br->b0 = mk_zero_dyv(num_terms);
  br->b0v = mk_constant_dyv(num_terms,100.0);
  br->xlo = 0.0;
  br->xhi = 1.0;
  return(br);
}

void free_brdat(brdat *br)
{
  free_dyv(br->x);
  free_dyv(br->w);
  free_dyv(br->y);
  free_dyv(br->b0);
  free_dyv(br->b0v);
  AM_FREE(br,brdat);
}

bool update_double_from_string_array(double *x,char *name,string_array *sa,
                                     char **r_errmess)
{
  bool result = string_array_size(sa) == 3 &&
                eq_string(string_array_ref(sa,0),name);
  if ( result )
  {
    if ( is_a_number(string_array_ref(sa,2)) )
      *x = atof(string_array_ref(sa,2));
    else
    {
      char buff[1000];
      sprintf(buff,"%s is not a number",string_array_ref(sa,2));
      *r_errmess = mk_copy_string(buff);
    }
  }
  return(result);
}

bool update_dyv_from_string_array(dyv *x,char *name,string_array *sa,
                                  char **r_errmess)
{
  bool result = string_array_size(sa) == 3;
  *r_errmess = FALSE;
  if ( result )
  {
    int i;
    result = FALSE;
    for ( i = 0 ; !result && i < dyv_size(x) ; i++ )
    {
      char buff[1000];
      double z = dyv_ref(x,i);
      sprintf(buff,"%s[%d]",name,i);
      if ( update_double_from_string_array(&z,buff,sa,r_errmess) )
      {
        result = TRUE;
        if ( *r_errmess == NULL )
          dyv_set(x,i,z);
      }
    }
  }
  return(result);
}

void update_brdat_from_string_array(brdat *br,string_array *sa,
                                    char **r_errmess)
{
  *r_errmess = NULL;
  if ( string_array_size(sa) >= 3 &&
       eq_string(string_array_ref(sa,0),"add") )
  {
    char *s1 = string_array_ref(sa,1);
    char *s2 = string_array_ref(sa,2);
    if ( is_a_number(s1) && is_a_number(s2) )
    {
      double w = ( string_array_size(sa) > 3 &&
                   is_a_number(string_array_ref(sa,3)) ) ?
                 atof(string_array_ref(sa,3)) : 1.0;
      add_to_dyv(br->x,atof(s1));
      add_to_dyv(br->y,atof(s2));
      add_to_dyv(br->w,w);
    }
  }
  else if ( !(update_dyv_from_string_array(br->x,"x",sa,r_errmess) ||
         update_dyv_from_string_array(br->y,"y",sa,r_errmess) ||
         update_dyv_from_string_array(br->w,"w",sa,r_errmess) ||
         update_dyv_from_string_array(br->b0,"b0",sa,r_errmess) ||
         update_dyv_from_string_array(br->b0v,"b0v",sa,r_errmess) ||
         update_double_from_string_array(&br->cx,"cx",sa,r_errmess) ||
         update_double_from_string_array(&br->rx,"rx",sa,r_errmess) ||
         update_double_from_string_array(&br->n0,"n0",sa,r_errmess) ||
         update_double_from_string_array(&br->ss0,"ss0",sa,r_errmess) ||
         update_double_from_string_array(&br->xlo,"xlo",sa,r_errmess) ||
         update_double_from_string_array(&br->xhi,"xhi",sa,r_errmess)) )
  {
    if ( *r_errmess == NULL )
      *r_errmess = mk_copy_string("Unrecognized variable");
  }
}

dyv *mk_zdyv_from_x(double x,double cx,double rx,int regtype)
{
  int num_inputs = 1;
  int num_terms = regtype_to_num_terms(num_inputs,regtype);
  dyv *z = mk_dyv(num_terms);
  double u = 1.0;
  double zval = (x - cx) / rx;
  int i;
  for ( i = 0 ; i < num_terms ; i++ )
  {
    dyv_set(z,i,u);
    u *= zval;
  }
  return(z);
}

void make_regress(brdat *br,dym **r_ztwz,dyv **r_ztwy,double *r_ytwy)
{
  int num_inputs = 1;
  int m = regtype_to_num_terms(num_inputs,br->regtype);
  int k;
  int i,j;

  *r_ztwz = mk_zero_dym(m,m);
  *r_ztwy = mk_zero_dyv(m);
  *r_ytwy = 0.0;

  for ( k = 0 ; k < dyv_size(br->x) ; k++ )
  {
    dyv *zk = mk_zdyv_from_x(dyv_ref(br->x,k),br->cx,br->rx,br->regtype);
    double wk = dyv_ref(br->w,k);
    dyv *wzk = mk_dyv_scalar_mult(zk,wk);
    double yk = dyv_ref(br->y,k);

    for ( i = 0 ; i < m ; i++ )
      for ( j = 0 ; j <= i ; j++ )
        dym_increment(*r_ztwz,i,j,dyv_ref(zk,i) * dyv_ref(wzk,j));

    for ( i = 0 ; i < m ; i++ )
      dyv_increment(*r_ztwy,i,dyv_ref(wzk,i) * yk);

    *r_ytwy += wk * yk * yk;
    free_dyv(zk);
    free_dyv(wzk);
  }

  for ( i = 0 ; i < m ; i++ )
    for ( j = i+1 ; j < m ; j++ )
      dym_set(*r_ztwz,i,j,dym_ref(*r_ztwz,j,i));
}

double y_mean_cdf_inv(brdat *br,double x,double prob,
                      dyv *beta_hat,dym *vbeta,
                      double nu,double sigma_n_sqd)
{
  dyv *z = mk_zdyv_from_x(x,br->cx,br->rx,br->regtype);
  double y_cdf_inv = t_cdf_inv(prob,
                               dyv_scalar_product(z,beta_hat),
                               sqrt(sigma_n_sqd * dym_xt_a_x_value(z,vbeta)),
                               nu);
  free_dyv(z);
  return(y_cdf_inv);
}


void print_mean_confidence_interval(brdat *br,
                                    double x,dyv *beta_hat,dym *vbeta,
                                    double nu,double sigma_n_sqd)
{
  dyv *ps = mk_dyv_5(0.01,0.025,0.5,0.975,0.99);
  int i;
  for ( i = 0 ; i < dyv_size(ps) ; i++ )
  {
    double prob = dyv_ref(ps,i);
    double y_cdf_inv = y_mean_cdf_inv(br,x,prob,beta_hat,vbeta,nu,sigma_n_sqd);
    printf("Prob(y < %9g | x = %9g) = %9g\n",y_cdf_inv,x,prob);
  }
}

double sample_from_gamma(double alpha,double beta)
{
  double prob = range_random(0.0,1.0);
  double result = gamma_cdf_inv(prob,alpha,beta);
  return(result);
}

double sample_from_chi_squared(double nu)
{
  double alpha = nu / 2.0;
  double beta = 0.5;
  double result = (nu > 1e4) ? nu : sample_from_gamma(alpha,beta);
  return(result);
}

double sample_from_scaled_inverse_chi_squared(double nu,double s_squared)
{
  double x = sample_from_chi_squared(nu);
  double result = nu * s_squared / x;
  return(result);
}

dyv *mk_vector_of_unit_normal_samples(int size)
{
  dyv *result = mk_dyv(size);
  int i;
  for ( i = 0 ; i < size ; i++ )
    dyv_set(result,i,gen_gauss());
  return(result);
}

dyv *mk_sample_from_normal(dyv *mu,dym *covariance)
{
  dym *lower_triangular = mk_dym(dyv_size(mu),dyv_size(mu));
  dyv *result = mk_vector_of_unit_normal_samples(dyv_size(mu));
  bool ok = attempt_cholesky_decomp(covariance,lower_triangular);
  if ( !ok )
  {
    void enforce_dym_symmetry(dym *d);
    dym *cov_copy = mk_copy_dym(covariance);
    int i;
    for (i = 0 ; i < dym_rows(covariance) ; i++ )
      dym_increment(cov_copy,i,i,0.01 * (1+fabs(dym_ref(covariance,i,i))));
    enforce_dym_symmetry(cov_copy);

    ok = attempt_cholesky_decomp(cov_copy,lower_triangular);

    if ( !ok )
    {
      fprintf_dym(stdout,"covariance",covariance,"\n");
      fprintf_dym(stdout,"cov_copy",cov_copy,"\n");
      my_error("mk_sample_from_normal. Matrix not SPD");
    }
    free_dym(cov_copy);
  }
  dym_times_dyv(lower_triangular,result,result);
  dyv_plus(result,mu,result);
  free_dym(lower_triangular);
  return(result);
}

dyv *mk_sample_beta(dyv *beta_hat,dym *vbeta,double nu,double sigma_n_sqd)
{
  double sigma_squared = sample_from_scaled_inverse_chi_squared(nu,sigma_n_sqd);
  dyv *zero = mk_zero_dyv(dyv_size(beta_hat));
  dyv *result = mk_sample_from_normal(zero,vbeta);
  dyv_scalar_mult(result,sqrt(sigma_squared),result);
  dyv_plus(result,beta_hat,result);
  free_dyv(zero);
  IFVERB(15) printf("Sampled beta0 = %g\n",dyv_ref(result,0));
  return(result);
}

void print_brdat(brdat *br)
{
  int i,k;
  int num_inputs = 1;
  int num_terms = regtype_to_num_terms(num_inputs,br->regtype);

  printf("%9s %9s %9s %9s ","k","x","y","w");
  for ( i = 0 ; i < num_terms ; i++ )
    printf("%8s%d","z",i);
  printf("\n");
  for ( k = 0 ; k < dyv_size(br->x) ; k++ )
  {
    dyv *zk = mk_zdyv_from_x(dyv_ref(br->x,k),br->cx,br->rx,br->regtype);
    printf("%9g %9g %9g %9g ",(double)k,dyv_ref(br->x,k),dyv_ref(br->y,k),
           dyv_ref(br->w,k));
    for ( i = 0 ; i < dyv_size(zk) ; i++ )
      printf("%9g ",dyv_ref(zk,i));
    printf("\n");
    free_dyv(zk);
  }

  printf("\n");

  printf("%9s %9s %9s %9s ","b0","","","");
  for ( i = 0 ; i < dyv_size(br->b0) ; i++ )
    printf("%9g ",dyv_ref(br->b0,i));
  printf("\n");

  printf("%9s %9s %9s %9s ","b0v","","","");
  for ( i = 0 ; i < dyv_size(br->b0v) ; i++ )
    printf("%9g ",dyv_ref(br->b0v,i));
  printf("\n");

  printf("cx = %g, rx = %g\n",br->cx,br->rx);
}

void print_regress(dyv *beta_hat,
                   dym *vbeta,
                   double nu,double sigma_n_sqd)
{
  int i;
  printf("nu = %g, sigma_n_sqd = %g\n",nu,sigma_n_sqd);
  printf("%9s %9s %9s %9s ","betahat","","","");
  for ( i = 0 ; i < dyv_size(beta_hat) ; i++ )
    printf("%9g ",dyv_ref(beta_hat,i));
  printf("\n");

  for ( i = 0 ; i < dyv_size(beta_hat) ; i++ )
  {
    int j;
    printf("%9s %9s %9s %9s ",(i==0)?"vbeta":"","","","");
    for ( j = 0 ; j < dyv_size(beta_hat) ; j++ )
      printf("%9g ",dym_ref(vbeta,i,j));
    printf("\n");
  }
}

/* Marginal on specified coefficient is a t-distribution
   with mean *r_mu, variance *r_sigma_squared and dof *r_nu */
void coefficient_marginal(dyv *beta_hat,dym *vbeta,
                          double nu,double sigma_n_sqd,
                          int coeff_index,
                          double *r_nu,double *r_mu,double *r_sigma_squared)
{
  *r_nu = nu;
  *r_mu = dyv_ref(beta_hat,coeff_index);
  *r_sigma_squared = dym_ref(vbeta,coeff_index,coeff_index) * sigma_n_sqd;
}

/* Marginal on prediction (at scaled polynomial point z) is a t-distribution
   with mean *r_mu, variance *r_sigma_squared and dof *r_nu */
/*              mean = poi -> beta_hat . z
            variance = z^T poi->vbeta z sigma_n_sqd
                 dof = nu
*/
void predict_marginal(dyv *beta_hat,dym *vbeta,double nu,double sigma_n_sqd,
                      dyv *z,
                      double *r_nu,double *r_mu,double *r_sigma_squared)
{
  *r_nu = nu;
  *r_mu = dyv_scalar_product(beta_hat,z);
  *r_sigma_squared = sigma_n_sqd * dym_xt_a_x_value(z,vbeta);
}

double t_half_95ci_width(double nu,double mu,double sigma_squared)
{
  return(t_cdf_inv(0.975,mu,sqrt(sigma_squared),nu) - mu);
}

double scaled_inverse_chi_squared_cdf_inv(double prob,double nu,
                                          double s_squared)
{
  double from_gamma = gamma_cdf_inv(1-prob,nu/2,0.5);
  double from_sics = nu * s_squared / from_gamma;
  return(from_sics);
}

void print_marginals(brdat *br,dyv *beta_hat,dym *vbeta,
                     double nu,double sigma_n_sqd)
{
  int i;
  int num_inputs = 1;
  int num_terms = regtype_to_num_terms(num_inputs,br->regtype);

  for ( i = 0 ; i < num_terms ; i++ )
  {
    double mnu,mmu,msigsqd,h;

    coefficient_marginal(beta_hat,vbeta,nu,sigma_n_sqd,i,&mnu,&mmu,&msigsqd);
    h = t_half_95ci_width(mnu,mmu,msigsqd);
    printf("b%d ~ t(dof=%g,mean=%g,sigsqd=%g) 95%%CI : [%g,%g]\n",
           i,mnu,mmu,msigsqd,mmu-h,mmu+h);
  }
  printf("sigma_sqd ~ Inv-Chi^2(dof=%g,s^2=%g).\n",nu,sigma_n_sqd);
  printf("sigma_sqd mean ");

  if ( nu <= 2.0 )
    printf("undefined ");
  else
    printf("%g ",nu * sigma_n_sqd / (nu - 2));

  printf("mode=%g 95%%CI : [%g,%g]\n",nu*sigma_n_sqd/(nu+2),
         scaled_inverse_chi_squared_cdf_inv(0.025,nu,sigma_n_sqd),
         scaled_inverse_chi_squared_cdf_inv(0.975,nu,sigma_n_sqd));
}

int dyv_num_nonzero(dyv *x)
{
  int result = 0;
  int i;
  for ( i = 0 ; i < dyv_size(x) ; i++ )
  {
    double v = dyv_ref(x,i);
    if ( v != 0.0 ) result += 1;
  }
  return(result);
}


/* In the following , i,j index over elements of input points (x's)
                   they also index over z vectors (vectors of basis function evaluations)

              k indexes over datapoints

              l indexes over outputs.


   Define v[l] as a version of y[l] scaled thusly:   v[l] = ( y[l] - ymid[l] ) / yrange[l]
                                               so    y[l] = ymid[l] + yrange[l] * v[l]

   A z vector corresponding to an x input in a quadratic or linear regression is
   defined shortly (assume x has components x[0] ... x[M-1])

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
     z[M + M(M+1)/2 - 1] = u[M-2] * u[M-1]
     z[M + M(M+1)/2] = u[M-1] * u[M-1]

   where the quadratic terms exist if and only if regtype==LINEAR_REGTYPE.

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

dyv *mk_u_from_x(reginfo *ri,dyv *x)
{
  int usin_size = dyv_size(x);
  dyv *result = mk_dyv_subtract(x,ri->mid);
  int i;
  for ( i = 0 ; i < usin_size ; i++ )
    dyv_set(result,i,dyv_ref(result,i) / dyv_ref(ri->range,i));
  return(result);
}


dyv *mk_z_from_u_and_regtype(dyv *u,int regtype)
{
  int usin_size = dyv_size(u);
  int z_size = regtype_to_num_terms(usin_size,regtype);
  int i,j;
  dyv *z = mk_dyv(z_size);
  dyv_set(z,0,1.0);
  for ( i = 0 ; i < usin_size ; i++ )
    dyv_set(z,i+1,dyv_ref(u,i));

  switch ( regtype )
  {
    case LINEAR_REGTYPE:
    {
      /* Skip */
    }
    break;
    case CIRCLE_REGTYPE:
    {
      double sum_uu = 0.0;
      for ( i = 0 ; i < usin_size ; i++ )
      {
        double ui = dyv_ref(u,i);
        sum_uu += ui * ui;
      }
      dyv_set(z,usin_size+1,sum_uu);
    }
    break;
    case ELLIPSE_REGTYPE:
    {
      for ( i = 0 ; i < usin_size ; i++ )
      {
        double ui = dyv_ref(u,i);
        dyv_set(z,usin_size+1+i,ui * ui);
      }
    }
    break;
    case QUADRATIC_REGTYPE:
    {
      int z_ptr = usin_size + 1;

      for ( i = 0 ; i < usin_size ; i++ )
        for ( j = i ; j < usin_size ; j++ )
        {
          dyv_set(z,z_ptr,dyv_ref(u,i) * dyv_ref(u,j));
          z_ptr += 1;
        }
      if ( z_ptr != z_size ) my_error("assertion");
    }
    break;
    default: my_error("bad regtype");
  }

  return(z);
}

dyv *mk_z_from_u(reginfo *ri,dyv *u)
{
  dyv *z = mk_z_from_u_and_regtype(u,ri->regtype);
  return(z);
}

dyv *mk_z_from_usin(reginfo *ri,dyv *x)
{
  dyv *u = mk_u_from_x(ri,x);
  dyv *z = mk_z_from_u(ri,u);
  free_dyv(u);
  return(z);
}

double v_from_y(reginfo *ri,double y,int out_index)
{
  double ymid = dyv_ref(ri->ymid,out_index);
  double yrange = dyv_ref(ri->yrange,out_index);
  return((y - ymid)/yrange);
}

double y_from_v(reginfo *ri,double v,int out_index)
{
  double ymid = dyv_ref(ri->ymid,out_index);
  double yrange = dyv_ref(ri->yrange,out_index);
  return(ymid + v * yrange);
}

dyv *mk_v_from_usout(reginfo *ri,dyv *y)
{
  dyv *v = mk_dyv(dyv_size(y));
  int i;
  for ( i = 0 ; i < dyv_size(v) ; i++ )
    dyv_set(v,i,v_from_y(ri,dyv_ref(y,i),i));
  return(v);
}

void print_cba(char *cname,char *bname,char *aname,
               double c,dyv *b,dym *a)
{
  printf("Quadratic form %s + %s^T x + 0.5 x^T %s x where\n",cname,bname,aname);
  printf("%s = %g\n",cname,c);
  fprintf_dyv(stdout,bname,b,"\n");
  fprintf_dym(stdout,aname,a,"\n");
  wait_for_key();
}

int reginfo_num_inputs(reginfo *ri)
{
  return(dyv_size(ri->mid));
}

int reginfo_num_terms(reginfo *ri)
{
  return(regtype_to_num_terms(reginfo_num_inputs(ri),reginfo_regtype(ri)));
}

/*** bpriors are used to make posters ***/

bprior *mk_empty_bprior()
{
  bprior *bp = AM_MALLOC(bprior);
  bp -> non_informative = FALSE;
  bp -> known_sdev = FALSE;
  bp -> noise_frac_prior = -777.777;
  bp -> nu0 = -777.777;
  bp -> constant_coeff_prior_sdev = -777.777;
  bp -> linear_coeff_prior_sdev = -777.777;
  bp -> quadratic_coeff_prior_sdev = -777.777;
  return bp;
}

bprior *mk_non_informative_bprior()
{
  bprior *bp = mk_empty_bprior();
  bp->non_informative = TRUE;
  bp->known_sdev = FALSE;
  return(bp);
}

bprior *mk_known_sdev_non_informative_coeffs_prior(double noise_frac_prior)
{
  bprior *bp = mk_empty_bprior();
  bp->non_informative = FALSE;
  bp->known_sdev = TRUE;
  bp->noise_frac_prior = noise_frac_prior;
  return(bp);
}

bprior *mk_informative_bprior(double nu0,double noise_frac_prior,
                              double constant_coeff_prior_sdev,
                              double linear_coeff_prior_sdev,
                              double quadratic_coeff_prior_sdev)
{
  bprior *bp = mk_empty_bprior();
  bp->non_informative = FALSE;
  bp->known_sdev = FALSE;
  bp->nu0 = nu0;
  bp->noise_frac_prior = noise_frac_prior;
  bp->constant_coeff_prior_sdev = constant_coeff_prior_sdev;
  bp->linear_coeff_prior_sdev = linear_coeff_prior_sdev;
  bp->quadratic_coeff_prior_sdev = quadratic_coeff_prior_sdev;
  return(bp);
}

void free_bprior(bprior *bp)
{
  AM_FREE(bp,bprior);
}

void fprintf_bprior(FILE *s,char *m1,bprior *bp,char *m2)
{
  fprintf(s,"%s -> non_informative = %d%s",m1,bp->non_informative,m2);
  fprintf(s,"%s -> known_sdev = %d%s",m1,bp->known_sdev,m2);
  fprintf(s,"%s -> noise_frac_prior = %g%s",m1,bp->noise_frac_prior,m2);
  fprintf(s,"%s -> nu0 = %g%s",m1,bp->nu0,m2);
  fprintf(s,"%s -> constant_coeff_prior_sdev = %g%s",m1,bp->constant_coeff_prior_sdev,m2);
  fprintf(s,"%s -> linear_coeff_prior_sdev = %g%s",m1,bp->linear_coeff_prior_sdev,m2);
  fprintf(s,"%s -> quadratic_coeff_prior_sdev = %g%s",m1,bp->quadratic_coeff_prior_sdev,m2);
}

bprior *mk_copy_bprior(bprior *bp)
{
  bprior *newbp = AM_MALLOC(bprior);
  newbp -> non_informative = bp -> non_informative;
  newbp -> known_sdev = bp -> known_sdev;
  newbp -> noise_frac_prior = bp -> noise_frac_prior;
  newbp -> nu0 = bp -> nu0;
  newbp -> constant_coeff_prior_sdev = bp -> constant_coeff_prior_sdev;
  newbp -> linear_coeff_prior_sdev = bp -> linear_coeff_prior_sdev;
  newbp -> quadratic_coeff_prior_sdev = bp -> quadratic_coeff_prior_sdev;
  return(newbp);
}

reginfo *posters_reginfo(posters *po)
{
  return(po->ri);
}

int reginfo_regtype(reginfo *ri)
{
  return(ri->regtype);
}

int posters_num_inputs(posters *po)
{
  return(reginfo_num_inputs(posters_reginfo(po)));
}

void make_cpbpap_from_beta(reginfo *ri,dyv *beta,
                           double *cp,dyv **bp,dym **ap)
{
  int i,j;
  int usin_size = reginfo_num_inputs(ri);
  int regtype = reginfo_regtype(ri);

  *cp = dyv_ref(beta,0);

  *bp = mk_dyv(usin_size);
  for ( i = 0 ; i < usin_size ; i++ )
    dyv_set(*bp,i,dyv_ref(beta,i+1));

  if ( regtype == LINEAR_REGTYPE )
    *ap = NULL;
  else
  {
    int z_ptr = usin_size + 1;
    *ap = mk_zero_dym(usin_size,usin_size);

    switch ( regtype )
    {
      case CIRCLE_REGTYPE:
      {
        double a = dyv_ref(beta,z_ptr);
        for ( i = 0 ; i < usin_size ; i++ )
          dym_set(*ap,i,i,2 * a);
      }
      break;
      case ELLIPSE_REGTYPE:
      {
        for ( i = 0 ; i < usin_size ; i++ )
          dym_set(*ap,i,i,2 * dyv_ref(beta,z_ptr+i));
      }
      break;
      case QUADRATIC_REGTYPE:
      {
        for ( i = 0 ; i < usin_size ; i++ )
          for ( j = i ; j < usin_size ; j++ )
          {
            double aij = dyv_ref(beta,z_ptr);
            if ( i == j )
              dym_set(*ap,i,j,2 * aij);
            else
            {
              dym_set(*ap,i,j,aij);
              dym_set(*ap,j,i,aij);
            }
            z_ptr += 1;
          }
      }
      break;
      default: my_error("Bad regtype");
    }
  }
#ifdef DEBUG
  print_cba("cp","bp","ap",*cp,*bp,*ap);
#endif
}

dyv *mk_r_from_reginfo(reginfo *ri)
{
  dyv *r = mk_dyv(reginfo_num_inputs(ri));
  int i;
  for ( i = 0 ; i < dyv_size(r) ; i++ )
    dyv_set(r,i,1.0 / dyv_ref(ri->range,i));
  return(r);
}

void make_cba_from_cpbpap(reginfo *ri,double cp,dyv *bp,dym *ap,int out_index,
                          double *c,dyv **b,dym **a)
{
  int usin_size = dyv_size(bp);
  dyv *r = mk_r_from_reginfo(ri);
  dyv *r_m = mk_diag_times_dyv(r,ri->mid);
  dyv *ap_r_m = (ap==NULL) ? mk_zero_dyv(dyv_size(r_m)) :
                             mk_dym_times_dyv(ap,r_m);
  double ym = dyv_ref(ri->ymid,out_index);
  double yr = dyv_ref(ri->yrange,out_index);

  *c = ym + yr * ( cp -
                   dyv_scalar_product(bp,r_m) +
                   0.5 * dyv_scalar_product(r_m,ap_r_m) );

  *b = mk_dyv_subtract(bp,ap_r_m);
  diag_times_dyv(r,*b,*b);
  dyv_scalar_mult(*b,yr,*b);

  if ( ri->regtype == LINEAR_REGTYPE )
    *a = NULL;
  else
  {
    int i,j;
    *a = mk_dym(usin_size,usin_size);
    for ( i = 0 ; i < usin_size ; i++ )
      for ( j = 0 ; j < usin_size ; j++ )
        dym_set(*a,i,j,yr * dyv_ref(r,i) * dyv_ref(r,j) * dym_ref(ap,i,j));
  }
  free_dyv(r);
  free_dyv(r_m);
  free_dyv(ap_r_m);
}

void make_cba_from_beta(reginfo *ri,dyv *beta,int out_index,
                        double *c,dyv **b,dym **a)
{
  double cp;
  dyv *bp;
  dym *ap;
  IFVERB(15) fprintf_dyv(stdout,"beta",beta,"\n");
  make_cpbpap_from_beta(ri,beta,&cp,&bp,&ap);
  make_cba_from_cpbpap(ri,cp,bp,ap,out_index,
                       c,b,a);
  free_dyv(bp);
  if ( ap != NULL ) free_dym(ap);
#ifdef DEBUG
  print_cba("c","b","a",*c,*b,*a);
#endif
}

bool dyv_is_outside_region(dyv *x,dyv *lo,dyv *hi)
{
  return( !(dyv_weakly_dominates(hi,x) &&
            dyv_weakly_dominates(x,lo)) );
}


/* rows may be NULL denotinf "use all datapoints" */
/* Uses only the points that are in "rows" and that are
   not outside the hyper-rectangle
   defined by xlo and xhi.
   xlo and xhi may be NULL denoting "infinitely large hyperrect" */
dym *mk_ztz_from_usins(reginfo *ri,dym *usins,ivec *rows,dyv *xlo,dyv *xhi)
{
  int num_points = (rows==NULL)?dym_rows(usins) : ivec_size(rows);
  int usin_size = dym_cols(usins);
  int z_size = regtype_to_num_terms(usin_size,reginfo_regtype(ri));
  dym *ztz = mk_zero_dym(z_size,z_size);
  int n;

  for ( n = 0 ; n < num_points ; n++ )
  {
    int row = (rows==NULL) ? n : ivec_ref(rows,n);
    dyv *usin = mk_dyv_from_dym_row(usins,row);
    if ( xlo==NULL || !dyv_is_outside_region(usin,xlo,xhi) )
    {
      dyv *z_dyv = mk_z_from_usin(ri,usin);
      int i,j;

      for ( i = 0 ; i < z_size ; i++ )
      {
        double zki = dyv_ref(z_dyv,i);
        for ( j = 0 ; j < z_size ; j++ )
          dym_increment(ztz,i,j,zki * dyv_ref(z_dyv,j));
      }
      free_dyv(z_dyv);
    }
    free_dyv(usin);
  }
  return(ztz);
}

/* rows may be NULL denotinf "use all datapoints" */
/* Uses only the points that are in "rows" and that are
   not outside the hyper-rectangle
   defined by xlo and xhi.
   xlo and xhi may be NULL denoting "infinitely large hyperrect" */
void make_reg_matrices_from_usins_and_usouts(reginfo *ri,dym *usins,
                                             dym *usouts,int out_index,ivec *rows,
                                             dyv *xlo,dyv *xhi,
                                             dyv **ztv,double *vtv)
{
  int num_points = (rows==NULL) ? dym_rows(usins) : ivec_size(rows);
  int usin_size = dym_cols(usins);
  int z_size = regtype_to_num_terms(usin_size,reginfo_regtype(ri));
  int j;

  *ztv = mk_zero_dyv(z_size);
  *vtv = 0.0;

  for ( j = 0 ; j < num_points ; j++ )
  {
    int row = (rows==NULL) ? j : ivec_ref(rows,j);
    dyv *usin = mk_dyv_from_dym_row(usins,row);
    if ( xlo == NULL || !dyv_is_outside_region(usin,xlo,xhi) )
    {
      dyv *z_dyv = mk_z_from_usin(ri,usin);
      double vk = v_from_y(ri,dym_ref(usouts,row,out_index),out_index);
      int i;

      *vtv += vk * vk;

      for ( i = 0 ; i < z_size ; i++ )
      {
        double zki = dyv_ref(z_dyv,i);
        dyv_increment(*ztv,i,vk * zki);
      }
      free_dyv(z_dyv);
    }
    free_dyv(usin);
  }
}

poster *mk_informative_poster(reginfo *ri,
                    dym *ztz,dyv *ztv,double vtv,
                    double n0,double noise_frac_prior,
                    double constant_coeff_prior_sdev,
                    double linear_coeff_prior_sdev,
                    double quadratic_coeff_prior_sdev)
{
  poster *p = AM_MALLOC(poster);
  double sigsqd0 = real_square(noise_frac_prior);
  double sigsqd_beta_constant =
    real_square(constant_coeff_prior_sdev) / sigsqd0;
  double sigsqd_beta_linear =
    real_square(linear_coeff_prior_sdev) / sigsqd0;
  double sigsqd_beta_quadratic =
    real_square(quadratic_coeff_prior_sdev) / sigsqd0;
  int usin_size = reginfo_num_inputs(ri);
  int regtype = reginfo_regtype(ri);
  int z_size = dym_rows(ztz);
  dyv *beta_0 = mk_zero_dyv(z_size);
  dyv *beta_0_variance = mk_dyv(z_size);
  int i;

  dyv_set(beta_0_variance,0,sigsqd_beta_constant);
  for ( i = 1 ; i <= usin_size ; i++ )
    dyv_set(beta_0_variance,i,sigsqd_beta_linear);
  if ( regtype != LINEAR_REGTYPE )
    for ( ; i < z_size ; i++ )
      dyv_set(beta_0_variance,i,sigsqd_beta_quadratic);

  IFVERB(15) {
    fprintf_double(stdout,"n",dym_ref(ztz,0,0),"\n");
    fprintf_dyv(stdout,"beta_0",beta_0,"\n");
    fprintf_dyv(stdout,"beta_0_variance",beta_0_variance,"\n");
    fprintf_double(stdout,"n0",n0,"\n");
    fprintf_double(stdout,"sigsqd0",sigsqd0,"\n");
  }

  bayesian_regression(dym_ref(ztz,0,0),ztz,ztv,vtv,
                      beta_0,beta_0_variance,
                      n0,sigsqd0,
                      &p->beta_hat,&p->vbeta,&p->nu,&p->sigma_n_sqd);

  free_dyv(beta_0);
  free_dyv(beta_0_variance);

  return(p);
}

/* Robust inversion. tries to invert a (which must be
   symmetric positive definite, but may be singular).

   If it succeeds returns a^-1 and sets *r_diag = 0.

   If fails, returns (a + D)^-1 where D = Diag(*r_diag,*r_diag, ... *r_diag)
     where *r_diag is a small positive real number designed to make
     a + D non-singular.
*/
dym *mk_robust_invert_spd(dym *a,double *r_diag)
{
  dym *a_inv = mk_invert_spd_cholesky(a);
  double delta = 1e-4;

  *r_diag = 0.0;

  while ( a_inv == NULL && delta < 1e4 )
  {
    void enforce_dym_symmetry(dym *a);
    dym *a_copy = mk_copy_dym(a);
    int i;
    enforce_dym_symmetry(a_copy);
    for ( i = 0 ; i < dym_rows(a_copy) ; i++ )
      dym_increment(a_copy,i,i,delta);
    a_inv = mk_invert_spd_cholesky(a_copy);
    *r_diag = delta;
    delta *= 2.0;
    free_dym(a_copy);
  }

  if ( a_inv == NULL )
  {
    fprintf_dym(stderr,"a",a,"\n");
    my_error("even with extreme measures (adding 10000 to main diagonal) "
             "I could not invert this matrix with cholesky decomposition");
  }

  return(a_inv);
}

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
double compute_sum_sqd_residuals(dym *ztz,dyv *ztv,double vtv,dyv *beta)
{
  double aaa = dyv_scalar_product(beta,ztv);
  double bbb = dym_xt_a_x_value(beta,ztz);
  double result = vtv - 2.0 * aaa + bbb;
  double magnitude = real_max(real_max(fabs(vtv),fabs(aaa)),fabs(bbb));

  if ( result < 0.0 )
  {
    magnitude = real_max(magnitude,dyv_max(beta));
    magnitude = real_max(magnitude,-dyv_min(beta));
  }

  if ( result/real_max(1e-5,magnitude) < -2e-3 )
  {
    fprintf_dym(stderr,"ztz",ztz,"\n");
    fprintf_dyv(stderr,"ztv",ztv,"\n");
    fprintf_dyv(stderr,"beta",beta,"\n");
    printf("vtv = %g, sum_squared_residuals = %g\n",vtv,result);
    my_error("-ve sum squared residuals. This can't happen. It is NOT any fault of the calculation of beta. It should NOT be caused by ztz being singular.");
  }
  else if ( result < 0.0 )
    result = 0.0;

  return result;
}

poster *mk_non_informative_poster(dym *ztz,dyv *ztv,double vtv)
{
  double n = dym_ref(ztz,0,0);
  double k = dym_rows(ztz);
  double d;
  dym *vbeta = mk_robust_invert_spd(ztz,&d);
    /* vbeta = (ztz + d)^-1 where d is hopefully 0 or at least very small. */
  double sum_sqd_residuals;
  poster *po = AM_MALLOC(poster);
  po -> vbeta = vbeta;
  po -> beta_hat = mk_dym_times_dyv(vbeta,ztv);
  po -> nu = real_max(1.0,n - k);

  if ( !(dym_ref(vbeta,0,0) > -1e20 && dym_ref(vbeta,0,0) < 1e20) )
  {
    dym *dummy;
    double d2;
    fprintf_dym(stdout,"ztz",ztz,"\n");
    fprintf_dym(stdout,"vbeta",vbeta,"\n");

    printf("mk_non_informative_poster: bad value. Will call again.");
    really_wait_for_key();
    dummy = mk_robust_invert_spd(ztz,&d2);
    fprintf_dym(stdout,"dummy",dummy,"\n");
    my_error("mk_non_informative_poster: bad value.");
  }

  sum_sqd_residuals = compute_sum_sqd_residuals(ztz,ztv,vtv,po->beta_hat);
  if ( d > 0.0 )
    sum_sqd_residuals = n * 0.5; /* HACK! I fakely add noise if
                                    regression was singular */
  if ( sum_sqd_residuals < -1e-3 )
  {
    bool is_dym_sym_pos_def_and_nonsingular(dym *a);
    dym *id = mk_dym_mult(vbeta,ztz);
    fprintf_dym(stdout,"ztz",ztz,"\n");
    fprintf_dym(stdout,"vbeta",vbeta,"\n");
    fprintf_dym(stdout,"id",id,"\n");
    printf("vtv = %g\n",vtv);
    fprintf_dyv(stdout,"beta_hat",po->beta_hat,"\n");
    fprintf_dyv(stdout,"ztv",ztv,"\n");
    printf("Is it SPD? Answer = %d\n",is_dym_sym_pos_def_and_nonsingular(ztz));
    wait_for_key();
    my_error("-ve sum sqd residuals");
  }
  else if ( sum_sqd_residuals < 0.0 )
    sum_sqd_residuals = 0.0;

  po -> sigma_n_sqd = sum_sqd_residuals / real_max(1e-3,n - k);

  return(po);
}

poster *mk_known_sdev_non_informative_coeffs_poster(dym *ztz,dyv *ztv,
                                                   double vtv,
                                                   double noise_frac_prior)
{
  /* Use Eq 8.11 and 8.12 of Gelman et al with Sigma_y = nfp^2 I */
  double d;
  dym *vbeta = mk_robust_invert_spd(ztz,&d);
  poster *po = AM_MALLOC(poster);
  po -> vbeta = vbeta;
  po -> beta_hat = mk_dym_times_dyv(vbeta,ztv);
  po -> nu = 1e9; /* Infinitely many points seen as far as noise
                     estimation is concerned */
  po -> sigma_n_sqd = noise_frac_prior * noise_frac_prior;
  return(po);
}

poster *mk_poster(reginfo *ri,dym *ztz,dyv *ztv,double vtv,bprior *bp)
{
  poster *po = ( bp->non_informative ) ?
               mk_non_informative_poster(ztz,ztv,vtv) :
               ( bp -> known_sdev ) ?
               mk_known_sdev_non_informative_coeffs_poster(ztz,ztv,vtv,
                                                    bp->noise_frac_prior) :
               mk_informative_poster(ri,ztz,ztv,vtv,bp->nu0,
                                bp->noise_frac_prior,
                                bp->constant_coeff_prior_sdev,
                                bp->linear_coeff_prior_sdev,
                                bp->quadratic_coeff_prior_sdev);
  return(po);
}

void dyv_mid_and_range_helper(dyv *lo,dyv *hi,dyv **mid,dyv **range)
{
  int i;
  *mid = mk_dyv_plus(lo,hi);
  dyv_scalar_mult(*mid,0.5,*mid);
  *range = mk_dyv_subtract(hi,lo);
  dyv_scalar_mult(*range,0.5,*range);
  for ( i = 0 ; i < dyv_size(*range) ; i++ )
    dyv_set(*range,i,real_max(0.005,dyv_ref(*range,i)));
}

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
                    int regtype)
{
  reginfo *ri = AM_MALLOC(reginfo);
  dyv *xmid,*xrange,*ymid,*yrange;

  dyv_mid_and_range_helper(xlo,xhi,&xmid,&xrange);
  dyv_mid_and_range_helper(ylo,yhi,&ymid,&yrange);

  ri -> regtype = regtype;
  ri -> mid = mk_copy_dyv(xmid);
  ri -> range = mk_copy_dyv(xrange);
  ri -> ymid = mk_copy_dyv(ymid);
  ri -> yrange = mk_copy_dyv(yrange);

  free_dyv(xmid);
  free_dyv(ymid);
  free_dyv(xrange);
  free_dyv(yrange);

  return(ri);
}

void make_lo_hi_from_mid_range(dyv *mid,dyv *range,dyv **r_lo,dyv **r_hi)
{
  *r_lo = mk_dyv_subtract(mid,range);
  *r_hi = mk_dyv_plus(mid,range);
}

posters *mk_posters_from_reginfo_and_matrices(dym *ztz,dym *ztv,dym *vtv,
                                              reginfo *ri,bprior *bp)
{
  posters *po = AM_MALLOC(posters);
  int l;
  int usout_size = dym_cols(ztv);

  po -> num_posters = usout_size;
  po -> ri = mk_copy_reginfo(ri);

  for ( l = 0 ; l < usout_size ; l++ )
  {
    dyv *ztv_l = mk_dyv_from_dym_col(ztv,l);
    double vtv_l = dym_ref(vtv,l,l);

    po -> ps[l] = mk_poster(po->ri,ztz,ztv_l,vtv_l,bp);

    free_dyv(ztv_l);
  }

  return(po);
}

posters *mk_posters_from_reginfo(dym *x,dym *y,ivec *rows,
				 reginfo *ri,bprior *bp)
{
  int usout_size = dym_cols(y);
  int l;
  dym *ztz;
  dym *ztv;
  dym *vtv;
  int z_size;
  posters *po = NULL;

  ztz = mk_ztz_from_usins(ri,x,rows, /*xlo=*/NULL, /*xhi=*/NULL);
  z_size = dym_rows(ztz);

  ztv = mk_dym(z_size,usout_size);
  vtv = mk_zero_dym(usout_size,usout_size);

  for ( l = 0 ; l < usout_size ; l++ )
  {
    dyv *ztv_l;
    double vtv_l;

    make_reg_matrices_from_usins_and_usouts(ri,x,y,l,rows,
					    /*xlo=*/NULL, /*xhi=*/NULL,
					    &ztv_l,&vtv_l);

    copy_dyv_to_dym_col(ztv_l,ztv,l);
    dym_set(vtv,l,l,vtv_l);

    free_dyv(ztv_l);
  }

  po = mk_posters_from_reginfo_and_matrices(ztz,ztv,vtv,ri,bp);

  free_dym(ztz);
  free_dym(ztv);
  free_dym(vtv);
  return(po);
}

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
posters *mk_posters(dym *x,dym *y,
                    ivec *rows,
                    dyv *xlo,dyv *ylo,
                    dyv *xhi,dyv *yhi,
                    int regtype,bprior *bp)
{
  reginfo *ri = mk_reginfo(xlo,ylo,xhi,yhi,regtype);
  posters *po = mk_posters_from_reginfo(x,y,rows,ri,bp);
  free_reginfo(ri);
  return(po);
}

poster *posters_ref(posters *pos,int out_index)
{
  if ( out_index < 0 || out_index >= pos->num_posters )
    my_error("posters_ref");
  return(pos->ps[out_index]);
}

/* Draw the expected (maximum aposteriaiai) values of the
   coefficients.
     po --- created with mk_posters
     out_index --- the output we are predicting (0 <= out_index < dym_cols(y))
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
void make_posters_mean_cba(posters *po,int out_index,double *c,dyv **b,dym **a)
{
  poster *p = posters_ref(po,out_index);
  make_cba_from_beta(posters_reginfo(po),p->beta_hat,out_index,c,b,a);
}

/* Draw a random sample quadratic form from the posteriors of a fit.
     po --- created with mk_posters
     out_index --- the output we are predicting (0 <= out_index < dym_cols(y))
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
void make_posters_sample_cba(posters *po,int out_index,double *c,dyv **b,dym **a)
{
  poster *p = posters_ref(po,out_index);
  dyv *beta = mk_sample_beta(p->beta_hat,p->vbeta,p->nu,p->sigma_n_sqd);
  make_cba_from_beta(posters_reginfo(po),beta,out_index,c,b,a);
  free_dyv(beta);
}

/* Returns the estimated noise. Looks at the posterior distribution
   on the noise, and returns three values.

    The low limit of the 95% confidence interval on the noise
    The median value
    The high limit.

   These values are all in original variables (y-values) coordinate system.

   NOTE THEY ARE STANDARD DEVIATIONS NOT VARIANCES.
*/
void posters_noise_sdev(posters *po,int out_index,double *lo95,double *median,double *hi95)
{
  poster *p = posters_ref(po,out_index);
  reginfo *ri = posters_reginfo(po);
  double r = dyv_ref(ri->yrange,out_index);

  *lo95 = r * sqrt(scaled_inverse_chi_squared_cdf_inv(0.025,p->nu,p->sigma_n_sqd));
  *median = r * sqrt(scaled_inverse_chi_squared_cdf_inv(0.500,p->nu,p->sigma_n_sqd));
  *hi95 = r * sqrt(scaled_inverse_chi_squared_cdf_inv(0.975,p->nu,p->sigma_n_sqd));
}

#ifdef FANCY
void draw_cba(double c,dyv *b,dym *a)
{
  int steps = 128;
  int i;
  for ( i = 0 ; i < steps ; i++ )
  {
    double x1 = i * 512.0 / steps;
    double x2 = (i+1) * 512.0 / steps;
    double a00 = (a == NULL) ? 0.0 : dym_ref(a,0,0);
    double y1 = c + dyv_ref(b,0) * x1 + 0.5 * a00 * x1 * x1;
    double y2 = c + dyv_ref(b,0) * x2 + 0.5 * a00 * x2 * x2;
    ag_line(x1,y1,x2,y2);
  }
}
#endif

void free_poster(poster *p)
{
  free_dym(p->vbeta);
  free_dyv(p->beta_hat);
  AM_FREE(p,poster);
}

void free_reginfo(reginfo *ri)
{
  free_dyv(ri->mid);
  free_dyv(ri->range);
  free_dyv(ri->ymid);
  free_dyv(ri->yrange);
  AM_FREE(ri,reginfo);
}

reginfo *mk_copy_reginfo(reginfo *ri)
{
  reginfo *newri = AM_MALLOC(reginfo);
  newri->regtype = ri->regtype;
  newri->mid = mk_copy_dyv(ri->mid);
  newri->range = mk_copy_dyv(ri->range);
  newri->ymid = mk_copy_dyv(ri->ymid);
  newri->yrange = mk_copy_dyv(ri->yrange);
  return(newri);
}

/* Free p and all its subcomponents */
void free_posters(posters *p)
{
  int i;
  for ( i = 0 ; i < p->num_posters ; i++ )
  {
    poster *po = posters_ref(p,i);
    if ( po != NULL ) free_poster(po);
  }
  free_reginfo(p->ri);
  AM_FREE(p,posters);
}

void fprintf_poster(FILE *s,char *m1,poster *p,char *m2)
{
  char buff[1000];
  sprintf(buff,"%s->vbeta",m1);
  fprintf_dym(s,buff,p->vbeta,m2);
  sprintf(buff,"%s->beta_hat",m1);
  fprintf_dyv(s,buff,p->beta_hat,m2);
  sprintf(buff,"%s->nu",m1);
  fprintf_double(s,buff,p->nu,m2);
  sprintf(buff,"%s->sigma_n_sqd",m1);
  fprintf_double(s,buff,p->sigma_n_sqd,m2);
}

void fprintf_reginfo(FILE *s,char *m1,reginfo *p,char *m2)
{
  char buff[1000];

  sprintf(buff,"%s->regtype",m1);
  fprintf_string(s,buff,regtype_to_name(p->regtype),m2);
  sprintf(buff,"%s->mid",m1);
  fprintf_dyv(s,buff,p->mid,m2);
  sprintf(buff,"%s->range",m1);
  fprintf_dyv(s,buff,p->range,m2);
  sprintf(buff,"%s->ymid",m1);
  fprintf_dyv(s,buff,p->ymid,m2);
  sprintf(buff,"%s->yrange",m1);
  fprintf_dyv(s,buff,p->yrange,m2);
}

void fprintf_posters(FILE *s,char *m1,posters *p,char *m2)
{
  char buff[1000];
  int i;
  sprintf(buff,"%s->num_posters",m1);
  fprintf_int(s,buff,p->num_posters,m2);
  sprintf(buff,"%s->ri",m1);
  fprintf_reginfo(s,buff,p->ri,m2);
  for ( i = 0 ; i < p->num_posters ; i++ )
  {
    sprintf(buff,"%s->ps[%d]",m1,i);
    fprintf_poster(s,buff,p->ps[i],m2);
  }
}

dyv *reginfo_middle(reginfo *ri)
{
  return ri->mid;
}

dyv *reginfo_range(reginfo *ri)
{
  return ri->range;
}

dyv *reginfo_ymiddle(reginfo *ri)
{
  return ri->ymid;
}

dyv *reginfo_yrange(reginfo *ri)
{
  return ri->yrange;
}

double reginfo_dsqd(reginfo *ri,dyv *x1,dyv *x2)
{
  double result = 0.0;
  int i;
  dyv *range = reginfo_range(ri);
  for ( i = 0 ; i < dyv_size(x1) ; i++ )
  {
    double di = dyv_ref(x1,i) - dyv_ref(x2,i);
    double sdi = di / (2 * dyv_ref(range,i));
    result += sdi * sdi;
  }
  return(result);
}

int reginfo_num_outputs(reginfo *ri)
{
  return(dyv_size(ri->ymid));
}

dyv *posters_middle(posters *po)
{
  return  reginfo_middle(posters_reginfo(po));
}

int posters_num_outputs(posters *po)
{
  return reginfo_num_outputs(posters_reginfo(po));
}

int posters_regtype(posters *po)
{
  return reginfo_regtype(posters_reginfo(po));
}

#ifdef FANCY
/* The animate posters test must be run under X-windows.

   On each cycle the following things are drawn:

    The x-region (horizontal extent of the yellow rectangle)
    The y-region (vertical extent of it)
    The datapoints seen so far
    "num_samples" samples from the posterior of the linear or quadratic
       function that generated the data
    The MAP (also mean) estimate of the coefficients

    The console displays some of the numbers involved, including
    the MAP coefficients and the estimated noise.

      The user may then click a mouse button

     LEFT = Add a point here
     MID  = Change the range (yellow rectangle)
             First click: choose bottom left
             Next click:  choose top right
     RIGHT = Exit program

    "regtype" and "frac_noise_prior" have the meanings described
    in the comments for mk_posters.
*/
void animate_posters(int regtype,double frac_noise_prior,int num_samples)
{
  dym *x = mk_dym(0,1);
  dym *y = mk_dym(0,1);
  bool finished = FALSE;
  dyv *xlo = mk_dyv_1(0.0);
  dyv *xhi = mk_dyv_1(512.0);
  dyv *ylo = mk_dyv_1(0.0);
  dyv *yhi = mk_dyv_1(512.0);
  bprior *bp = mk_informative_bprior(
                    1.0, /* nu0 */
                    frac_noise_prior,
                    2.0, /* constant_coeff_prior_sdev */
                    5.0, /*linear_coeff_prior_sdev */
                    5.0 /*quadratic_coeff_prior_sdev */);
  ag_on("");

  while ( !finished )
  {
    int i;
    dym *a;
    dyv *b;
    double c;
    posters *po = NULL;
    double u,v;
    double lo95,median,hi95;
    int button;

    ag_on("");
    ag_set_pen_color(AG_YELLOW);
    ag_box(dyv_ref(xlo,0),dyv_ref(ylo,0),dyv_ref(xhi,0),dyv_ref(yhi,0));


    po = mk_posters(x,y,/*rows=*/NULL,xlo,ylo,xhi,yhi,regtype,bp);

    fprintf_posters(stdout,"po",po,"\n");

    ag_set_pen_color(AG_RED);

    for ( i = 0 ; i < num_samples ; i++ )
    {
      make_posters_sample_cba(po,0,&c,&b,&a);
      draw_cba(c,b,a);
      free_dyv(b);
      if ( a != NULL ) free_dym(a);
    }

    make_posters_mean_cba(po,0,&c,&b,&a);
    ag_set_pen_color(AG_BLUE);

    draw_cba(c,b,a);

    ag_set_pen_color(AG_BLACK);
    for ( i = 0 ; i < dym_rows(x) ; i++ )
      ag_disc(dym_ref(x,i,0),dym_ref(y,i,0),5.0);
    ag_set_pen_color(AG_WHITE);
    for ( i = 0 ; i < dym_rows(x) ; i++ )
      ag_disc(dym_ref(x,i,0),dym_ref(y,i,0),2.0);

    if ( regtype != LINEAR_REGTYPE )
      printf("MAP estimate of function is y = %g + %g * x + 0.5 %g x^2\n",
             c,dyv_ref(b,0),dym_ref(a,0,0));
    else
      printf("MAP estimate of function is y = %g + %g * x\n",
             c,dyv_ref(b,0));

    posters_noise_sdev(po,0,&lo95,&median,&hi95);

    printf("Noise: 95lo = %g, median = %g, 95hi = %g\n",lo95,median,hi95);
    free_dyv(b);
    if ( a != NULL ) free_dym(a);
    free_posters(po);

    button = ag_get_xy(&u,&v);
    if ( button == 3 )
      finished = TRUE;
    else if ( button == 2 )
    {
      ag_set_pen_color(AG_GREEN);
      ag_line(0.0,v,512.0,v);
      ag_line(u,0.0,u,512.0);
      dyv_set(xlo,0,u);
      dyv_set(ylo,0,v);
      ag_print(100.0,450.0,"Click Top Right of ROI");
      (void) ag_get_xy(&u,&v);
      dyv_set(xhi,0,u);
      dyv_set(yhi,0,v);
    }

    else if ( button == 1 )
    {
      add_row(x);
      dym_set(x,dym_rows(x)-1,0,u);
      add_row(y);
      dym_set(y,dym_rows(y)-1,0,v);
    }
  }
  free_dyv(xlo);
  free_dyv(xhi);
  free_dyv(ylo);
  free_dyv(yhi);
  free_bprior(bp);
  free_dym(x);
  free_dym(y);
}

/* Runs the animate test taking the following parameters from
   the command line:

     [regtype linear|circle|ellipse|quadratic]
     [fnp <double>] (the fraction of the y-range used as the prior
                     estimate of noise)
     [samples n]  (the number of sampled curves from the posterior) */

void posters_main(int argc,char *argv[])
{
  int regtype = regtype_from_args("regtype",argc,argv,QUADRATIC_REGTYPE);
  double frac_noise_prior = double_from_args("fnp",argc,argv,0.5);
  int num_samples = int_from_args("samples",argc,argv,20);
  animate_posters(regtype,frac_noise_prior,num_samples);
}
#endif

void explain_cba(char *message,double c,dyv *b,dym *a)
{
  char buff[100];

  printf("c%s = %g\n",message,c);

  sprintf(buff,"b%s",message);
  fprintf_dyv(stdout,buff,b,"\n");

  if ( a == NULL )
    printf("A%s = 0\n",message);
  else
  {
    sprintf(buff,"A%s",message);
    fprintf_dym(stdout,buff,a,"\n");
  }
}

void explain_posters_basic(posters *po, int out_index)
{
  dym *a;
  dyv *b;
  double c;

  make_posters_mean_cba(po,out_index,&c,&b,&a);
  explain_cba("",c,b,a);

  if ( a != NULL ) free_dym(a);
  free_dyv(b);
}

void explain_posters(posters *po, int out_index)
{
  printf("MAP fit is y = c + b^T x + 0.5 x^T A x where...\n");
  explain_posters_basic(po,out_index);
}

double coefficient_confidence(posters *po,int out_index,int coeff_num)
{
  poster *p = posters_ref(po,out_index);
  double nu,mu,sigsqd;
  double mean = dyv_ref(p->beta_hat,coeff_num);
  double top_of_interval;

  coefficient_marginal(p->beta_hat,p->vbeta,p->nu,p->sigma_n_sqd,
		       coeff_num,&nu,&mu,&sigsqd);

  top_of_interval = t_cdf_inv(0.99,mu,sigsqd,nu);
  return top_of_interval - mean;
}

void explain_posters_detailed(posters *po, int out_index)
{
  poster *p = posters_ref(po,out_index);
  int num_coeffs = dyv_size(p->beta_hat);
  double cmean,cmax,cdiff;
  dyv *bmean,*bmax,*bdiff;
  dym *amean,*amax,*adiff;
  dyv *coeffs_max = mk_dyv(num_coeffs);
  int i;
  reginfo *ri = posters_reginfo(po);

  for ( i = 0 ; i < num_coeffs ; i++ )
    dyv_set(coeffs_max,i,dyv_ref(p->beta_hat,i) +
                       coefficient_confidence(po,out_index,i));

  make_cba_from_beta(ri,p->beta_hat,out_index,&cmean,&bmean,&amean);
  make_cba_from_beta(ri,coeffs_max,out_index,&cmax,&bmax,&amax);

  cdiff = cmax - cmean;
  bdiff = mk_dyv_subtract(bmax,bmean);
  adiff = (amean==NULL)?NULL:mk_dym_subtract(amax,amean);

  cdiff = fabs(cdiff);
  for ( i = 0 ; i < dyv_size(bdiff) ; i++ )
    dyv_set(bdiff,i,fabs(dyv_ref(bdiff,i)));
  if ( adiff != NULL )
  {
    int j;
    for ( i = 0 ; i < dyv_size(bdiff) ; i++ )
      for ( j = 0 ; j < dyv_size(bdiff) ; j++ )
        dym_set(adiff,i,j,fabs(dym_ref(adiff,i,j)));
  }

  printf("MAP fit is y = c + b^T x + 0.5 x^T A x where...\n");
  explain_cba("",cmean,bmean,amean);

  printf("\n\nNow for the confidence intervals on the above...\n");

  explain_cba("_confidence",cdiff,bdiff,adiff);

  if ( amean != NULL ) free_dym(amean);
  free_dyv(bmean);
  if ( amax != NULL ) free_dym(amax);
  free_dyv(bmax);
  if ( adiff != NULL ) free_dym(adiff);
  free_dyv(bdiff);
  free_dyv(coeffs_max);
}

/* Same function as mk_reginfo() above except includes all
   datapoints mentioned in rows...not only those inside a hypercube.
   It's okay for rows to be NULL, denoting "use all datapoints" */
reginfo *mk_reginfo_from_all_datapoints(dym *x,dym *y,ivec *rows,int regtype)
{
  dyv *inlo,*inhi,*outlo,*outhi;
  reginfo *ri;
  make_sensible_vector_limits(x,rows,&inlo,&inhi);
  make_sensible_vector_limits(y,rows,&outlo,&outhi);
  ri = mk_reginfo(inlo,outlo,inhi,outhi,regtype);

  free_dyv(inlo);
  free_dyv(inhi);
  free_dyv(outlo);
  free_dyv(outhi);

  return(ri);
}

/* Same function as mk_posters() above except includes all
   datapoints...not only those inside a hypercube.
   rows may be NULL denoting "use all datapoints" */
posters *mk_posters_from_all_datapoints(dym *x,dym *y,ivec *rows,
                    int regtype,bprior *bp)
{
  reginfo *ri = mk_reginfo_from_all_datapoints(x,y,rows,regtype);
  posters *po = mk_posters_from_reginfo(x,y,rows,ri,bp);
  free_reginfo(ri);
  return(po);
}

/* The prediction of the out_index'th output attribute has a Prob Density Function
   the is a t-distribution. This function computes and returns the
   dof, mean, variance of this t-distribution */
void posters_predict_marginal(posters *po,int out_index,dyv *in,
                              double *r_dof,double *r_mean,double *r_variance)
{
  dyv *z = mk_z_from_usin(posters_reginfo(po),in);
  poster *poi = posters_ref(po,out_index);
  double scaled_y_mu;
  double scaled_y_var;
  double nu;
  double ymid = dyv_ref(posters_reginfo(po)->ymid,out_index);
  double yrange = dyv_ref(posters_reginfo(po)->yrange,out_index);

  predict_marginal(poi->beta_hat,poi->vbeta,poi->nu,poi->sigma_n_sqd,z,
                   &nu,&scaled_y_mu,&scaled_y_var);

  *r_dof      = nu;
  *r_mean     = ymid + scaled_y_mu * yrange;
  *r_variance = yrange * yrange * scaled_y_var;

  free_dyv(z);
}

/* Predicts the value of the given input. Also gives the
   upper and lower XXX-%ile bounds on this prediction
   where XXX is supplied by the called (e.g. 0.975 for 95%-ile bounds). */
double posters_predict(posters *po,int out_index,dyv *in,
                       double confidence_level, /* eg 0.975 for 95% bounds */
                       double *r_low_bound,double *r_high_bound)
{
  double dof,mean,var;

  posters_predict_marginal(po,out_index,in,&dof,&mean,&var);

  if ( r_high_bound != NULL )
  {
    *r_high_bound = t_cdf_inv(confidence_level,mean,sqrt(var),dof);
    if ( r_low_bound == NULL ) my_error("both should be null or non-null");
    *r_low_bound = 2 * mean - *r_high_bound;   /* Using symmetry of t-dist */
  }
  else if ( r_low_bound != NULL )
    my_error("hi null so should be lo null");

  return(mean);
}

/* What is the probability we'd have discovered such a high-magnitude
   sample mean if the true mean was really 0.0? */
double t_significance_level(double dof,double mean,double variance)
{
  return 2 * t_cdf(fabs(mean),0.0,sqrt(variance),dof) - 1.0;
}

void posters_coefficient_marginal(posters *po,int coeff_index,int out_index,
				  double *r_nu,double *r_mu,
				  double *r_sigma_squared)
{
  poster *p = posters_ref(po,out_index);
  coefficient_marginal(p->beta_hat,p->vbeta,p->nu,p->sigma_n_sqd,
		       coeff_index,r_nu,r_mu,r_sigma_squared);
}

/* What is the probability we'd have discovered such a high-magnitude
   fitted coefficient if the true coefficient was really 0.0? */
double posters_coefficient_significance_level(posters *po,
					      int coeff_index,int out_index)
{
  double dof,mean,variance;
  posters_coefficient_marginal(po,coeff_index,out_index,
			       &dof,&mean,&variance);
  return t_significance_level(dof,mean,variance);

}

/* What is the probability we'd have discovered such a high-magnitude
   fitted slope if the true slope was really 0.0? */
double posters_slope_significance_level(posters *po,
					int in_index,int out_index)
{
  return posters_coefficient_significance_level(po,in_index+1,out_index);
}

/* rows may be NULL meaning "use all rows" */
double posters_sum_sqd_residuals(posters *po,dym *x,dym *y,
				 ivec *rows,int out_index)
{
  int i;
  double ssr = 0.0;
  int num_rows = (rows == NULL) ? dym_rows(y) : ivec_size(rows);
  for ( i = 0 ; i < num_rows ; i++ )
  {
    int row = (rows == NULL) ? i : ivec_ref(rows,i);
    dyv *x_i = mk_dyv_from_dym_row(x,row);
    double true_y_i = dym_ref(y,row,out_index);
    double predict_y_i = posters_predict(po,out_index,x_i,0.9,NULL,NULL);
    ssr += real_square(true_y_i - predict_y_i);
    free_dyv(x_i);
  }
  return ssr;
}

/* rows may be NULL meaning "use all rows" */
dyv *mk_dyv_from_subset_of_dym_rows(dym *x,ivec *rows,int col)
{
  int num_rows = (rows == NULL) ? dym_rows(x) : ivec_size(rows);
  dyv *result = mk_dyv(num_rows);
  int i;

  for ( i = 0 ; i < num_rows ; i++ )
  {
    int row = (rows == NULL) ? i : ivec_ref(rows,i);
    dyv_set(result,i,dym_ref(x,row,col));
  }

  return result;
}

/* rows may be NULL meaning "use all rows" */
double posters_out_variance(dym *y,ivec *rows,int out_index)
{
  dyv *ysub = mk_dyv_from_subset_of_dym_rows(y,rows,out_index);
  double result = real_square(dyv_sdev(ysub));
  free_dyv(ysub);
  return result;
}

/* rows may be NULL meaning "use all rows" */
double posters_fraction_variance_explained(posters *po,dym *x,dym *y,
					   ivec *rows,int out_index)
{
  int num_rows = (rows == NULL) ? dym_rows(x) : ivec_size(rows);
  double sum_sqd_residuals = posters_sum_sqd_residuals(po,x,y,rows,out_index);
  double out_variance = posters_out_variance(y,rows,out_index);
  double result = 1.0 - sum_sqd_residuals /
                        real_max(1e-30,num_rows * out_variance);
  return real_max(0.0,real_min(1.0,result));
}

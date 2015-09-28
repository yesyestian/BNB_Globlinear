/*
   File:        stats.c
   Author:      Andrew W. Moore
   Created:     Sat Jun 17 10:17:26 EDT 1995
   Description: Cached cumulative denisty functions for fast statistics

   Copyright (C) 1995, Andrew W. Moore
*/

#include "stats.h"

/******* Simple little numerical statistics utilities.... *****/

/**** erf and lgamma don't exist on msdos. Here are our re-implementations,
      courtesy of Numerical Recipes in C.  Programmers should call these
      functions rather than erf and lgamma.
      - Mary Soon Lee, 2 Nov 1995
****/

/* This should be called instead of lgamma on Windows NT platforms.
   The code was taken from Numerical Recipes in C (where it is called
   gammln) , amended to use doubles not floats. 

   N.B. I have added a check to make sure the argument is in the correct
   range.
*/
double am_lgamma(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	if (xx <= 0.0)
	  my_error("lgamma given a non-positive argument");
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
	int i;
	double an,b,c,d,del,h;

	*gln=am_lgamma(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) my_error("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
	int n;
	double sum,del,ap;

	*gln=am_lgamma(a);
	if (x <= 0.0) {
		if (x < 0.0) my_error("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		my_error("a too large, ITMAX too small in routine gser");
		return;
	}
}

#undef ITMAX
#undef EPS
#undef FPMIN

double numrec_gammp(double a, double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) my_error("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
      
/* This should be called instead of erf on Windows NT platforms.  The
   code was taken from Numerical Recipes in C, amended to use doubles
   not floats. */
double am_erf(double x)
{
	return x < 0.0 ? -numrec_gammp(0.5,x*x) : numrec_gammp(0.5,x*x);
}



void fprintf_integ(FILE *s,char *m1,integ *it,char *m2)
{
  char buff[50];
  sprintf(buff,"%s -> integral",m1); fprintf_dyv(s,buff,it->integral,m2);
  sprintf(buff,"%s -> xlo",m1); fprintf_double(s,buff,it->xlo,m2);
  sprintf(buff,"%s -> xhi",m1); fprintf_double(s,buff,it->xhi,m2);
  sprintf(buff,"%s -> parameter",m1); fprintf_double(s,buff,it->parameter,m2);
  sprintf(buff,"%s -> constant",m1); fprintf_double(s,buff,it->constant,m2);
}

void free_integ(integ *it)
{
  free_dyv(it->integral);
  AM_FREE(it,integ);
}

integ *mk_integ(
    double (*h)(double parameter,double constant,double x),
    double xlo,
    double xhi,
    double parameter,
    double constant,
    int size
  )
/*
   Returns an it in which
   it->integral[i] = integal_from_xlo_to(x_lo + h*i) of h(parameter,x) dx
                     ------------------------------------------------
                     integal_from_xlo_to_x_hi of h(parameter,x) dx
*/
{
  integ *it = AM_MALLOC(integ);
  dyv *dig = mk_dyv(size);
  int i;
  double sum = 0.0;
  double last_pdf = 0.0;
  double delta = (xhi - xlo) / (size-1);

  if ( h(parameter,constant,xhi) > 1e-6 )
    my_error("Hmm... I was really hoping h(parameter,xhi) == 0");

  dyv_set(dig,0,0.0);

  for ( i = 1 ; i < size ; i++ )
  {
    double x = xlo + i * delta;
    double this_pdf = h(parameter,constant,x);
    if (i == 1) sum += delta * this_pdf;
    else        sum += delta * (this_pdf + last_pdf) / 2.0;
    dyv_set(dig,i,sum);
    last_pdf = this_pdf;  /* added 2/26/97  JGS */
  }

  dyv_scalar_mult(dig,1.0 / sum,dig);

  it -> integral = dig;
  it -> xlo = xlo;
  it -> xhi = xhi;
  it -> parameter = parameter;
  it -> constant = constant;

  return(it);
}

void get_index_and_fraction(
    double x,
    double lo,
    double hi,
    int size,
    int *r_index,
    double *r_fraction
  )
/*
   Assume we have a linear map, A --> B
      such that 0 --> lo
                (size-1) --> hi

    what values i (integer) and f (0 <= f < 1) would cause (i+f) --> x ?
*/
{
  double z = (size - 1) * (x - lo) / (hi - lo);
  *r_index = (int) floor(z);
  *r_fraction = z - *r_index;

  if ( *r_index == size-1 && *r_fraction < 1e-6 )
  {
    *r_fraction = 1.0;
    *r_index -= 1;
  }

  if ( *r_index < 0 || *r_index >= size-1 )
    my_error("oqaisnbcpwinc");
}

double two_interpolate(double y0,double y1,double fraction)
/*
   Given a linear map such that 0 ---> y0
                                1 ---> y1
    find result such that  f ---> result
*/
{
  return( (1 - fraction) * y0 + fraction * y1 );
}

double integ_cdf(integ *it,double x)
{
  double result;

  if ( x <= it->xlo )
    result = 0.0;
  else if ( x >= it->xhi )
    result = 1.0;
  else
  {
    int index;
    double fraction;
    get_index_and_fraction(x,it->xlo,it->xhi,
                           dyv_size(it->integral),&index,&fraction
                          );
    result = two_interpolate(dyv_ref(it->integral,index),
                             dyv_ref(it->integral,index+1),
                             fraction
                            );
  }

  return(result);
}

double integ_cdf_inv(integ *it,double prob)
/*
   Find x s.t. integ_cdf(in,x) = prob
*/
{
  double result;

  if ( prob < 0.0 || prob > 1.0 )
  {
    result = 0.0;
    printf("****** prob = %g (should be between 0 and 1)\n",prob);
    my_error("integ_cdf_inv: illegal prob");
  }
  else if ( prob == 0.0 )  /* == with doubles usually dodgy, but here
                                    harmless */
    result = it->xlo;
  else if ( prob == 1.0 )
    result = it->xhi;
  else
  {
    int lo = 0;
    int hi = dyv_size(it->integral);
    double result_index;

    while ( lo < hi-1 )
    {
      int mid = (lo + hi)/2;
      double value = dyv_ref(it->integral,mid);
      if ( value < prob )
        lo = mid;
      else
        hi = mid;
    }
    
    if ( hi - lo != 1 ) my_error("ouvbobvrlobfv");

      /* If c(x) is cdf of index x, then the local
         linear behaviour is c(x) = y1 + (x - x1) * (y2 - y1) / (x2 - x1) 

         If x1 = (double) lo, x2 = (double) hi, then x2 - x1 = 1.

         If we want x such that c(x) = prob then we need

             x = x1 + (prob - y1) / (y2 - y1)
      */
     
    result_index = ((double) lo) + 
                   ( (prob - dyv_ref(it->integral,lo)) /
                     (real_max(1e-10,dyv_ref(it->integral,hi) -
                                     dyv_ref(it->integral,lo)
                              )
                     )
                   );

    result = it->xlo + (it->xhi - it->xlo) * 
             result_index / ( dyv_size(it->integral) - 1 );
  }

  return(result);
}

void fprintf_minteg(FILE *s,char *m1,minteg *mit,char *m2)
{
  char buff[50];
  int i;

  sprintf(buff,"%s -> integ_array_size",m1);
  fprintf_int(s,buff,mit->integ_array_size,m2);
  sprintf(buff,"%s -> integ_size",m1);
  fprintf_int(s,buff,mit->integ_size,m2);

  sprintf(buff,"%s -> parameter_lo",m1);
  fprintf_double(s,buff,mit->parameter_lo,m2);
  sprintf(buff,"%s -> parameter_hi",m1);
  fprintf_double(s,buff,mit->parameter_hi,m2);

  for ( i = 0 ; i < mit->integ_array_size ; i++ )
  {
    sprintf(buff,"%s -> its[%d]",m1,i); 
    if ( mit->integ_array[i] == NULL )
      fprintf(s,"%s = NULL%s",buff,m2);
    else
      fprintf_integ(s,buff,mit->integ_array[i],m2);
  }
}

void free_minteg(minteg *mit)
{
  int i;
  for ( i = 0 ; i < mit->integ_array_size ; i++ )
  {
    if ( mit->integ_array[i] != NULL )
      free_integ(mit->integ_array[i]);
  }

  AM_FREE_ARRAY(mit->integ_array,integ_ptr,mit->integ_array_size);
  AM_FREE(mit,minteg);
}

minteg *mk_minteg(
    double (*h)(double parameter,double constant,double x),
    double (*compute_constant)(double parameter),
    double (*xlo)(double parameter,double constant),
    double (*xhi)(double parameter,double constant),
    double parameter_lo,
    double parameter_hi,
    int integ_array_size,
    int integ_size
  )
{
  minteg *mit = AM_MALLOC(minteg);
  int i;

  mit -> integ_array_size = integ_array_size;
  mit -> integ_size = integ_size;
  mit -> h = h;
  mit -> compute_constant = compute_constant;
  mit -> xlo = xlo;
  mit -> xhi = xhi;
  mit -> parameter_lo = parameter_lo;
  mit -> parameter_hi = parameter_hi;

  mit -> integ_array = AM_MALLOC_ARRAY(integ_ptr,integ_array_size);

  for ( i = 0 ; i < integ_array_size ; i++ )
    mit -> integ_array[i] = NULL;

  return(mit);
}

void maybe_compute_minteg_entry(minteg *mit,int i)
{
  if ( mit -> integ_array[i] == NULL )
  {
    double parameter = mit->parameter_lo +
                  i * (mit->parameter_hi - mit->parameter_lo) / 
                  (mit->integ_array_size - 1);
    double constant = mit->compute_constant(parameter);
    double xlo = mit->xlo(parameter,constant);
    double xhi = mit->xhi(parameter,constant);

    mit -> integ_array[i] = 
      mk_integ(mit->h,xlo,xhi,parameter,constant,mit->integ_size);
  }
}

double minteg_pdf(minteg *mit,double parameter,double x)
{
  double constant = mit->compute_constant(parameter);
  double result = mit->h(parameter,constant,x);
  return(result);
}

void get_pair_of_integs(
    minteg *mit,
    double parameter,
    integ **r_it0,
    integ **r_it1,
    double *r_fraction
  )
{
  int i;

  if ( parameter < mit->parameter_lo ||
       parameter > mit->parameter_hi
     )
  {
    printf("In one of the minteg statistics calls, you have asked for\n");
    printf("a value of the parameter (param = %g) that is out of range.\n",
           parameter
          );
    fprintf_minteg(stderr,"mit",mit,"\n");
    my_error("See above (minteg.c)");
  }

  get_index_and_fraction(parameter,mit->parameter_lo,mit->parameter_hi,
                         mit->integ_array_size,
                         &i,r_fraction
                        );

  maybe_compute_minteg_entry(mit,i);
  maybe_compute_minteg_entry(mit,i+1);

  *r_it0 = mit->integ_array[i];
  *r_it1 = mit->integ_array[i+1];
}

double minteg_cdf(minteg *mit,double parameter,double x)
{
  integ *it0;
  integ *it1;
  double fraction;
  double result;

  get_pair_of_integs(mit,parameter,&it0,&it1,&fraction);

  result = two_interpolate(integ_cdf(it0,x),integ_cdf(it1,x),fraction);

  return(result);
}

double minteg_cdf_inv(minteg *mit,double parameter,double prob)
/*
   Find x s.t. integ_cdf(in,x) = prob
*/
{
  integ *it0;
  integ *it1;
  double fraction;
  double result;

  get_pair_of_integs(mit,parameter,&it0,&it1,&fraction);

  result = 
    two_interpolate(integ_cdf_inv(it0,prob),integ_cdf_inv(it1,prob),fraction);

  return(result);
}

/*
  The parameter is the standard deviation of the gamma distribution.
  The distribution is assumed to have mean 1.

  A gamma dist with parameters alpha, beta has mean alpha/beta,
  variance alpha / beta^2.

  From this it is easy to see that a gamma dist with mean mu and std dev sigma
    has parameters alpha = (mu / sigma)^2 and beta = mu / sigma^2

  It can be shown that if X has a gamma dist with parameters alpha, beta
    then cX has a gamma dist with parameters alpha,beta/c

  This means that a gamma dist with mean mu and std dev sigma is
  distributed as X / mu where

    X has a gamma dist with mean 1 and std dev sigma/mu


  So that's why we will has sigma as our parameter, and mean 1.
  That means,

   alpha = 1/sigma^2 and beta = 1/sigma^2

  Note too that a gamma dist has p.d.f.

   f(x|a,b) = ( b^a / gamma(a)  ) x^(a-1) exp(- b x)

   logf = -b x + (a-1) log(x) + a log(b) - lgamma(a)

  
  For us, we'll use constant = a log(b) - lgamma(a). Note that since
    a = b = 1 / sigma^2, we have constant = b log(b) - lgamma(b)
*/   

double sd_gamma_h(double std_dev,double constant,double x)
{
  double a = 1 / (std_dev * std_dev);
  double b = a;
  double z = real_max(x,1e-6);
  double logf = -b * z + (a-1) * log(z) + constant;
  double result = exp(logf);
  return(result);
}

double sd_gamma_compute_constant(double std_dev)
{
  double b = 1 / (std_dev * std_dev);
  double c = b * log(b) - am_lgamma(b);
  return(c);
}

double sd_gamma_xlo(double std_dev,double constant)
{
  return(0.0);
}

double sd_gamma_xhi(double std_dev,double constant)
{
  double x = 1.0;
  bool ok = FALSE;
  int i;

  for ( i = 0 ; i < 50 && !ok ; i++ )
  {
    double pdf = sd_gamma_h(std_dev,constant,x);
    if ( pdf < 1e-6 )
      ok = TRUE;
    else
      x *= 2;
  }

  if ( !ok )
  {
    fprintf(stderr,"Convergence problem in minteg.c\n");
    fprintf(stderr,"Function is sd_gamma_xhi\n");
    fprintf(stderr,"std_dev = %g\n",std_dev);
    fprintf(stderr,"constant = %g\n",constant);
    fprintf(stderr,"i = %d\n",i);
    fprintf(stderr,"x = %g\n",x);
    fprintf(stderr,"ok = %d\n",ok);
  }
    

  return(x);
}

#define MAX_GAMMA_STD_DEV 50.0
#define GAMMA_ARRAY_SIZE 257

minteg *mk_sd_gamma_minteg(double std_dev_hi)
{
  int integ_array_size = GAMMA_ARRAY_SIZE;
  int integ_size = 12900; /* was 129 */
  minteg *mit = 
    mk_minteg(sd_gamma_h,
              sd_gamma_compute_constant,
              sd_gamma_xlo,
              sd_gamma_xhi,
              0.005,
              std_dev_hi,
              integ_array_size,
              integ_size
             );
  return(mit);
}

minteg *Sd_gamma_minteg = NULL;

minteg *maybe_make_sd_gamma_minteg()
{
  if ( Sd_gamma_minteg == NULL )
    Sd_gamma_minteg =
      mk_sd_gamma_minteg(MAX_GAMMA_STD_DEV);

  return(Sd_gamma_minteg);
}

double sd_gamma_pdf(double x,double std_dev)
/* pdf assuming mean 1, std dev "std_dev" */
{
  minteg *mit = maybe_make_sd_gamma_minteg();
  return(minteg_pdf(mit,std_dev,x));
}
 
double sd_gamma_cdf(double x,double std_dev)
/* cdf assuming mean 1, std dev "std_dev" */
{
  minteg *mit = maybe_make_sd_gamma_minteg();
  return(minteg_cdf(mit,std_dev,x));
}
 
double sd_gamma_cdf_inv(double prob,double std_dev)
/* cdf_inv assuming mean 1, std dev "std_dev" */
{
  minteg *mit = maybe_make_sd_gamma_minteg();
  return(minteg_cdf_inv(mit,std_dev,prob));
}
 
/*
   If X has gamma dist with params alpha,beta
  then
    beta/alpha X has gamma dist with params alpha,beta * alpha/beta =
                     gamma dist with params alpha,alpha
  which is a gamma dist with mean 1 and variance 1/alpha
  which is a gamma dist with mean 1 and std dev 1/sqrt(alpha)
  modified 9-5-95 by JS
  The random variable x*beta/alpha does have the distribution described.
  However, what we want is the distribution of x.  Reading the pdf values
  of x from gamma(alpha,beta) to be the same as x*beta/alpha from
  gamma(alpha,alpha) shrinks or expands the pdf function along the x
  axis by a factor of beta/alpha.  This causes the integral of the function
  to shrink or expand by that factor and must be corrected.  This correction
  only applies to the pdf, not the cdf.
*/
double gamma_pdf(double x,double alpha,double beta)
{
/*  return(sd_gamma_pdf(x * beta/alpha,1.0 / sqrt(alpha))); */
  return((beta/alpha)*sd_gamma_pdf(x * beta/alpha,1.0 / sqrt(alpha)));
}

double gamma_cdf(double x,double alpha,double beta)
{
  return(sd_gamma_cdf(x * beta/alpha,1.0 / sqrt(alpha)));
}

double gamma_cdf_inv(double prob,double alpha,double beta)
{
  return(alpha/beta * 
         sd_gamma_cdf_inv(prob,1.0 / sqrt(alpha))
        );
}

/*
  The parameter is the square root of degrees of freedom of the t distribution.
  The distribution is assumed to have mean 0 and variance 1.

  A t dist with nu d.o.f's has the following pdf

  f(x|nu) = gamma((nu+1)/2)
            ----------------------- ( 1 + x^2/nu )^(-(nu+1)/2)
            sqrt(nu pi) gamma(nu/2)

  Thus logf = constant - (nu+1)/2 * log(1 + x^2 / nu)

  where constant = lgamma((nu+1)/2) - 1/2 log(nu) - 1/2 log(pi) - lgamma(nu/2)

  where nu = parameter^2
*/   

double roodof_t_h(double root_nu,double constant,double x)
{
  double nu = root_nu * root_nu;
  double logf = constant - 0.5 * (nu+1.0) * log(1.0 + x * x / nu);
  double result = exp(logf);
  return(result);
}

double roodof_t_compute_constant(double root_nu)
{
  double nu = root_nu * root_nu;
  double constant = am_lgamma((nu+1.0)/2) - am_lgamma(nu/2.0) - 
                    (log(nu) + log(PI)) / 2.0;
  return(constant);
}

double roodof_t_xhi(double root_nu,double constant)
{
  double x = 1.0;
  bool ok = FALSE;
  int i;

  for ( i = 0 ; i < 50 && !ok ; i++ )
  {
    double pdf = roodof_t_h(root_nu,constant,x);
    if ( pdf < 1e-6 )
      ok = TRUE;
    else
      x *= 2;
  }

  if ( !ok ) my_error("tttttklaowidn");

  return(x);
}

double roodof_t_xlo(double root_nu,double constant)
{
  return(-roodof_t_xhi(root_nu,constant));
}

#define MAX_DOF 400.0
/* Modified 9-6-95 by JS.  Integ_array_size and integ_size used to be 
 * 65 and 129.  That is probably too coarse.  For example, compare the
 * cdf generated by calls to t_cdf with that generated by integrating
 * the calls to t_pdf with mu=0 precision=1 dof=1 and 1000 samples
 * ranging from dist->xmin to dist->xmax
 */
minteg *mk_roodof_t_minteg()
{
  int integ_array_size = 650;  /* was 65  */
  int integ_size = 1299;       /* was 129 */
  minteg *mit = 
    mk_minteg(roodof_t_h,
              roodof_t_compute_constant,
              roodof_t_xlo,
              roodof_t_xhi,
              1e-3,
              sqrt(MAX_DOF),
              integ_array_size,
              integ_size
             );
  return(mit);
}

minteg *Roodof_t_minteg = NULL;

minteg *maybe_make_roodof_t_minteg()
{
  if ( Roodof_t_minteg == NULL )
    Roodof_t_minteg =
      mk_roodof_t_minteg();

  return(Roodof_t_minteg);
}

double roodof_t_pdf(double x,double root_nu)
/* pdf assuming mean 0, std dev 1, nu = root_nu^2 */
{
  minteg *mit = maybe_make_roodof_t_minteg();
  double result;

  result = minteg_pdf(mit,root_nu,x);
  return(result);
}
 
double roodof_t_cdf(double x,double root_nu)
/* cdf assuming mean 0, std dev 1, nu = root_nu^2 */
/* We know it's symmetric, so we increase accuracy by forcing symmetry */
{
  minteg *mit = maybe_make_roodof_t_minteg();
  double r1 = minteg_cdf(mit,root_nu,x);
  double r2 = 1.0 - minteg_cdf(mit,root_nu,- x);
  return(0.5 *(r1 + r2));
}
 
double roodof_t_cdf_inv(double prob,double root_nu)
/* cdf_inv assuming mean 0, std dev 1, nu = root_nu^2 */
/* We know it's symmetric, so we increase accuracy by forcing symmetry */
{
  minteg *mit = maybe_make_roodof_t_minteg();
  double r1 = minteg_cdf_inv(mit,root_nu,prob);
  double r2 = -minteg_cdf_inv(mit,root_nu,1.0 - prob);
  return(0.5 * (r1 + r2));
}

/* modified 9-5-95 by JS, see comments on gamma_pdf */
double t_pdf(double x,double mu,double std_dev,double dof)
{
  if ( dof > MAX_DOF ) dof = MAX_DOF; /* Won't make any difference. At this
                                         magnitude it's a gaussian */
  return(roodof_t_pdf((x - mu)/std_dev,sqrt(dof))/std_dev);
}

double t_cdf(double x,double mu,double std_dev,double dof)
{
  if ( dof > MAX_DOF ) dof = MAX_DOF; /* Won't make any difference. At this
                                         magnitude it's a gaussian */
  return(roodof_t_cdf((x - mu)/std_dev,sqrt(dof)));
}

double t_cdf_inv(double prob,double mu,double std_dev,double dof)
{
  if ( dof > MAX_DOF ) dof = MAX_DOF; /* Won't make any difference. At this
                                         magnitude it's a gaussian */
  return(mu + std_dev * roodof_t_cdf_inv(prob,sqrt(dof)));
}

void close_statistics()
{
  if ( Sd_gamma_minteg != NULL )
  {
    free_minteg(Sd_gamma_minteg);
    Sd_gamma_minteg = NULL;
  }
  if ( Roodof_t_minteg != NULL )
  {
    free_minteg(Roodof_t_minteg);
    Roodof_t_minteg = NULL;
  }
}

void stats_malloc_report()
{
  if ( Sd_gamma_minteg != NULL || Roodof_t_minteg != NULL )
  {
   fprintf(stdout,"#\n#IMPORTANT IMPORTANT IMPORTANT ************\n# The cached data structures used by DAMUT statistics are still\n");
   fprintf(stdout,"# allocated. Call close_statistics() [prototype in stats.h] if\n");
   fprintf(stdout,"# you wish to deallocate them.\n");
  }
}

double beta_that_forces_requested_median(double alpha,double median)
/*
   If Gamma(alpha,1) has median K
   Then Gamma(alpha,beta) has median K/beta.

   We want median = K/beta, so we put beta = K/median where
     K = Median(Gamma(alpha,1))
*/
{
  double k = gamma_cdf_inv(0.5,alpha,1.0);
  double result = k / median;
  
  {
/*
    double check = gamma_cdf_inv(0.5,alpha,result);
    printf("alpha=%g, goal_median=%g, chosen_beta=%g,actual_median=%g\n",
           alpha,median,result,check
          );
    wait_for_key();
*/
  }

  return(result);
}

/**** Stuff follows from Numerical Recipes In C
      (Press, Teukolsky Vetterling and Flannery:
        Cambridge University Press)

  We are using the following functions:
    gammln
    gcf
    gser
    gammq
    chsone

  See the book (and BUY THE BOOK if you don't have it)
  for details. The book is very very good.

  Note: All floats turned into doubles. nrerror replaced with my_error.
*****/
double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double gammq(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) my_error("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

void chsone(double bins[], double ebins[], int nbins, int knstrn, double *df,
	double *chsq, double *prob)
{
	double gammq(double a, double x);
	int j;
	double temp;

	*df=nbins-knstrn;
	*chsq=0.0;
	for (j=1;j<=nbins;j++) {
		if (ebins[j] <= 0.0) my_error("Bad expected number in chsone");
		temp=bins[j]-ebins[j];
		*chsq += temp*temp/ebins[j];
	}
	*prob=gammq(0.5*(*df),0.5*(*chsq));
}

/********** End Numerical Recipes in C stuff *********/

/* As below, except we insist as a pre-condition that
   no entry in hypothesized_dist has an expected value of
   zero.
*/
double chi_squared_prob_helper(dyv *actual_dist,dyv *hypothesized_dist,int dof)
{
  int size = dyv_size(actual_dist);
  int nr_size = size + 1;
  double *bins = AM_MALLOC_ARRAY(double,nr_size);
  double *ebins = AM_MALLOC_ARRAY(double,nr_size);
  int i;
  int knstrn = size - dof;
  double df,chsq,prob;
 
  if ( size != dyv_size(hypothesized_dist) )
    my_error("chi_squared_prob");

  for ( i = 0 ; i < size ; i++ )
  {
    bins[i+1] = (double) dyv_ref(actual_dist,i);
    ebins[i+1] = (double) dyv_ref(hypothesized_dist,i);
  }

  chsone(bins,ebins,size,knstrn,&df,&chsq,&prob);

  AM_FREE_ARRAY(bins,double,nr_size);
  AM_FREE_ARRAY(ebins,double,nr_size);

  return prob;
}

/* 
   PRE: size of actual_dist is same as size of hypothesized_dist.
        Any entry in which hypothesized_dist has a value of
        zero must have an actual_dist value of zero (i.e.
          forall i, hy_dist[i]==0 => ac_dist[i] == 0

   Given two distributions represented as histograms 
   (actual_dist and hypothesized_dist), how much evidence is there that they are
   from the same distribution? 
   Note that these things must be counts. Each element of actual_dist must
   be an integer. Each element of hypothesized_dist may be non-integer
   because we're talking expected counts there.

   The prob returned by this function answers that question
   using a standard chi-squared test. If it is low (e.g. < 0.05), then it is
   unlikely that they are the same. 

   The "dof" parameter is the mysterious "Degrees Of Freedom" that haunts
   any use of the word "Chi". 
   
       If it is possible for any entry in the dist
       to take any value, then set dof==size.

       If the sum of values is constrained to a certain value
       then set dof==size-1.

       If there are more constraints than that, then subtract
       more from size.
*/
double chi_squared_prob(dyv *actual_dist,dyv *hypothesized_dist,int dof)
{
  double result = -1.0;
  double min_hyp_dist = dyv_min(hypothesized_dist);
  if ( min_hyp_dist < 0.0 )
    my_error("chi_squared_prob: -ve count in hypothesized_dist");
  else if ( min_hyp_dist > 0.0 )
    result = chi_squared_prob_helper(actual_dist,hypothesized_dist,dof);
  else
  {
    dyv *copy_ad = mk_dyv(0);
    dyv *copy_hd = mk_dyv(0);
    int i;
    for ( i = 0 ; i < dyv_size(actual_dist) ; i++ )
    {
      if ( dyv_ref(hypothesized_dist,i) > 0.0 )
      {
        add_to_dyv(copy_ad,dyv_ref(actual_dist,i));
        add_to_dyv(copy_hd,dyv_ref(hypothesized_dist,i));
        dof -= 1;
      }
      else if ( dyv_ref(actual_dist,i) > 0.0 )
        my_error("chi_squared_prob: actual_dist value must be zero if hyp dist value is zero");
    }
    dof = int_max(2,dof);
    result = chi_squared_prob_helper(copy_ad,copy_hd,dof);
    free_dyv(copy_ad);
    free_dyv(copy_hd);
  }

  return result;
}

/******************** Now some utilities that can help with k-fold-crossvalidation.. *******/

/* Define an ivec as being a natset if and only if it 
   contains a subset of the natural numbers in a
   strictly increasing order. */

/* PRE: a and b are natsets.
   PRE: b is a subset of a

   POST: returns a natset c such that C intersect B = empty
                             C union B     = A

       (i.e. n is in result if and only if n is in
             A but not in B)

   This function is useful if you have a set of rows (record numbers)
   to use as a test set (in b) and you have all the rows for the test
   and train set in a, and you want to get the set of remaining
   rows to use as a training set. */
ivec *mk_ivec_set_difference(ivec *a,ivec *b)
{
  int a_size = ivec_size(a);
  int b_size = ivec_size(b);
  int c_size = a_size - b_size;
  int a_ptr = 0;
  int b_ptr = 0;
  int c_ptr = 0;
  ivec *c = mk_ivec(c_size);

  for ( a_ptr = 0 ; a_ptr < a_size ; a_ptr++ )
  {
    int a_val = ivec_ref(a,a_ptr);
    if ( b_ptr == b_size || ivec_ref(b,b_ptr) > a_val )
    {
      ivec_set(c,c_ptr,a_val);
      c_ptr += 1;
    }
    else if ( ivec_ref(b,b_ptr) == a_val )
      b_ptr += 1;
    else
      my_error("mk_ivec_set_difference: b not a subset of a");
  }

  return c;
}

/* If you've got a dataset where you are training and testing
   using the set of rows specified in "train_and_test_rows" and you
   wish to do k-fold cross-validation (k == num_folds) then what 
   subset of train_and_test_rows should you use as the train set
   fold_num'th fold? And which should be the test set? 
   This function tells you.
   It uses a deterministic shuffle and then selects the given rows. 

   PRE: train_and_test_rows is a natset (see above) (it means
        strictly increasing set of natural numbers) (in which case
        num_rows is ignored).

        OR train_and_test_rows may be NULL, indicating "use all
        rows in the set { 0 , 1 , ... num_rows-1 }

   POST: *r_train_set is the fold_num'th training set
         *r_test_set is the fold_num'th test set.
         (Note that the union of these two is the original set of rows).
         The resulting sets are both returned in strictly increasing
         order.

    Note that if you called this function repeatedly with the same number
    of folds, but you varied fold_num from 0 up to num_folds on each
    call, then the union of all the resulting test_sets would contain the
    same set as train_and_test_rows. But all pairs of resulting test_sets 
    would have an empty intersection.
*/
void make_kfold_rows(ivec *train_and_test_rows,int num_rows,int num_folds,int fold_num,
                     ivec **r_train_rows,ivec **r_test_rows)
{
  int save_seed = int_random(300000);
  ivec *srows = (train_and_test_rows==NULL) ? mk_identity_ivec(num_rows) :
                                              mk_copy_ivec(train_and_test_rows);
  int srows_size = ivec_size(srows);
  int start_i = (int) floor(fold_num * num_rows / (double) num_folds);
  int end_i = (int) floor((fold_num+1) * num_rows / (double) num_folds);
  int i;
  *r_train_rows = mk_ivec(srows_size - (end_i - start_i));
  *r_test_rows = mk_ivec(end_i - start_i);
  am_srand(12345);
  shuffle_ivec(srows);

  for ( i = 0 ; i < srows_size ; i++ )
  {
    ivec *update_me = (i >= start_i && i < end_i) ? *r_test_rows : *r_train_rows;
    int update_index = (i < start_i) ? i :
                       (i < end_i) ? i - start_i : i - (end_i - start_i);
    ivec_set(update_me,update_index,ivec_ref(srows,i));
  }
  free_ivec(srows);
  am_srand(save_seed);
  ivec_sort(*r_train_rows,*r_train_rows);
  ivec_sort(*r_test_rows,*r_test_rows);
}

/* If you've got a dataset where you are training and testing
   using the set of rows specified in "train_and_test_rows" and you
   wish to make a test-set of size "num_test_rows" and a training-set
   of size "num_training_row (= num_tain_and_test - num_test)"
   this returns the training and test rows you should use, obtained by
   shuffling train_and_test.

   PRE: train_and_test_rows is a natset (see above) (it means
        strictly increasing set of natural numbers) (in which case
        num_rows is ignored).

        OR train_and_test_rows may be NULL, indicating "use all
        rows in the set { 0 , 1 , ... num_rows-1 }

   POST: *r_train_set is the training set
         *r_test_set is the test set.
         (Note that the union of these two is the original set of rows).
         The resulting sets are both returned in strictly increasing
         order.
*/
void make_train_and_test_rows(ivec *train_and_test_rows,
			      int num_rows,int num_test_rows,
			      ivec **r_train_rows,ivec **r_test_rows)
{
  int save_seed = int_random(300000);
  ivec *srows = (train_and_test_rows==NULL) ? 
                mk_identity_ivec(num_rows) :
                mk_copy_ivec(train_and_test_rows);
  int srows_size = ivec_size(srows);
  int i;
  *r_train_rows = mk_ivec(srows_size - num_test_rows);
  *r_test_rows = mk_ivec(num_test_rows);
  am_srand(12345);
  shuffle_ivec(srows);

  for ( i = 0 ; i < num_test_rows ; i++ )
    ivec_set(*r_test_rows,i,ivec_ref(srows,i));

  for ( i = num_test_rows ; i < srows_size ; i++ )
    ivec_set(*r_train_rows,i-num_test_rows,ivec_ref(srows,i));

  free_ivec(srows);
  am_srand(save_seed);
  ivec_sort(*r_train_rows,*r_train_rows);
  ivec_sort(*r_test_rows,*r_test_rows);
}

/* Makes a dym consisting of a subset of the rows in x. The members of
   of the subset are those rows mentioned in "rows".
   Result will this have "ivec_size(rows)" rows and dym_cols(x) columns */
dym *mk_dym_from_subset_of_rows(dym *x,ivec *rows)
{
  int num_rows = ivec_size(rows);
  int i;
  dym *result = mk_dym(num_rows,dym_cols(x));

  for ( i = 0 ; i < num_rows ; i++ )
  {
    int row = ivec_ref(rows,i);
    dyv *vec = mk_dyv_from_dym_row(x,row);
    copy_dyv_to_dym_row(vec,result,i);
    free_dyv(vec);
  }
  return result;
}

/* Creates two new dyms, both with the same number of columns as the
   input x. *r_test will have "num_test_rows" rows and *r_train will
   have "dym_rows(x) - num_test_rows" rows. *r_test will consist of 
   "num_test_rows" rows of x selected at random without replacement. 
   *r_train will contain the rest. The rows will appear in the same order
   as they did in the original. */
void break_dym_into_train_and_test(dym *x,int num_test_rows,
				   dym **r_train,dym **r_test)
{
  ivec *train_rows;
  ivec *test_rows;
  make_train_and_test_rows(/* train_and_test_rows = */ NULL /* denoting ALL */,
			   dym_rows(x),num_test_rows,&train_rows,&test_rows);

  *r_train = mk_dym_from_subset_of_rows(x,train_rows);
  *r_test = mk_dym_from_subset_of_rows(x,test_rows);

  free_ivec(train_rows);
  free_ivec(test_rows);
}


/*
   File:         hrect.cpp
   Author:       Andrew W. Moore
   Created:      Jan 13th, 1998
   Description:  The representation of motor control hrects.
*/

#include "hrect.h"
#include "./draw/amgr.h"

/* A hrect represents a hyperrectangle in k-dimensional space.
   This is easily represented as two vectors: the values of
   all the low coordinates and high coordinates. For example
   this 2-d rectangle

  y = 7        +-------+
               |       |
  y = 4        +-------+

               ^       ^
               x = -1  x = 3

   can be represented by lo = ( -1 , 4 )
                         hi = (  3 , 7 )
*/

/* makes a hyper rectangle with lowest value in ith dimension
   dyv_ref(lo,i) and highest value dyv_ref(hi,i) */
hrect *mk_hrect(dyv *lo,dyv *hi)
{
  hrect *hr = AM_MALLOC(hrect);
  hr -> lo = mk_copy_dyv(lo);
  hr -> hi = mk_copy_dyv(hi);

  am_assert( dyv_size(hr->lo) == dyv_size(hr->hi) );

  return hr;
}

void copy_hrect(hrect *src,hrect *dest)
{
  copy_dyv(src->lo,dest->lo);
  copy_dyv(src->hi,dest->hi);
}

hrect *mk_copy_hrect(hrect *hr)
{
  return mk_hrect(hr->lo,hr->hi);
}

/* This parses two strings (which must be space-separated ASCII
   representations of numbers). It turns the first into the low
   vector for the hyperrectnagle and the second into the high vector.
   Clearly both vectors must be the same length (i.e. have the same
   number of numbers). The number of numbers is the dimensionality of
   the hyper-rectangle */
hrect *mk_hrect_from_strings(char *lo_string,char *hi_string)
{
  dyv *lo = mk_dyv_from_string_simple(lo_string);
  dyv *hi = mk_dyv_from_string_simple(hi_string);
  hrect *hr = mk_hrect(lo,hi);

  free_dyv(lo);
  free_dyv(hi);

  return hr;
}

/* The dimensionality of the hyper-rectangle */
int hrect_size(const hrect *hr)
{
  return dyv_size(hr->lo);
}

/* MAKES and returns the center of the hyper-rectangle. (This is
   of course just the mean of lo and hi) */
dyv *mk_hrect_middle(hrect *hr)
{
  dyv *sum = mk_dyv_plus(hr->lo,hr->hi);
  dyv_scalar_mult(sum,0.5,sum);
  return sum;
}

/* If x is inside the hyper-rectangle it is unchanged.
   If x is outside the hyper-rectangle it is changed to
   become the closest point to x that is in (or on the surface
   of) the hyper-rectangle */
void clip_to_inside_hrect(const hrect *hr, dyv *x)
{
  int i;
  for ( i = 0 ; i < dyv_size(x) ; i++ )
  {
    double xi = dyv_ref(x,i);
    if ( xi < dyv_ref(hr->lo,i) )
      dyv_set(x,i,dyv_ref(hr->lo,i));
    else if ( xi > dyv_ref(hr->hi,i) )
      dyv_set(x,i,dyv_ref(hr->hi,i));
  }
}

/* Returns TRUE if and only if query_point is inside or on
   the edge of the hyper-rectangle */
bool is_in_hrect(const hrect *hyp, const dyv *query_point)
{
  return dyv_weakly_dominates(hyp->hi,query_point) &&
         dyv_weakly_dominates(query_point,hyp->lo);
}

/* am_frees the hyper-rect. (And, in accordance with AM conventions,
   it frees all the substructures). */
void free_hrect(hrect *hyp)
{
  free_dyv(hyp->lo);
  free_dyv(hyp->hi);
  AM_FREE(hyp,hrect);
}

/* Prints using the AM-FPRINT-CONVENTION described in motutils.h */
void fprint_hrect(FILE *s,char *s1,char *s2,hrect *hyp)
{
  char prefix[1000];
  sprintf(prefix,"%s -> %s",s1,s2);

  fprint_dyv(s,prefix,"lo",hyp->lo);
  fprint_dyv(s,prefix,"hi",hyp->hi);
}

/* Prints using the AM-FPRINTF-CONVENTION */
void fprintf_hrect(FILE *s,char *m1,hrect *hyp,char *m2)
{
  char buff[500];
  sprintf(buff,"%s -> lo",m1);
  fprintf_dyv(s,buff,hyp->lo,m2);
  sprintf(buff,"%s -> hi",m1);
  fprintf_dyv(s,buff,hyp->hi,m2);
}

/* Suppose you want to draw using ag_dot, ag_line etc (defined in
   [x]damut/amgr.h). But suppose you want it to be the case
   that ag_dot(p1,p2) will draw on the bottom left hand corner
   of the graphics display, where
      p1 = lo value of i1'th dimension of hyper-rect
      p2 = lo value of i2'th dimension of hyper-rect

 And suppose you want it to be the case
   that ag_dot(q1,q2) will draw on the top right hand corner
   of the graphics display, where
      q1 = hi value of i1'th dimension of hyper-rect
      q2 = hi value of i2'th dimension of hyper-rect

   Then this is the function for you! */
void hrect_set_ag_frame(hrect *hr,int i1,int i2)
{
  set_ag_frame(dyv_ref(hr->lo,i1),dyv_ref(hr->lo,i2),
               dyv_ref(hr->hi,i1),dyv_ref(hr->hi,i2));
}

/* As above, except shrinks the frame so that on all four sides
   of the ag_ window there's a margin of width 5% window width. */
void hrect_set_ag_frame_with_border(hrect *hr,int i1,int i2)
{
  double w1 = hrect_width_ref(hr,i1);
  double w2 = hrect_width_ref(hr,i2);
  set_ag_frame(dyv_ref(hr->lo,i1)-w1/20.0,dyv_ref(hr->lo,i2)-w2/20.0,
               dyv_ref(hr->hi,i1)+w1/20.0,dyv_ref(hr->hi,i2)+w2/20.0);
}

void ag_hrect(hrect *hr)
{
  ag_rectangle(dyv_ref(hr->lo,0),dyv_ref(hr->lo,1),
         dyv_ref(hr->hi,0),dyv_ref(hr->hi,1));
}

/* rows may be NULL denoting "use all rows" */
hrect *mk_hrect_bounding_dym_rows(dym *x,ivec *rows)
{
  dyv *lo,*hi;
  hrect *hr;
  make_vector_limits(x,rows,&lo,&hi);
  hr = mk_hrect(lo,hi);
  free_dyv(lo);
  free_dyv(hi);
  return hr;
}

hrect *mk_hrect_bounding_dyv_array(const dyv_array *da)
{
  dyv *lo,*hi;
  hrect *hr;
  make_vector_limits_from_dyv_array(da,&lo,&hi);
  hr = mk_hrect(lo,hi);
  free_dyv(lo);
  free_dyv(hi);
  return hr;
}

int hrect_widest_dim(hrect *hr)
{
  int result = -1;
  double wmax = -1.0;
  int i;
  for ( i = 0 ; i < hrect_size(hr) ; i++ )
  {
    double w = hrect_width_ref(hr,i);
    if ( i == 0 || w > wmax )
    {
      wmax = w;
      result = i;
    }
  }
  return result;
}

double hrect_middle_ref(hrect *hr,int index)
{
  return 0.5 * ( dyv_ref(hr->hi,index) + dyv_ref(hr->lo,index) );
}

double hrect_width_ref(hrect *hr,int index)
{
  return dyv_ref(hr->hi,index) - dyv_ref(hr->lo,index);
}

/* x_minlik and x_maxlik may be NULL.

 *distmin is set to be min_x_in_hr( [x - center]^2 )
 *distmax is set to be max_x_in_hr( [x - center]^2 )

   if x_minlik is non-null, it gets filled with
     argmax_x_in_hr ( [x - center]^2 )
   if x_maxlik is non-null, it gets filled with
     argmin_x_in_hr ( [x - center]^2 )

  Note that "x_minlik" means the point least likely to have been created
  by center, which is why it is defined as an argmax.
*/
void closest_and_furthest_hrect_distance(hrect *hr,dyv *center,
                                         double *distmin,double *distmax,
                                         dyv *x_minlik,dyv *x_maxlik)
{
  int i;
  double sum_max_dsqd = 0.0;
  double sum_min_dsqd = 0.0;

  for ( i = 0 ; i < dyv_size(center) ; i++ )
  {
    double lo = hrect_lo_ref(hr,i);
    double hi = hrect_hi_ref(hr,i);
    double x = dyv_ref(center,i);
    double dmin,dmax;

    if ( x <= lo )
    {
      dmin = lo - x;
      dmax = hi - x;
      if ( x_minlik != NULL ) dyv_set(x_minlik,i,hi);
      if ( x_maxlik != NULL ) dyv_set(x_minlik,i,lo);
    }
    else if ( x >= hi )
    {
      dmin = x - hi;
      dmax = x - lo;
      if ( x_minlik != NULL ) dyv_set(x_minlik,i,lo);
      if ( x_maxlik != NULL ) dyv_set(x_minlik,i,hi);
    }
    else
    {
      dmin = 0.0;
      dmax = ( x < (lo+hi)/2.0 ) ? hi-x : x-lo;
      if ( x_minlik != NULL ) dyv_set(x_minlik,i,(x<(lo+hi)/2) ? hi : lo);
      if ( x_maxlik != NULL ) dyv_set(x_minlik,i,x);
    }

    sum_max_dsqd += dmax * dmax;
    sum_min_dsqd += dmin * dmin;
  }

  *distmax = sqrt(sum_max_dsqd);
  *distmin = sqrt(sum_min_dsqd);
}

/* Returns a new hrect translated by subtracting the vector
   x from the lower limits and the higher limits of hr */
hrect *mk_hrect_subtract_dyv(hrect *hr,dyv *x)
{
  hrect *hr_trans = mk_copy_hrect(hr);
  dyv_subtract(hr->lo,x,hr_trans->lo);
  dyv_subtract(hr->hi,x,hr_trans->hi);
  return hr_trans;
}

hrect *mk_zero_hrect(int size)
{
  dyv *zeroes = mk_zero_dyv(size);
  hrect *hr = mk_hrect(zeroes,zeroes);
  free_dyv(zeroes);
  return hr;
}

/* Makes a slightly larger hrect that encloses hr, but has "sensible" (i.e. round)
   numbers at the bootom and top of all its dimensions */
hrect *mk_hrect_sensible_limits(hrect *hr)
{
  hrect *newhr = mk_copy_hrect(hr);
  set_vector_limits_sensible(newhr->lo,newhr->hi);
  return newhr;
}

bool hrect_equal(hrect *hr1,hrect *hr2)
{
  return dyv_equal(hr1->lo,hr2->lo) &&
         dyv_equal(hr1->lo,hr2->lo);
}

/* Generates a vector randomly uniformly from hrect */
dyv *mk_sample_from_hrect(hrect *hr)
{
  int size = hrect_size(hr);
  int i;
  dyv *x = mk_dyv(size);
  for ( i = 0 ; i < size ; i++ )
    dyv_set(x,i,range_random(hrect_lo_ref(hr,i),hrect_hi_ref(hr,i)));
  return x;
}

/* Returns the area (volume?) of the hrect, i.e. the product of all
   its widths. */
double hrect_area(hrect *hr)
{
  /* go through each dimension, extract the values of the low and high
     points in that dimension, calculate d_low - d_high) and multiply
     that by the current area calculation */
  double area = 1.0;
  int i;
  for (i = 0; i < hrect_size(hr); i++)
  {
    double lo = hrect_lo_ref(hr,i);
    double hi = hrect_hi_ref(hr,i);
    area *= hi - lo;
  }

  return area;
}


/* 
   File:         hrect.h
   Author:       Andrew W. Moore
   Created:      Jan 13th, 1998
   Description:  The representation of motor control hrects.
*/

#ifndef HRECT_H
#define HRECT_H

#include "motutils.h"

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

typedef struct hrect
{
  /* Both the following vectors should be treated as private */
  dyv *lo;
  dyv *hi;
} hrect;

#define hrect_lo_ref(hr,i) (dyv_ref((hr)->lo,(i)))
#define hrect_hi_ref(hr,i) (dyv_ref((hr)->hi,(i)))
#define hrect_lo_set(hr,i,v) (dyv_set((hr)->lo,(i),(v)))
#define hrect_hi_set(hr,i,v) (dyv_set((hr)->hi,(i),(v)))

/* makes a hyper rectangle with lowest value in ith dimension
   dyv_ref(lo,i) and highest value dyv_ref(hi,i) */
hrect *mk_hrect(dyv *lo,dyv *hi);

void copy_hrect(hrect *src,hrect *dest);
hrect *mk_copy_hrect(hrect *hr);

/* This parses two strings (which must be space-separated ASCII
   representations of numbers). It turns the first into the low
   vector for the hyperrectnagle and the second into the high vector.
   Clearly both vectors must be the same length (i.e. have the same
   number of numbers). The number of numbers is the dimensionality of
   the hyper-rectangle */
hrect *mk_hrect_from_strings(char *lo_string,char *hi_string);


/* The dimensionality of the hyper-rectangle */
int hrect_size(const hrect *hr);


/* MAKES and returns the center of the hyper-rectangle. (This is
   of course just the mean of lo and hi) */
dyv *mk_hrect_middle(hrect *hr);


/* If x is inside the hyper-rectangle it is unchanged.
   If x is outside the hyper-rectangle it is changed to
   become the closest point to x that is in (or on the surface
   of) the hyper-rectangle */
void clip_to_inside_hrect(const hrect *hr,dyv *x);


/* Returns TRUE if and only if query_point is inside or on
   the edge of the hyper-rectangle */
bool is_in_hrect(const hrect *hyp, const dyv *query_point);


/* am_frees the hyper-rect. (And, in accordance with AM conventions,
   it frees all the substructures). */
void free_hrect(hrect *hyp);


/* Prints using the AM-FPRINT-CONVENTION described in motutils.h */
void fprint_hrect(FILE *s,char *s1,char *s2,hrect *hyp);

/* Prints using the AM-FPRINTF-CONVENTION */
void fprintf_hrect(FILE *s,char *m1,hrect *hyp,char *m2);

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
void hrect_set_ag_frame(hrect *hr,int i1,int i2);

/* As above, except shrinks the frame so that on all four sides
   of the ag_ window there's a margin of width 5% window width. */
void hrect_set_ag_frame_with_border(hrect *hr,int i1,int i2);

void ag_hrect(hrect *hr);

/* rows may be NULL denoting "use all rows" */
hrect *mk_hrect_bounding_dym_rows(dym *x,ivec *rows);
int hrect_widest_dim(hrect *hr);
double hrect_middle_ref(hrect *hr,int index);
double hrect_width_ref(hrect *hr,int index);

hrect *mk_hrect_bounding_dyv_array(const dyv_array *da);

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
                                         dyv *x_minlik,dyv *x_maxlik);

/* Returns a new hrect translated by subtracting the vector
   x from the lower limits and the higher limits of hr */
hrect *mk_hrect_subtract_dyv(hrect *hr,dyv *x);

hrect *mk_zero_hrect(int size);

/* Makes a slightly larger hrect that encloses hr, but has "sensible" (i.e. round)
   numbers at the bootom and top of all its dimensions */
hrect *mk_hrect_sensible_limits(hrect *hr);

bool hrect_equal(hrect *hr1,hrect *hr2);

/* Generates a vector randomly uniformly from hrect */
dyv *mk_sample_from_hrect(hrect *hr);

/* Returns the area (volume?) of the hrect, i.e. the product of all
   its widths. */
double hrect_area(hrect *hr);

#endif /* #ifndef HRECT_H */

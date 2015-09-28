/*
   File:        qopt.h
   Author:      Andrew W. Moore
   Created:     Tue Sep 24 14:42:15 EDT 1996
   Description: Box-constrained quadratic optimization

   Copyright (C) 1996, A. W. Moore
*/

#ifndef QOPT_H
#define QOPT_H

#include "./utils/amut.h"

#define ZD_MINIMUM 0
#define ZD_MAXIMUM 1
#define ZD_SADDLE  2

/*
Glossary:

   In the following functions, fixed_atts is an ivec saying which
   attribute values are frozen (fixed) at specific values.

  Virtual points: In the original space some attributes are fixed and
  others are variable. Virtual space is the smaller but equivalent problem
  in which all attributes are variable. Simple example

  Original problem: y = 3 x1^2 + 9 x2^2 - x1 x2 + 2 x1 + 3, with x1 constrained to be 7;

  Virtual problem: y = 3 * 7^2 + 9x2^2 - 7 x2 + 14 + 3,
               <=> y = 161 - 7 z + 9 z^2
          where (z) = mk_virtual_from_real(x1,x2)

  zdpoint: The value x at which a quadratic form has zero derivative in all directions.

  zd_type:

  As defined in qopt.h, a zd_type (pronounced "zero_derivative type")
   can be ZD_MINIMUM, ZD_MAXIMUM or ZD_SADDLE point. This is a descriptor
   for the type of quadratic form.
*/

/* Declared but not defined - sir 30/5/2000
ivec *mk_variable_atts(int num_atts,ivec *fixed_atts);
/ * dv += m * x (matrix times vector) * /
void dyv_add_mx(dyv *dv,dym *m,dyv *x);

/ * dv += m^T * x (matrix times vector) * /
void dyv_add_mtx(dyv *dv,dym *m,dyv *x);


/ * Given that some of the variables must be fixed in the
   quadratic form y = x^T A x + b x , make a smaller quadform
   that includes only the non-fixed attributes. Returns the result
   in *r_va, *r_vb *r_vc * /
void make_virtual_quad(dym *a,dyv *b,ivec *fixed_atts,dyv *fixed_vals,
                       dym **r_va,dyv **r_vb,double *r_vc);

/ * As defined in qopt.h, a zd_type (pronounced "zero_derivative type")
   can be ZD_MINIMUM, ZD_MAXIMUM or ZD_SADDLE point. This is a descriptor
   for the type of quadratic form * /
int zd_type_from_atpa(dym *atpa);
*/



/* Computes x^T A x + b^T x + c */
double eval_quadform(dym *a,dyv *b,double c,dyv *x);

/* Computes the point x for which the gradient of f(x) is zero in
   all directions.

    f(x) = x^T A x + b^T x + c

    so its a quadratic so (unless A is 0) such a point exists and is unique.
*/
dyv *mk_zdpoint_fast(dym *a,dyv *b);

/* virtual_x is a vector that only includes the non-fixed components
   of the original_x.
*/

/* Declared but not defined - sir 30/5/2000
dyv *mk_x_from_virtual_x(dyv *vx,ivec *fixed_atts,dyv *fixed_vals);


/ * Write the vector x = (vx,fx), i.e. broken into its fixed and
   variable components.

     Write y(vx) = (vx,fx)^T A (vx,fx) + b^T (vx,fx) + c

   This finds the (vx,fx) that makes dy/dvx(i) zero for all
   variable components.

   Note that the above description assumes the variable attributes
   appear before the fixed attributes. The truth is that they can
   appear in any order. The fixed attribute indexes are defined in fixed_atts
   and the fixed values in fixed_vals.
 * /
dyv *mk_zdpoint_with_fixed(dym *a,dyv *b,ivec *fixed_atts,dyv *fixed_vals,
			   double *r_zdpoint_val,int *r_zdpoint_type);
*/



#define LOLIM_ATT_TYPE         0
#define HILIM_ATT_TYPE         1
#define UNCONSTRAINED_ATT_TYPE 2
#define CONSTRAINED_ATT_TYPE   3


/* Write the vector x = (vx,fx), i.e. broken into its fixed and
   variable components.

     Write y(vx) = 0.5 * (vx,fx)^T A (vx,fx) + b^T (vx,fx) + c

   This finds the (vx,fx) that globalled maximizes y(vx,fx) for all
   variable components, constrained so that
      lo <= (vx,fx) <= hi
   where <= denotes vector dominance.

   Note that the above description assumes the variable attributes
   appear before the fixed attributes. The truth is that they can
   appear in any order. The fixed attribute indexes are defined in fixed_atts
   and the fixed values in fixed_vals.
*/
dyv *mk_constrained_quadopt(dym *a,dyv *b,double c,ivec *fixed_atts,dyv *fixed_vals,
			    dyv *lo,dyv *hi,double *r_opt_value);


typedef struct quadargs
{
  dym *a;
  dyv *b;
  ivec *fixed_atts;
  dyv *fixed_vals;
  dyv *lo;
  dyv *hi;
} quadargs;

bool is_dym_sym_pos_def_and_nonsingular(dym *a);
bool is_dym_sym_neg_def_and_nonsingular(dym *a);

#endif /* #ifndef QOPT_H */

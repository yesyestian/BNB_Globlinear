/*
   File:        qopt.c
   Author:      Andrew W. Moore
   Created:     Tue Sep 24 14:42:15 EDT 1996
   Description: Box-constrained quadratic optimization

   Copyright (C) 1996, A. W. Moore
*/

#include "qopt.h"

dym *mk_dym_from_file(FILE *s,char **r_errmess);

#define ZD_MINIMUM 0
#define ZD_MAXIMUM 1
#define ZD_SADDLE  2

/* 
Glossary:

   In the following functions, fixed_atts is an ivec saying which
   attribute values are frozen (fixed) at specific values.


  zdpoint: The value x at which a quadratic form has zero derivative in all directions.

  zd_type:

  As defined in qopt.h, a zd_type (pronounced "zero_derivative type")
   can be ZD_MINIMUM, ZD_MAXIMUM or ZD_SADDLE point. This is a descriptor
   for the type of quadratic form.

  Quadratic forms are used of the form

  y(x) = 0.5 x^T A x + b^T x + c where A is symmetric
     (but not necessarily positive definite).
*/

/* Computes 0.5 x^T A x + b^T x + c. Assume A is symmetric. */
double eval_quadform(dym *a,dyv *b,double c,dyv *x)
{
  double result = 0.5 * dym_xt_a_x_value(x,a) + 
                  dyv_scalar_product(x,b) + c;
  return(result);
}

/* Computes the point x for which the gradient of f(x) is zero in
   all directions.

    f(x) = 0.5 * x^T A x + b^T x + c

   If this point is a maximum, returns x.
   
   If it's not a maximum, doesn't bother to return x.

  Assume A is symmetric.
*/
dyv *mk_zdpoint_fast(dym *a,dyv *b)
{
  dyv *zdpoint = NULL;

  if ( !is_dym_symmetric(a) )
    my_error("mk_zdpoint_fast: precondition violated");
  else
  {
    int size = dyv_size(b);
    dym *minus_a = mk_dym_scalar_mult(a,-1.0);
    dym *lot = mk_dym(size,size);
    bool minus_a_is_spd = attempt_cholesky_decomp(minus_a,lot);

    if ( minus_a_is_spd )
    {
      /* We know that lot lot^T = a , and lot is lower triangular */
      zdpoint = mk_solve_lot_equation(lot,b);
      /* If lot is singular, zdpoint is now NULL.
         If lot nonsingular, we now have lot^T lot zdpoint = b */
    }
  
    if ( Verbosity > 30.0 )
    {
      fprintf_dym(stdout,"a",a,"\n");
      fprintf_dyv(stdout,"b",b,"\n");
      if ( zdpoint == NULL )
        printf("The quadratic y(x) = x^t a x + b^t x has no maximum\n");
      else
      {
        dym *lott = mk_dym_transpose(lot);
        dym *lot_lott = mk_dym_mult(lot,lott);
        dyv *minus_a_zdpoint = mk_dym_times_dyv(minus_a,zdpoint);
        fprintf_dym(stdout,"minus_a",minus_a,"\n");
        fprintf_dym(stdout,"lott",lott,"\n");
        fprintf_dym(stdout,"lot_lott",lot_lott,"\n");
        fprintf_dyv(stdout,"minus_a_zdpoint",minus_a_zdpoint,"\n");
        fprintf_dym(stdout,"lot",lot,"\n");
        fprintf_dyv(stdout,"max of y(x) = x^t a x + b^t x",zdpoint,"\n");
        free_dym(lott);
        free_dym(lot_lott);
        free_dyv(minus_a_zdpoint);
      }
      wait_for_key();
    }
  
    free_dym(minus_a);
    free_dym(lot);
  }
  
  return(zdpoint);
}

/*
"subject to the constraint xj=q,

 x1*,x2*, ... xj*, ... xN* is the maximizer of
 the following expression, and maximum value is v,

       [x1]T [a11 .... a1j .... a1n] [x1]   [b1]T [x1]
       [x2]  [a21 .... a2j .... a2n] [x2]   [b2]  [x2]
       [ .]  [.         .        . ] [ .]   [ .]  [ .]
       [ .]  [.         .        . ] [ .]   [ .]  [ .]
 0.5 * [xj]  [aj1 .... ajj .... ajn] [xj] + [bj]  [xj] + c
       [ .]  [.         .        . ] [ .]   [ .]  [ .]
       [ .]  [.         .        . ] [ .]   [ .]  [ .]
       [xn]  [an1 .... anj .... ann] [xn]   [bn]  [xn]
"

The above is an equivalent problem statement to

"maximize the following expression, and find its
 maximum value v..."
       [x1]T   [a11 .... a1,j-1   a1,j+1 .... a1n   ] [x1  ] 
       [x2]    [a21 .... a2,j-1   a2,j+1 .... a2n   ] [x2  ] 
       [ .]    [.         .      .              .   ] [ .  ]   
       [ .]    [.         .      .              .   ] [ .  ]   
 0.5 * [xj-1]  [aj-1,1 ..aj-1,j-1 aj-1,j+1 .. aj-1,n] [xj-1]
       [xj+1]  [aj+1,1 ..aj+1,j-1 aj+1,j+1 .. aj+1,n] [xj+1]
       [ .]    [.         .      .              .   ] [ .  ]   
       [ .]    [.         .      .              .   ] [ .  ]   
       [xn]    [an1 .... an,j-1   an,j+1 .... ann   ] [xn  ] 
 
      +

[b1   + q a1j   ]T [x1  ]
[b2   + q a2j   ]  [x2  ]
[ .             ]  [ .  ]
[ .             ]  [ .  ]
[bj-1 + q aj-1,j]  [xj-1] + c + 0.5 * q^2 ajj + qbj
[ .             ]  [ .  ]
[ .             ]  [ .  ]
[bn   + q  anj  ]  [xn  ]

  Assume A is symmetric.
*/

/* The following functions operationalize the above 
   equivalence, turning an N-dimensional 
   quadratic optimization problem over an N-dimensional
   hyper-rectangle with one equality constraint
   into an N-1 dimensional problem over an N-1 dimensional
   box.
*/
void quadopt_regular_to_eqv(int j,double q,
                            dym *a,dyv *b,double c,
                            dyv *lo,dyv *hi,ivec *att_types,
                            dym **eqv_a,dyv **eqv_b,double *eqv_c,
                            dyv **eqv_lo,dyv **eqv_hi,ivec **eqv_att_types)
{
  int size = dyv_size(b);
  int eqv_size = size-1;
  int eqv_i;

  *eqv_a = mk_dym(eqv_size,eqv_size);
  *eqv_b = mk_dyv(eqv_size);
  *eqv_lo = mk_dyv(eqv_size);
  *eqv_hi = mk_dyv(eqv_size);
  *eqv_att_types = mk_ivec(eqv_size);

  for ( eqv_i = 0 ; eqv_i < eqv_size ; eqv_i++ )
  {
    int i = (eqv_i < j) ? eqv_i : eqv_i + 1;
    double eqv_b_eqv_i = dyv_ref(b,i) + q * dym_ref(a,i,j);
    int eqv_k;

    for ( eqv_k = 0 ; eqv_k < eqv_size ; eqv_k++ )
    {
      int k = (eqv_k < j) ? eqv_k : eqv_k + 1;
      dym_set(*eqv_a,eqv_i,eqv_k,dym_ref(a,i,k));
    }
    dyv_set(*eqv_b,eqv_i,eqv_b_eqv_i);
    dyv_set(*eqv_lo,eqv_i,dyv_ref(lo,i));
    dyv_set(*eqv_hi,eqv_i,dyv_ref(hi,i));
    ivec_set(*eqv_att_types,eqv_i,ivec_ref(att_types,i));
  }

  /* eqv_c = c + 0.5 * q^2 ajj + qbj */

  *eqv_c = c + 0.5 * q * q * dym_ref(a,j,j) + q * dyv_ref(b,j);
}

void xopt_regular_to_eqv(int j,dyv *xopt,dyv **eqv_xopt)
{
  int size = dyv_size(xopt);
  int eqv_size = size-1;
  int eqv_i;

  *eqv_xopt = mk_dyv(eqv_size);

  for ( eqv_i = 0 ; eqv_i < eqv_size ; eqv_i++ )
  {
    int i = (eqv_i < j) ? eqv_i : eqv_i + 1;
    dyv_set(*eqv_xopt,eqv_i,dyv_ref(xopt,i));
  }
}

void xopt_eqv_to_regular(int j,double q,dyv *eqv_xopt,dyv **xopt)
{
  int eqv_size = dyv_size(eqv_xopt);
  int size = eqv_size + 1;
  int eqv_i;

  *xopt = mk_dyv(size);

  for ( eqv_i = 0 ; eqv_i < eqv_size ; eqv_i++ )
  {
    int i = (eqv_i < j) ? eqv_i : eqv_i + 1;
    dyv_set(*xopt,i,dyv_ref(eqv_xopt,eqv_i));
  }

  dyv_set(*xopt,j,q);
}



/* Assume y(x) = 0.5 * x^T a x + b^T x. Assume A is symmetric.

   What is dy/x_k (ie the partial derivative of y with respect to k'th
   component of x?)
   
   It's a linear function of x. Write dy/dx_k = p + q^T x
   
   Then p = b[k]
   
   q[i] = a[k,i]
   
     (Proof: if i != k, then the two relevant terms in the above definition
       of y(x) are 
y(x) =...+ 0.5 a[i,k] * x[i] * x[k] +...+ 0.5 a[k,i] * x[k] * x[i] +..
the partial derivative wrt x[k] of which is 
   0.5 * (a[i,k] + a[k,i]) * x[i]
   = a[i,k] * x[i] because a[ik]==a[ki]
       
       and if i == k, the relevant term is 0.5 * a[k,k] * x[k]^2
       the partial derivative of which is a[k,k] * x[k]

   This function Makes and returns the vector q, and it returns p in the 
   third argument.
*/
dyv *mk_derivative_function(dym *a,dyv *b,int k,double *r_p)
{
  /* int num_atts = dyv_size(b); */
  dyv *q = mk_dyv_from_dym_row(a,k);
  
  *r_p = dyv_ref(b,k);

  return(q);
}

dyv *mk_rectangle_center(dyv *lo,dyv *hi)
{
  dyv *x = mk_dyv_plus(lo,hi);
  dyv_scalar_mult(x,0.5,x);
  return(x);
}

/* Assume y(x) = 0.5 * x^T a x + b^T x. Assume A symmetric.

    Assume that x is allowed to vary within the hyperrectangle defined
    by lo and hi.

    Then what is the minimum value of dy/dx_k subject to contraints?
     and what is the maximum value of dy/dx_k subject to contraints?
*/
void compute_derivative_limits(dym *a,dyv *b,dyv *lo,dyv *hi,int k,
                               double *r_deriv_lo,double *r_deriv_hi)
{
  double p;
  dyv *q = mk_derivative_function(a,b,k,&p);
    /* dy/dx[k] = p + q^T x */
  int att_num;
  int num_atts = dyv_size(b);

  *r_deriv_lo = p;
  *r_deriv_hi = p;

  for ( att_num = 0 ; att_num < num_atts ; att_num++ )
  {
    double lo_val = dyv_ref(lo,att_num);
    double hi_val = dyv_ref(hi,att_num);
    double lo_contribution;
    double hi_contribution;
    double qatt = dyv_ref(q,att_num);

    if ( qatt < 0.0 )
    {
      lo_contribution = hi_val * qatt;
      hi_contribution = lo_val * qatt;
    }
    else
    {
      lo_contribution = lo_val * qatt;
      hi_contribution = hi_val * qatt;
    }

    *r_deriv_lo += lo_contribution;
    *r_deriv_hi += hi_contribution;
  }

  free_dyv(q);
}

#define LOLIM_ATT_TYPE         0
#define HILIM_ATT_TYPE         1
#define UNCONSTRAINED_ATT_TYPE 2
#define CONSTRAINED_ATT_TYPE   3

/* See compute_constrained_quadopt comments for meaning of above */

/* It may be the case that we can deduce that one of the CONSTRAINED_ATT_TYPE
   variables will be driven to an upper or lower limit. This'll happen if the 
   partial derivative wrt to some attribute is seen to be non-negative or non-positive
   throughout the hyperrectangle.

   This function tries to find such a variable.

   If there is no such variable it returns -1, and sets *r_known_att_type to -1.

   If there is such a variable it returns the att_num (ie index) of that variable
   and returns either LOLIM_ATT_TYPE or HILIM_ATT_TYPE in *r_known_att_type
   depending on whether the variable is known driven to its low or high limit
   respectively.
*/
int find_known_lo_or_hi_att_num(dym *a,dyv *b,dyv *lo,dyv *hi,ivec *att_types,
                                int *r_known_att_type)
{
  int num_atts = dyv_size(b);
  int att_num;
  int result_att_num = -1;
  *r_known_att_type = -1;
  for ( att_num = 0 ; result_att_num < 0 && att_num < num_atts ; att_num++ )
  {
    int att_type = ivec_ref(att_types,att_num);
    if ( att_type == CONSTRAINED_ATT_TYPE )
    {
      double deriv_lo,deriv_hi;
      compute_derivative_limits(a,b,lo,hi,att_num,
                                &deriv_lo,&deriv_hi);
      if ( deriv_lo >= 0.0 || deriv_hi <= 0.0 )
      {
        result_att_num = att_num;
        *r_known_att_type = (deriv_lo >= 0.0) ? HILIM_ATT_TYPE : LOLIM_ATT_TYPE;
      }
    }
  }

  return(result_att_num);
}

/* After calling
    *r_lo_deriv[k] = min value of partial derivative of y wrt x_k over the
                     whole rectangle.

    *r_hi_deriv[k] = corresponding hi value.
*/
void compute_var_att_deriv_limits(dym *a,dyv *b,dyv *lo,dyv *hi,
                                  dyv **r_lo_deriv,dyv **r_hi_deriv)
{
  int num_atts = dyv_size(b);
  int att_num;

  *r_lo_deriv = mk_dyv(num_atts);
  *r_hi_deriv = mk_dyv(num_atts);

  for ( att_num = 0 ; att_num < num_atts ; att_num++ )
  {
    double deriv_lo,deriv_hi;

    compute_derivative_limits(a,b,lo,hi,att_num,
                              &deriv_lo,&deriv_hi);
    
    dyv_set(*r_lo_deriv,att_num,deriv_lo);
    dyv_set(*r_hi_deriv,att_num,deriv_hi);
  }
}

/* Computes the max possible value of y(x) = xt a x + bt x subject to regular
   constraints, and using the knowlegde that
           var_att_min_deriv[k] is the minimum value of dy/dx_k 
           var_att_max_deriv[k] is the maximum value of dy/dx_k 

  Theorem:

           max             f(x)     is less than or equal to....
  x in [lo,hi] hyperect,

         f(xc) + SUM_k max( -hw_k dmin_k , hw_k + dmax_k ) 

     xc[k] = (hi[k] + lo[k]) / 2 

     hw_k = (hi[k] - lo[k) / 2
     dmin_k = var_att_min_deriv[k]
     dmax_k = var_att_max_deriv[k]
*/
double compute_max_possible(dym *a,dyv *b,double c,dyv *lo,dyv *hi,
                            dyv *var_att_min_deriv,
                            dyv *var_att_max_deriv)
{
  dyv *xc = mk_rectangle_center(lo,hi);
  double f_xc = eval_quadform(a,b,c,xc);
  double result = f_xc;
  int k;

  for ( k = 0 ; k < dyv_size(b) ; k++ )
  {
    double hw_k = (dyv_ref(hi,k) - dyv_ref(lo,k)) / 2.0;
    double dmin_k = dyv_ref(var_att_min_deriv,k);
    double dmax_k = dyv_ref(var_att_max_deriv,k);
    double max_contrib = real_max(-hw_k * dmin_k,hw_k * dmax_k);

    result += max_contrib;
  }
  free_dyv(xc);
  return(result);
}

bool dyv_strictly_outside_bounds(dyv *x,dyv *lo,dyv *hi)
{
  bool result = !dyv_weakly_dominates(x,lo) ||
                !dyv_weakly_dominates(hi,x);
  return(result);
}

int Num_maxes = 0; /* Just used for profiling */
int Num_nodes = 0;

/* Find xopt = the x in [lo,hi] that maximizes 

     y(x) = 0.5 a x^2 + b x + c.

   If *r_xopt is NULL or y(xopt) > *r_opt_value,
   set *r_xopt to be the singleton vector (x) and *r_opt_value := y(xopt).

    Else leave everything unchanged
*/
void constrained_quadopt_1d(double a,double b,double c,double lo,double hi,
                            dyv **r_xopt,double *r_opt_value)
{
  bool has_max = a < 0.0;
  bool found_max = FALSE;
  double xopt = 0.0;
  double value = 0.0;

  if ( has_max )
  {
    xopt = -b / a;
    found_max = xopt >= lo && xopt <= hi;
    value = 0.5 * a * xopt * xopt + b * xopt + c;
  }

  if ( !found_max )
  {
    double lo_value = 0.5 * a * lo * lo + b * lo + c;
    double hi_value = 0.5 * a * hi * hi + b * hi + c;
    if ( lo_value >= hi_value )
    {
      xopt = lo;
      value = lo_value;
    }
    else
    {
      xopt = hi;
      value = hi_value;
    }
  }

  if ( Verbosity >= 20.0 )
  {
    printf("1_d opt a=%g b=%g c=%g lo=%g hi=%g xopt=%g value=%g\n",
           a,b,c,lo,hi,xopt,value);
    if ( Verbosity > 20.0) wait_for_key();
  }

  if ( *r_xopt == NULL || value > *r_opt_value )
  {
    if ( *r_xopt == NULL ) 
      *r_xopt = mk_dyv_1(xopt);
    else
      dyv_set(*r_xopt,0,xopt);
    
    *r_opt_value = value;
  }
}

/*
   Let y(x) = 0.5 * x^T a x + b^t x + c.

   Let xopt be the value of x that maximizes y(x) subject to the following
   constraints.
  
   If the i'th attribute has att_type UNCONSTRAINED_ATT_TYPE, then x[i] may
    take any value between -infinity to +infinity

   If the i'th attribute has att_type CONSTRAINED_ATT_TYPE, then x[i] must 
   lie in the closed interval [ lo[i],hi[i] ]: lo[i] <= x[i] <= hi[i].

   If xopt lies in the closed hyperrectangle defined by lo and hi, then this
   functions MAKES  xopt. Otherwise, if xopt lies strictly outside this
   region, the function sets xopt == NULL.

  PRE: att_types[j] == CONSTRAINED_ATT_TYPE or UNCONSTRAINED_ATT_TYPE

  PRE: *r_xopt must enter with *r_xopt = NULL or *r_xopt = a legally allocated dyv.

  PRE: if *r_xopt is non-null on entry, *r_opt_value = y(*r_xopt) on entry.

  POST:
  If xopt gets set to NULL by the function (ie if maximum solution is strictly
  outside limits or is infinite), *r_xopt and *r_opt_value are left unchanged.
  
  If xopt gets set non-NULL and (*r_xopt==NULL or y(xopt) > y(*r_xopt)), then
  (if necessary) *r_xopt is freed and replaced by xopt. *r_opt_value is set
  to y(xopt).

  Else, if xopt is non-null, but *r_xopt is non-null and y(xopt) <= *r_opt_value,
  then *r_xopt and r_opt_value are left unchanged.

  Consequently:
    On exit, *r_xopt == NULL, or *r_opt_value == y(*r_xopt)

  do_cutoff_check: You can check for cutoffs by seeing if a variable has 
                   a strictly +ve or -ve derivative over whole region, and by
                   checking if we're negative definite. This costs though, and
                   sometimes you know in advance of calling that no cutoff will
                   be found (when? when you've changed an attribute from
                   CONSTRAINED to UNCONSTRAINED). The caller can signal whether its
                   worth doing a cutoff check by setting do_cutoff_check.

                   Top-level caller should set it to TRUE.
*/
void constrained_quadopt(dym *a,dyv *b,double c,dyv *lo,dyv *hi,ivec *att_types,
                         bool do_cutoff_check,
                         dyv **r_xopt,double *r_opt_value)
{
  /* Declare mutually recursive helper function ... */
  void constrained_quadopt_helper(dym *a,dyv *b,double c,dyv *lo,dyv *hi,ivec *att_types,
                                dyv **r_xopt,double *r_opt_value,
                                int k,int k_att_type);
  int num_atts = dyv_size(b);

  if ( num_atts == 1 )
    constrained_quadopt_1d(dym_ref(a,0,0),dyv_ref(b,0),c,dyv_ref(lo,0),dyv_ref(hi,0),
                           r_xopt,r_opt_value);
  else
  {
    int constrained_index = find_index_in_ivec(att_types,CONSTRAINED_ATT_TYPE);
  
    if ( Verbosity > 10.0 ) 
    {
      fprintf_ivec(stdout,"Enter cqh att_types",att_types,"\n");
      fprintf_dym(stdout,"a",a,"\n");
      fprintf_dyv(stdout,"b",b,"\n");
      printf("c = %g\n",c);
      if ( *r_xopt == NULL )
        printf("*r_xopt on entry = NULL\n");
      else
      {
        fprintf_dyv(stdout,"*r_xopt on entry",*r_xopt,"\n");
        printf("*r_opt_val on entry = %g\n",*r_opt_value);
      }
      if ( Verbosity > 20.0 ) wait_for_key();
    }
  
    Num_nodes += 1;
  
    if ( constrained_index < 0 ) /* There is no k s.t. att_types[k] == CONSTRAINED_ATT_TYPE */
    {
      dyv *xopt;
  
      if ( num_atts == 0 )
        xopt = mk_dyv(0);
      else
      {
        Num_maxes += 1;
        xopt = mk_zdpoint_fast(a,b);
      }
  
      if ( xopt == NULL )
      {
        if ( Verbosity > 10.0 ) printf("This quadratic has no maximum\n");
      }
      else
      {
        bool weakly_inside_bounds = !dyv_strictly_outside_bounds(xopt,lo,hi);
  
        if ( Verbosity > 10.0 )
        {
          printf("  ..when all else unconstrained, the surface has\n");
          printf("    a maximum at the following location (x)..\n");
          fprintf_dyv(stdout,"x",xopt,"\n");
          printf("We are%s weakly inside the bounds\n",(weakly_inside_bounds)?"":"n't");
        }
  
        if ( weakly_inside_bounds )
        {
          double value = eval_quadform(a,b,c,xopt);
          bool better_than_on_entry = ( *r_xopt == NULL || value > *r_opt_value );
  
          if ( Verbosity > 10.0 )
          {
            printf("This %s the best legal value I've seen so far.\n",
                   (better_than_on_entry) ? "is" : "isn't");
          }
  
          if ( better_than_on_entry )
          {
            if ( *r_xopt != NULL ) free_dyv(*r_xopt);
            *r_xopt = mk_copy_dyv(xopt);
            *r_opt_value = value;
          }
        }
  
        free_dyv(xopt);
      }
  
      if ( Verbosity > 20.0 ) wait_for_key();
    }
    else  /* att_types[constrained_index] == CONSTRAINED_ATT_TYPE */
    {
      bool do_all_three_children = TRUE;
      bool first_lim_att_type = LOLIM_ATT_TYPE;  /* If you do both LO and HI child, which first? */
      bool second_lim_att_type = HILIM_ATT_TYPE;
      
      if ( do_cutoff_check )
      {
        dyv *zdpoint = mk_zdpoint_fast(a,b);
       
        if ( zdpoint != NULL )
        {
          double value = eval_quadform(a,b,c,zdpoint);
          if ( !dyv_strictly_outside_bounds(zdpoint,lo,hi) )
          {
            if ( *r_xopt == NULL || *r_opt_value < value )
            {
              if ( *r_xopt != NULL ) free_dyv(*r_xopt);
              *r_xopt = mk_copy_dyv(zdpoint);
              *r_opt_value = value;
            }
            /* No point in exploring any children */
            do_all_three_children = FALSE;
          }
          else if ( *r_xopt != NULL && *r_opt_value >= value )
          {
            do_all_three_children = FALSE;
              /* No point in searching: none of our children can possibly do any
                 better than something we've already seen */
          }
          free_dyv(zdpoint);
        }
        
        if ( do_all_three_children )
        {
  
          /* I will now search to see if any of the dimensions give us a
             constrained cutoff */
          int known_att_type;
          int known_lo_or_hi_att_num = 
            find_known_lo_or_hi_att_num(a,b,lo,hi,att_types,&known_att_type);
      
          if ( known_lo_or_hi_att_num >= 0 )
          {
            constrained_quadopt_helper(a,b,c,lo,hi,att_types,
                                       r_xopt,r_opt_value,
                                       known_lo_or_hi_att_num,known_att_type);
            do_all_three_children = FALSE;
          }
          else
          {
            dyv *var_att_min_deriv;
            dyv *var_att_max_deriv;
            double max_possible_value;
  
            compute_var_att_deriv_limits(a,b,lo,hi,
                                         &var_att_min_deriv,
                                         &var_att_max_deriv);
  
          /* var_att_min_deriv[k] is the minimum value of dy/dx_k 
             over the whole rectangle. */
  
            max_possible_value = compute_max_possible(a,b,c,lo,hi,
                                                      var_att_min_deriv,
                                                      var_att_max_deriv);
  
            if ( *r_xopt != NULL && max_possible_value <= *r_opt_value )
              do_all_three_children = FALSE;
  
            free_dyv(var_att_min_deriv);
            free_dyv(var_att_max_deriv);
          }
        }
      }
  
      if ( do_all_three_children )
      {
  
        /* It is not worth doing the unconstrained child if the current constrained_index
           is the ONLY CONSTRAINED_ATT_TYPE. Why not? Because then the subnode will
           be a leaf node, which will repeat the calculation done in the cutoff
           check. So we will make sure we don't call the helper with UNCONSTRAINED_ATT_TYPE
           if constrained_index is the only att_num whose att_type is CONSTRAINED_ATT_TYPE */
  
        if ( num_of_given_value(att_types,CONSTRAINED_ATT_TYPE) > 1 )
          constrained_quadopt_helper(a,b,c,lo,hi,att_types,r_xopt,r_opt_value,
                                     constrained_index,UNCONSTRAINED_ATT_TYPE);
  
        constrained_quadopt_helper(a,b,c,lo,hi,att_types,r_xopt,r_opt_value,
                                   constrained_index,first_lim_att_type);
        constrained_quadopt_helper(a,b,c,lo,hi,att_types,r_xopt,r_opt_value,
                                   constrained_index,second_lim_att_type);
      }
    }
  
    if ( Verbosity > 10.0 )
    {
      fprintf_ivec(stdout,"Exit cqh att_types",att_types,"\n");
      if (  *r_xopt == NULL )
        printf("xopt = NULL\n");
      else
      {
        fprintf_oneline_dyv(stdout,"xopt",*r_xopt,"\n");
        printf("opt value = %g\n",*r_opt_value);
      }
      if ( Verbosity > 20.0 ) wait_for_key();
    }
  }
}
  
/*
   This function does a little preproccessing and then calls constrained_quadopt
   followed by tiny postprocessing.

   If k_att_type == LOLIM_ATT_TYPE or HILIM_ATT_TYPE, this function solves
   an equivalent, one-fewer-dimension problem, and then (if it has a better
   result than (*r_xopt and *r_opt_value on entry) it replaces (*r_xopy and
   *r_opt_value) accordingly.
 
   If k_att_type == UNCONSTRAINED_ATT_TYPE, calls constrained_quadopt
   with att_types[k] == UNCONSTRAINED_ATT_TYPE.

  PRE:
    all preconditions for constrained_quadopt

    att_types[k] = CONSTRAINED_ATT_TYPE;

    k_att_type is one of HILIM_ATT_TYPE, LOLIM_ATT_TYPE or UNCONSTRAINED_ATT_TYPE
*/
void constrained_quadopt_helper(dym *a,dyv *b,double c,dyv *lo,dyv *hi,ivec *att_types,
                                dyv **r_xopt,double *r_opt_value,
                                int k,int k_att_type)
{
  if ( Verbosity > 10.0 )
  {
    printf("I will now set attribute %d to %s value",
           k,(k_att_type==UNCONSTRAINED_ATT_TYPE) ? "an unconstrained" : 
             (k_att_type==LOLIM_ATT_TYPE) ? "its lowest" : "its highest");

    if ( k_att_type == LOLIM_ATT_TYPE )
      printf("of %g\n",dyv_ref(lo,k));
    else if ( k_att_type == HILIM_ATT_TYPE )
      printf("of %g\n",dyv_ref(hi,k));

    printf("\n");
    if ( Verbosity > 20.0 )
      wait_for_key();
  }

  if ( ivec_ref(att_types,k) != CONSTRAINED_ATT_TYPE ||
       !(k_att_type == LOLIM_ATT_TYPE || k_att_type == HILIM_ATT_TYPE || k_att_type == UNCONSTRAINED_ATT_TYPE )
     )
  {
    my_error("constrained_quadopt_helper: precondition violation");
  }

  if ( k_att_type == UNCONSTRAINED_ATT_TYPE )
  {
    ivec_set(att_types,k,UNCONSTRAINED_ATT_TYPE);
    constrained_quadopt(a,b,c,lo,hi,att_types,
                        FALSE, /* do_cutoff_check = no thanks */
                        r_xopt,r_opt_value);
    ivec_set(att_types,k,CONSTRAINED_ATT_TYPE);
  }
  else
  {
    dyv *lim = (k_att_type == LOLIM_ATT_TYPE) ? lo : hi;
    double q = dyv_ref(lim,k);
    dym *eqv_a;
    dyv *eqv_b,*eqv_lo,*eqv_hi,*eqv_xopt;
    double eqv_c;
    ivec *eqv_att_types;
    double value = *r_opt_value;

    quadopt_regular_to_eqv(k,q,
                           a,b,c,lo,hi,att_types,
                           &eqv_a,&eqv_b,&eqv_c,&eqv_lo,&eqv_hi,&eqv_att_types);
    if ( *r_xopt != NULL )
      xopt_regular_to_eqv(k,*r_xopt,&eqv_xopt);
    else
      eqv_xopt = NULL;

    constrained_quadopt(eqv_a,eqv_b,eqv_c,eqv_lo,eqv_hi,eqv_att_types,
                        TRUE, /* do_cutoff_check = yessirree */
                        &eqv_xopt,&value);

    /* Now we must be careful (Andrew blew 90 mins on a bug here).
       We must only update *r_xopt if value is better than *r_opt_value.
       (it is NOT true that it is harmless to rebuild *r_xopt in any case) */

    if ( eqv_xopt != NULL && (*r_xopt == NULL || value > *r_opt_value) )
    {
      if ( *r_xopt != NULL ) free_dyv(*r_xopt);
      xopt_eqv_to_regular(k,q,eqv_xopt,r_xopt);
      *r_opt_value = value;
    }

    free_dym(eqv_a);
    free_dyv(eqv_b);
    free_dyv(eqv_lo);
    free_dyv(eqv_hi);
    if ( eqv_xopt != NULL ) free_dyv(eqv_xopt);
    free_ivec(eqv_att_types);
  }
}

/* Write the vector x = (vx,fx), i.e. broken into its fixed and
   variable components. 

     Write y(vx) = 0.5 (vx,fx)^T A (vx,fx) + b^T (vx,fx) + c

   This finds the (vx,fx) that globally maximizes y(vx,fx) for all
   variable components, constrained so that 
      lo <= (vx,fx) <= hi
   where <= denotes vector dominance.

   Note that the above description assumes the variable attributes
   appear before the fixed attributes. The truth is that they can 
   appear in any order. The fixed attribute indexes are defined in fixed_atts
   and the fixed values in fixed_vals.
*/
dyv *mk_constrained_quadopt(dym *a,dyv *b,double c,ivec *fixed_atts,dyv *fixed_vals,
			    dyv *lo,dyv *hi,double *r_opt_value)
{
  dyv *temp_lo = mk_copy_dyv(lo);
  dyv *temp_hi = mk_copy_dyv(hi);
  dyv *xopt = NULL;
  int i;
  ivec *att_types = mk_constant_ivec(dyv_size(lo),CONSTRAINED_ATT_TYPE);

  *r_opt_value = -7777.7777;  /* debugging helper; this value is undefined */
  
  for ( i = 0 ; i < ivec_size(fixed_atts) ; i++ )
  {
    int att_num = ivec_ref(fixed_atts,i);
    dyv_set(temp_lo,att_num,dyv_ref(fixed_vals,i));
    dyv_set(temp_hi,att_num,dyv_ref(fixed_vals,i));
  }

  constrained_quadopt(a,b,c,temp_lo,temp_hi,att_types,
                      TRUE,  /* do_cutoff_check */
                      &xopt,r_opt_value);

  free_dyv(temp_lo);
  free_dyv(temp_hi);
  free_ivec(att_types);

  if ( xopt == NULL )
    my_error("Assert: xopt non null");

  return(xopt);
}
  
void save_quadargs_to_file(FILE *s, quadargs *qa)
{
  fprintf(s, "Quadargs\n");
  fprintf(s,"a\n");
  save_dym_to_file(s,qa->a);
  fprintf(s,"b\n");
  save_dyv_to_file(s,qa->b);
  fprintf(s,"fixed_atts\n");
  save_ivec_to_file(s,qa->fixed_atts);
  fprintf(s,"fixed_vals\n");
  save_dyv_to_file(s,qa->fixed_vals);
  fprintf(s,"lo\n");
  save_dyv_to_file(s,qa->lo);
  fprintf(s,"hi\n");
  save_dyv_to_file(s,qa->hi);
}

void free_quadargs(quadargs *qa)
{
  if ( qa->a != NULL ) free_dym(qa->a);
  if ( qa->b != NULL ) free_dyv(qa->b);
  if ( qa->fixed_atts != NULL ) free_ivec(qa->fixed_atts);
  if ( qa->fixed_vals != NULL ) free_dyv(qa->fixed_vals);
  if ( qa->lo != NULL ) free_dyv(qa->lo);
  if ( qa->hi != NULL ) free_dyv(qa->hi);
  AM_FREE(qa,quadargs);
}

quadargs *mk_empty_quadargs()
{
  quadargs *qa = AM_MALLOC(quadargs);
  qa->a = NULL;
  qa->b = NULL;
  qa->fixed_atts = NULL;
  qa->fixed_vals = NULL;
  qa->lo = NULL;
  qa->hi = NULL;
  return(qa);
}

/* If *r_errmess non-null, returns and does nothing.
   Else it checks that the next line in the file matches the string.
   If it doesn't, makes a complaining expo and stores it in *r_err_expo */
void check_next_line(FILE *s,char *string,char **r_errmess)
{
  if ( *r_errmess == NULL )
  {
    string_array *line = mk_string_array_from_line(s);
    if ( line == NULL )
      *r_errmess = mk_copy_string("File unexpectedly ended/empty");
    else if ( string_array_size(line) == 0 )
    {
      *r_errmess = mk_printf("Problem reading from file. "
                             "Expected the token %s, but instead saw an empty line.",
                             string);
    }
    else if ( !caseless_eq_string(string_array_ref(line,0),string) )
    {
      *r_errmess = mk_printf("Problem reading from file. "
                             "Expected the token %s, but instead saw the token %s.",
                             string,string_array_ref(line,0));
    }
    if ( line != NULL ) free_string_array(line);
  }
}

/* Returns NULL if no problem */
char *mk_check_quadargs_errmess(quadargs *qa)
{
  int size = dym_rows(qa->a);
  bool size_bad = dym_cols(qa->a) != size ||
                  dyv_size(qa->b) != size ||
                  dyv_size(qa->lo) != size ||
                  dyv_size(qa->hi) != size;
  char *errmess = NULL;

  if ( size_bad )
    errmess = mk_copy_string("Incompatible a rows / a cols / b size / lo size / hi size");
  if ( errmess == NULL && ivec_size(qa->fixed_atts) != dyv_size(qa->fixed_vals) )
    errmess = mk_copy_string("Incompatible fixed_atts / fixed vals sizes");
  if ( errmess == NULL && ivec_size(qa->fixed_atts) > 0 && ivec_max(qa->fixed_atts) >= size )
    errmess = mk_copy_string("One of the fixed_atts att numbers is too large");
  if ( errmess == NULL && ivec_size(qa->fixed_atts) > 0 && ivec_min(qa->fixed_atts) < 0 )
    errmess = mk_copy_string("One of the fixed_atts att numbers is negative");

  return(errmess);
}

quadargs *mk_quadargs_from_file(FILE *s, char **r_errmess)
{
  quadargs *qa = mk_empty_quadargs();

  *r_errmess = NULL;

  check_next_line(s,"Quadargs",r_errmess);
  check_next_line(s,"a",r_errmess);

  if ( *r_errmess == NULL ) qa->a = mk_dym_from_file(s,r_errmess);
  check_next_line(s,"b",r_errmess);
  if ( *r_errmess == NULL ) qa->b = mk_dyv_from_file(s,r_errmess);
  check_next_line(s,"fixed_atts",r_errmess);
  if ( *r_errmess == NULL ) qa->fixed_atts = mk_ivec_from_file(s,r_errmess);
  check_next_line(s,"fixed_vals",r_errmess);
  if ( *r_errmess == NULL ) qa->fixed_vals = mk_dyv_from_file(s,r_errmess);
  check_next_line(s,"lo",r_errmess);
  if ( *r_errmess == NULL ) qa->lo = mk_dyv_from_file(s,r_errmess);
  check_next_line(s,"hi",r_errmess);
  if ( *r_errmess == NULL ) qa->hi = mk_dyv_from_file(s,r_errmess);

  if ( *r_errmess == NULL )
    *r_errmess = mk_check_quadargs_errmess(qa);

  if ( *r_errmess != NULL )
  {
    free_quadargs(qa);
    qa = NULL;
  }

  return(qa);
}

dyv *mk_random_dyv(int size)
{
  dyv *x = mk_dyv(size);
  int i;
  for ( i = 0 ; i < size ; i++ )
    dyv_set(x,i,range_random(-1.0,1.0));
  return(x);
}

dym *mk_random_symmetric_dym(int rows,int cols)
{
  dym *x = mk_dym(rows,cols);
  dym *xt;
  int i,j;
  for ( i = 0 ; i < rows ; i++ )
    for ( j = 0 ; j <= i ; j++ )
    {
      double v = ((i==j)?3.0:1.0) * range_random(-5.0,1.0);
      dym_set(x,i,j,v);
      if ( j < i )
        dym_set(x,j,i,v);
    }

  xt = mk_dym_transpose(x);
  dym_plus(x,xt,x);
  free_dym(xt);
  return(x);
}

quadargs *mk_random_quadargs(int seed,int size)
{
  quadargs *qa = mk_empty_quadargs();
  int temp_seed = int_random(10000);
  am_srand(seed);
  qa -> a = mk_random_symmetric_dym(size,size);
  qa -> b = mk_random_dyv(size);
  qa -> fixed_atts = mk_ivec(0);
  qa -> fixed_vals = mk_dyv(0);
  qa -> lo = mk_constant_dyv(size,-1.0);
  qa -> hi = mk_constant_dyv(size,1.0);
  am_srand(temp_seed);
  return(qa);
}

/* kth row of result is result of doing optimization of
   quadargs in size dimensions, using seed+k. nodes[k] = number of
   search nodes visited */
dym *mk_quadres_dym(int start_seed,int size,int num_runs,ivec **nodes,double *r_mean_nodes,int *r_time)
{
  dym *result = mk_dym(num_runs,size);
  int k;
  int start_seconds = global_time();

  *nodes = mk_ivec(num_runs);
  for ( k = 0 ; k < num_runs ; k++ )
  {
    quadargs *qa = mk_random_quadargs(start_seed+k,size);
    dyv *xopt;
    double value;
    int j;

    Num_maxes = 0;
    Num_nodes = 0;

    xopt = mk_constrained_quadopt(qa->a,qa->b,0.0,qa->fixed_atts,qa->fixed_vals,
			          qa->lo,qa->hi,&value);

    printf("%2d %dd Maxs %3d Nds %4d Value %5.2f  Xopt ",
           start_seed+k,size,Num_maxes,Num_nodes,value);
    for ( j = 0 ; j < dyv_size(xopt) ; j++ )
      printf("%5.2f ",dyv_ref(xopt,j));
    printf("\n");

    copy_dyv_to_dym_row(xopt,result,k);
    ivec_set(*nodes,k,Num_nodes);
    free_quadargs(qa);
    free_dyv(xopt);
  }

  *r_mean_nodes = ivec_sum(*nodes) / (double) num_runs;
  *r_time = global_time() - start_seconds;

  return(result);
}
 
void run_quadres()
{
  ivec *nodes;
  int num_runs = 40;
  double mean_nodes;
  int time;
  dym *quadres = mk_quadres_dym(1,6,num_runs,&nodes,&mean_nodes,&time);
  char *resname = "sym";
  char fname[100];
  FILE *s;

  sprintf(fname,"%s.qrs",resname);
  s = safe_fopen(fname,"w");
  save_dym_to_file(s,quadres);
  fclose(s);
  printf("Saved %d xopts to file %s\n",num_runs,fname);

  sprintf(fname,"%s.mxs",resname);
  s = safe_fopen(fname,"w");
  save_ivec_to_file(s,nodes);
  fprintf(s,"Mean number of search nodes visited = %g\n",mean_nodes);
  fprintf(s,"Total number of seconds needed = %d\n",time);
  fprintf(stdout,"Mean number of search nodes visited = %g\n",mean_nodes);
  fprintf(stdout,"Total number of seconds needed = %d\n",time);

  fclose(s);
  printf("Saved %d number_of_nodes to file %s\n",num_runs,fname);

  free_ivec(nodes);
  free_dym(quadres);
}

void test_quadargs(quadargs *qa)
{
  dyv *xopt;
  double value;
  save_quadargs_to_file(stdout,qa);

  Num_maxes = 0;
  Num_nodes = 0;

  xopt = mk_constrained_quadopt(qa->a,qa->b,0.0,qa->fixed_atts,qa->fixed_vals,
			        qa->lo,qa->hi,&value);

  printf("The optimum has value %g and is at location xopt...\n",value);
  fprintf_dyv(stdout,"xopt",xopt,"\n");
  printf("Number of unconstrained maximizes = %d\n",Num_maxes);
  printf("Number of search nodes = %d\n",Num_nodes);

  free_dyv(xopt);
}

#ifdef NEVER
extern int runMLD;

void main(int argc,char *argv[])
{
#define QUADRES
  printf("Call int_from_args()...\n");
  runMLD = int_from_args("mld",argc,argv,0);
  printf("runMLD = %d\n",runMLD);
  printf("argc = %d\n",argc);

#ifdef QUADRES
  run_quadres();
#else
  char *filename = string_from_args("filename",argc,argv,"default.qa");
  FILE *s = fopen(filename,"r");
  quadargs *qa;

  Verbosity = input_realnum("Enter Verbosity (suggest 0 15 or 25)> ");

  if ( s == NULL )
  {
    int seed = 1234;
    int size = input_int("Can't load from file. Make random. What size?> ");
    qa = mk_random_quadargs(seed,size);
    s = safe_fopen("default.qa","w");
    save_quadargs_to_file(s,qa);
    fclose(s);
    printf("Saved your quadargs to file default.qa\n");
  }
  else
  {
    expo *ex;
    qa = mk_quadargs_from_file(s,&ex);
    fclose(s);
    if ( qa == NULL )
    {
      expo_render_plain(stdout,ex);
      free_expo(ex);
    }
  }

  if ( qa != NULL )
  {
    test_quadargs(qa);
    free_quadargs(qa);
  }
#endif
  wait_for_key();
  am_malloc_report();
  wait_for_key();
}

#endif


bool is_dym_sym_pos_def_and_nonsingular(dym *a)
{
  dym *z = mk_dym(dym_rows(a),dym_cols(a));
  bool result = attempt_cholesky_decomp(a,z);
  if ( result )
  {
    dyv *b = mk_constant_dyv(dym_rows(a),1.0);
    dyv *solve = mk_solve_lot_equation(z,b);
    result = solve != NULL;
    free_dyv(b);
    if ( solve != NULL ) free_dyv(solve);
  }
  free_dym(z);
  return(result);
}

bool is_dym_sym_neg_def_and_nonsingular(dym *a)
{
  dym *minus_a = mk_dym_scalar_mult(a,-1.0);
  bool result = is_dym_sym_pos_def_and_nonsingular(minus_a);
  free_dym(minus_a);
  return(result);
}

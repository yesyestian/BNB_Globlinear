
/*
   File:        amdyv.h
   Author:      Andrew W. Moore
   Created:     Thu Sep 15 21:01:12 EDT 1994
   Updated:     amdm was split into amdyv, amdym and svd by Frank Dellaert, Aug 14 1997
   Description: Header for Dynamically allocated and deallocated vectors

   Copyright 1996, Schenley Park Research
*/

#ifndef AMDYV_H
#define AMDYV_H

#include "ambs.h"      /* Very basic operations */

/*
%%%%%% DYVs (DYnamic Vectors)

dyvs are just one dimensional matrices. You can conveniently allocate
and free them too using the same conventions as for dynamic matrices.

dyv *mk_dyv(size) CREATES and returns one of the fellows.

int dyv_size(dyv *dv) returns the number of elements in the dyv

double dyv_ref(dv,i) returns the value of the i'th element.
  Indexing begins at 0: 0 <= i < dyv_size(dv).

All the numerical operations on dym's have counterpart numerical
operations on dyv's:

void constant_dyv(dyv *r_d,double v);
void zero_dyv(dyv *r_d);
dyv *mk_constant_dyv(int size,double v);
dyv *mk_zero_dyv(int size);
void dyv_scalar_mult(dyv *d, double alpha, dyv *r_d);
dyv *mk_dyv_scalar_mult(dyv *d,double alpha);
void dyv_scalar_add(dyv *d, double alpha, dyv *r_d);
dyv *mk_dyv_scalar_add(dyv *d,double alpha);
void copy_dyv(dyv *d, dyv *r_d);
dyv *mk_copy_dyv(dyv *d);
void dyv_plus(const dyv *d_1, const dyv *d_2, dyv *r_d);
dyv *mk_dyv_plus(dyv *a,dyv *b);
void dyv_subtract(dyv *d_1,dyv *d_2,dyv *r_d);
dyv *mk_dyv_subtract(dyv *a,dyv *b);
void dyv_sort(dyv *dv,dyv *r_dv);  It is fine, as always, if dy and r_dv are
                                   the same.
dyv *mk_dyv_sort(dyv *dv);


%%%%%%%%% Making small dyvs 
dyv *mk_dyv_1(double x) makes a 1-element dyv containing x as its only element
dyv *mk_dyv_2(double x,double y) 
         makes a 2-element dyv containing x as its 0th-indexed element
                                          y as its 1-index element
dyv *mk_dyv_3(double x,double y , double z) .... obvious.

%%%%%%%%% Complex vector operations

double dyv_scalar_product(dyv *a,dyv *b);
Returns a . b

double dyv_dsqd(dyv *a,dyv *b)
Returns (a - b).(a - b)

void fprintf_dyv(FILE *s,char *m1,const dyv *d,char *m2);

*/


typedef struct dyv_struct
{
  int dyv_code;
  int array_size;
  int size;
  double *farr;
} dyv, *dyv_ptr;

dyv *mk_dyv(int size);

void free_dyv(dyv *d);

void dyv_malloc_report(void);

int safe_dyv_size(const dyv *d);
double safe_dyv_ref(const dyv *d, int i);
void safe_dyv_set(dyv *d,int i,double value);
void safe_dyv_increment(dyv *d,int i,double value);

#ifdef AMFAST

#define dyv_size(d) ((d)->size)
#define dyv_ref(d,i) ((d)->farr[i])
#define dyv_set(d,i,v) (d)->farr[i] = (v)
#define dyv_increment(d,i,v) (d)->farr[i] += (v)

#else /* If not AMFAST */

#define dyv_size(d) (safe_dyv_size(d))
#define dyv_ref(d,i) (safe_dyv_ref(d,i))
#define dyv_set(d,i,v) (safe_dyv_set(d,i,v))
#define dyv_increment(d,i,v) (safe_dyv_increment(d,i,v))

#endif /* #ifdef AMFAST */

/*Added by Dan: Something I've wanted for a LONG time!*/
#define dyv_array_ref_ref(d,i,j) dyv_ref(dyv_array_ref(d,i),j);
#define dyv_array_ref_set(d,i,j,x) dyv_set(dyv_array_ref(d,i),j,x);

void dyv_increase_length(dyv *d,int extra_size);
void add_to_dyv(dyv *d,double new_val);

void dyv_remove(dyv *d,int index); /* Reduces size by 1, removes index'th 
                                      element. All elements to right of
                                      delete point copied one to left.
                                      See comments in amdm.c more details */
void dyv_remove_last_element(dyv *d); /* Reduce size by 1, remove 
                                        last rightmost elt */

/* Increases dv in length by 1 and shifts all elements
   with original index greater or equal to index one to the
   right and inserts val at index. */
void dyv_insert(dyv *dv,int index,double val);

void copy_dyv_to_farr(dyv *d, double *farr);

double *mk_farr_from_dyv(dyv *d);

void copy_farr_to_dyv(double *farr,int size,dyv *r_d);

dyv *mk_dyv_from_farr(double *farr,int size);

void copy_dyv_to_tdarr_row(dyv *dv,double **tdarr,int row);

void copy_dyv_to_tdarr_col(dyv *dv,double **tdarr,int col);

void copy_tdarr_row_to_dyv(double **tdarr,dyv *dv,int row);

dyv *mk_dyv_from_tdarr_row(double **tdarr,int row,int tdarr_cols);

void copy_tdarr_col_to_dyv(double **tdarr,dyv *dv,int col);

dyv *mk_dyv_from_tdarr_col(double **tdarr,int col,int tdarr_rows);


void constant_dyv(dyv *r_d,double v);

void zero_dyv(dyv *r_d);

bool zero_dyvp(dyv *d);

dyv *mk_constant_dyv(int size,double v);

dyv *mk_zero_dyv(int size);

void dyv_scalar_mult(const dyv *d, double alpha, dyv *r_d);

dyv *mk_dyv_scalar_mult(const dyv *d,double alpha);

void dyv_scalar_add(dyv *d, double alpha, dyv *r_d);

dyv *mk_dyv_scalar_add(dyv *d,double alpha);

void copy_dyv(const dyv *d, dyv *r_d);

dyv *mk_copy_dyv(const dyv *d);

void dyv_plus(const dyv *d_1, const dyv *d_2, dyv *r_d);

dyv *mk_dyv_plus(const dyv *a,const dyv *b);

void dyv_subtract(const dyv *d_1,const dyv *d_2,dyv *r_d);

dyv *mk_dyv_subtract(const dyv *a,const dyv *b);

double dyv_scalar_product(const dyv *a,const dyv *b);

double dyv_dsqd(const dyv *a,const dyv *b);

double dyv_magnitude(const dyv *a);

int index_in_sorted_dyv(dyv *a,double x);

void dyv_sort(dyv *dv,dyv *r_dv);

dyv *mk_dyv_sort(dyv *dv);


double dyv_sum(const dyv *dv);

double dyv_mean(const dyv *dv);
  /* Mean of all elements in dv */

double dyv_median(const dyv *dv); /* added by Artur Dubrawski on Aug 02 1996, efficiented by AWM */
  /* Median of all elements in dv */

double dyv_sdev(const dyv *dv);
  /* Sdev of all elements in dv */

double dyv_min(const dyv *dv);
  /* Min of all elements in dv. ERROR if dv sized zero */

double dyv_max(const dyv *dv);
  /* Max of all elements in dv. ERROR if dv sized zero */

int dyv_argmin(const dyv *dv);

int dyv_argmax(const dyv *dv);


dyv *mk_dyv_1(double x0);
dyv *mk_dyv_2(double x0,double x1);
dyv *mk_dyv_3(double x0,double x1,double x2);
dyv *mk_dyv_4(double x0,double x1,double x2,double x3);
dyv *mk_dyv_5(double x0,double x1,double x2,double x3,double x4);
dyv *mk_dyv_6(double x0,double x1,double x2,double x3,double x4,double x5);

double xt_diag_x_value(dyv *x, dyv *a);

dyv *mk_user_input_dyv(char *message,int dims);

dyv *mk_basic_dyv_from_args(char *name,int argc,char *argv[],int size);

dyv *mk_dyv_from_args(char *name,int argc,char *argv[],dyv *deflt);
/* COPIES in deflt (if so required) */
dyv *mk_dyv_from_args1(char *name,int argc,char *argv[],dyv *deflt);

double *mk_nrecipes_vector_from_dyv(dyv *d);

void copy_nrecipes_vector_to_dyv(double *nrfarr,dyv *d);

void free_nrecipes_vector(double *nrfarr,dyv *d);

void fprintf_dyv(FILE *s,char *m1,const dyv *d,char *m2);

#ifdef AMFAST
#define assert_dyv_shape(d,size,name)
#define check_dyv_code(d,name)
#else
void assert_dyv_shape(dyv *d,int size,char *name);
void check_dyv_code(const dyv *d, const char *name);
#endif /* #ifdef AMFAST */

/* Returns TRUE if any of x's elements are NaN */
/* Declared but not defined - sir 8/6/2000
bool dyv_isnan(dyv *x);*/

/* Returns TRUE if any elements are NaN or Inf */
bool dyv_is_ill_defined(dyv *x);

#define DYV_CODE 4509

#endif /* #ifndef AMDYV_H */

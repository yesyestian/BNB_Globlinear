/* 
   File:        amdmex.c
   Author:      Andrew W. Moore
   Created:     Mon Apr 10 21:39:26 EDT 1995
   Updated:     8 Dec 96
   Description: Extensions and advanced amdm stuff (no direct dym data access)

   Copyright 1996, Schenley Park Research

   This file contains advanced utility functions involving dyvs dyms and
   ivecs. It never accesses the data structures directly, so if the
   underlying representation of dyms and dyvs changes these won't need to.

   The prototypes of these functions are declared at the end of amdm.h
*/

#include "amdmex.h"

/* Makes a dyv of the same size as lo and hi (which must
   both be the same size such that
   dyv_ref(result,i) is uniformly randomly distributed between
   dyv_ref(lo,i) and dyv_ref(hi,i)
*/
dyv *mk_dyv_range_random(dyv *lo, dyv *hi)
{
  int i, n=dyv_size(lo);
  dyv *r = mk_dyv(n);

  assert_dyv_shape(hi,n,"mk_dyv_range_random"); 

  for (i=0; i<n; i++) {
    dyv_set(r,i,range_random(dyv_ref(lo,i),dyv_ref(hi,i)));
  }
  return r;
}

/* x := x with y appended on the end */
void append_to_dyv(dyv *x,dyv *y)
{
  int i;
  for ( i = 0 ; i < dyv_size(y) ; i++ )
    add_to_dyv(x,dyv_ref(y,i));
}

/* Return x with y appended on the end */
dyv *mk_dyv_append(dyv *x,dyv *y)
{
  dyv *z = mk_copy_dyv(x);
  append_to_dyv(z,y);
  return z;
}

/* x := x with y appended on the end */
void append_to_ivec(ivec *x,ivec *y)
{
  int i;
  for ( i = 0 ; i < ivec_size(y) ; i++ )
    add_to_ivec(x,ivec_ref(y,i));
}

/* Return x with y appended on the end */
ivec *mk_ivec_append(ivec *x,ivec *y)
{
  ivec *z = mk_copy_ivec(x);
  append_to_ivec(z,y);
  return z;
}

/* Forall (i,j) such that ilo <= i < ihi
                          jlo <= j < jhi (note the strict inequality on right)

   We do a[i][j] += delta
*/
void dym_increment_block(dym *a,int ilo,int jlo,int ihi,int jhi,double delta)
{
  int i,j;
  for ( i = ilo ; i < ihi ; i++ )
    for ( j = jlo ; j < jhi ; j++ )
      dym_increment(a,i,j,delta);
}

/* Forall i such that ilo <= i < ihi (note the strict inequality on right)

   We do a[i] += delta
*/
void dyv_increment_block(dyv *a,int ilo,int ihi,double delta)
{
  int i;
  for ( i = ilo ; i < ihi ; i++ )
    dyv_increment(a,i,delta);
}

void diag_times_dyv(dyv *a,dyv *b,dyv *x)
{
  int i;
  for ( i = 0 ; i < dyv_size(a) ; i++ )
    dyv_set(x,i,dyv_ref(a,i) * dyv_ref(b,i));
}

dyv *mk_diag_times_dyv(dyv *a,dyv *b)
{
  dyv *x = mk_dyv(dyv_size(a));
  diag_times_dyv(a,b,x);
  return(x);
}

void indices_of_sorted_dyv(const dyv *dv,ivec *iv)
/*
   NOTE: ivec structure (integer vectors) defined in sortind.ch
   PRE: dv and iv must be same size. iv's contents irrelevant.
   POST: iv contains sorted indexes into dv, so that
         forall iv[j] is the j'th smallest element of dv.

         thus forall i,j, (i < j) => dyv_ref(dv,iv[i]) <= dyv_ref(dv,iv[j])
          and iv contains a permutation of [0 ... iv->size-1]
*/
{
  int *iarr = am_malloc_ints(dyv_size(dv));
#ifndef AMFAST
  int i;
#endif
  indices_sort_realnums(dv->farr,dv->size,iarr);
  copy_iarr_to_ivec(iarr,dv->size,iv);
  am_free_ints(iarr,dv->size);
#ifndef AMFAST
  for ( i = 0 ; i < ivec_size(iv)-1 ; i++ )
  {
    int index1 = ivec_ref(iv,i);
    int index2 = ivec_ref(iv,i+1);
    if ( dyv_ref(dv,index1) > dyv_ref(dv,index2) )
    {
      fprintf_dyv(stdout,"dv",dv,"\n");
      fprintf_ivec(stdout,"iv",iv,"\n");
      printf("iv should be sorted indices of dyv, but consider\n"
	     "elements %d and %d of iv\n",i,i+1);
      my_error("mk_indices_of_sorted_dyv broken");
    }
  }
#endif
}

ivec *mk_indices_of_sorted_dyv(const dyv *dv)
{
  ivec *iv = mk_ivec(dyv_size(dv));
  indices_of_sorted_dyv(dv,iv);
  return(iv);
}

void indices_of_sorted_ivec(ivec *v,ivec *iv)
/*
   NOTE: ivec structure (integer vectors) defined in sortind.ch
   PRE: v and iv must be same size. iv's contents irrelevant.
   POST: iv contains sorted indexes into v, so that
         forall iv[j] is the j'th smallest element of v.

         thus forall i,j, (i < j) => ivec_ref(v,iv[i]) <= ivec_ref(v,iv[j])
          and iv contains a permutation of [0 ... iv->size-1]
*/
{
  int *iarr = am_malloc_ints(ivec_size(v));
#ifndef AMFAST
  int i;
#endif
  indices_sort_integers(v->iarr,v->size,iarr);
  copy_iarr_to_ivec(iarr,v->size,iv);
  am_free_ints(iarr,v->size);
#ifndef AMFAST
  for ( i = 0 ; i < ivec_size(iv)-1 ; i++ )
  {
    int index1 = ivec_ref(iv,i);
    int index2 = ivec_ref(iv,i+1);
    if ( ivec_ref(v,index1) > ivec_ref(v,index2) )
    {
      fprintf_ivec(stdout,"v",v,"\n");
      fprintf_ivec(stdout,"iv",iv,"\n");
      printf("iv should be sorted indices of ivec, but consider\n"
	     "elements %d and %d of iv\n",i,i+1);
      my_error("mk_indices_of_sorted_ivec broken");
    }
  }
#endif
}

ivec *mk_indices_of_sorted_ivec(ivec *v)
{
  ivec *iv = mk_ivec(ivec_size(v));
  indices_of_sorted_ivec(v,iv);
  return(iv);
}

void sorted_eigens_of_spd_dym(dym *d,dym *evectors,dyv *evalues)
/*
   PRE: d is a symmetric positive definite dym. (i.e. "spd")
        evectors is the same shape. evalues
        has size = d's width (and thus also height).

        initial contents of evectors and evalues irrelevant.

   POST: the i'th row of evectors contains the eigenvector with the
         'ith largest eigenvalue. Row 0 has the largest of all.

         Similarly, the i'th element of evalues contains the i'th
         largest eigenvalue.

   NOTE: This uses the fact that after SVD on a symmetric positive
         definite  matrix
         the eigenvectors are stored in the columns of the "v" matrix
         with the w's being the corresponding eigenvalues.

         SUBNOTE: If you have a symmetric matrix which is not posdef,
                  this function will return the correct thing except
                  you'll only have the magnitude of the eigenvalues---you'll 
                  have to work out the sign yourself (easily done).

           SUBSUBNOTE: What about non-symmetric matrixes. What will SVD do?
                       Good question.
*/
{
  dym *u,*v;
  dyv *w_vec;
  ivec *indices;
  int i;
  dym *evec_copy = mk_dym(dym_rows(d),dym_cols(d)); 
     /* In case d is same memory as evectors */

  if ( dym_cols(d) != dym_rows(d) )
    my_error("amdm.c, sorted_eigens_of_spd_dym: dym should be square");
/*
  if ( !is_dym_symmetric_positive_definite(d) )
    my_error("amdm.c, sorted_eigens_of_spd_dym: not spd");
*/
  make_svd_components(d,&u,&w_vec,&v);
  indices = mk_indices_of_sorted_dyv(w_vec);

  for ( i = 0 ; i < dym_rows(d) ; i++ )
  {
    int j = ivec_ref(indices,dym_rows(d) - i - 1);
              /* The index of the i'th LARGEST evalue */
    double ith_largest_eval = dyv_ref(w_vec,j);
    dyv *evec_of_ith_largest_eval = mk_dyv_from_dym_col(v,j);
    copy_dyv_to_dym_row(evec_of_ith_largest_eval,evec_copy,i);
    dyv_set(evalues,i,ith_largest_eval);
    free_dyv(evec_of_ith_largest_eval);
  }

  if ( dyv_min(evalues) < 0.0 )
  {
    fprintf_dym(stderr,"d",d,"\n");
    fprintf_dym(stderr,"evectors",evectors,"\n");
    fprintf_dyv(stderr,"evalues",evalues,"\n");
    fprintf(stderr,"Found a -ve eigenvalue, so d not sym pos def\n");
    my_error("amdm.c, sorted_eigens_of_spd_dym: not spd");
  }

  copy_dym(evec_copy,evectors);
  free_dym(u);
  free_dym(v);
  free_dyv(w_vec);
  free_ivec(indices);
  free_dym(evec_copy);
}

dyv *mk_sorted_eigens_of_spd_dym(dym *d)
/*
   PRE: d is a spd (symmetric positive definite).
   POST: the i'th element of result contains the i'th
         largest eigenvalue. dyv_ref(result,0) is the largest of all.
*/
{
  dym *evectors = mk_dym(dym_rows(d),dym_cols(d));
  dyv *evalues = mk_dyv(dym_rows(d));
  sorted_eigens_of_spd_dym(d,evectors,evalues);
  free_dym(evectors);
  return(evalues);
}

void random_unit_dyv(dyv *dv)
/*
   Fills up dv with a uniformly randomly chosen vector
   for which the magnitude is 1.
*/
{
  int i;
  double mag = 0.0;

  while ( mag < 1e-6 ) /* This loop will only be executed once, 999999
                          times out of 1000000 */
  {
    for ( i = 0 ; i < dyv_size(dv) ; i++ )
      dyv_set(dv,i,gen_gauss());
    mag = dyv_magnitude(dv);
  }

  dyv_scalar_mult(dv,1.0/mag,dv);
}

/* PRE: dyvs must be same size.
   Returns true if and only if all components of dx are >= corresponding
   components of dy
*/
bool dyv_weakly_dominates(const dyv *dx, const dyv *dy)
{
  dyv *diff = mk_dyv_subtract(dx,dy);
  bool result = (dyv_size(dx) == 0) ? TRUE : dyv_min(diff) >= 0.0;
  free_dyv(diff);
  return(result);
}

/* Sensible if args are NULL. False if different size */
bool dyv_equal(dyv *x1,dyv *x2)
{
  bool result = TRUE;

  if ( EQ_PTR(x1,x2) )
    result = TRUE;
  else if ( x1 == NULL || x2 == NULL )
    result = FALSE;
  else if ( dyv_size(x1) != dyv_size(x2) )
    result = FALSE;
  else
  {
    int i;
    for ( i = 0 ; result && i < dyv_size(x1) ; i++ ) 
      result = result && dyv_ref(x1,i) == dyv_ref(x2,i); /* Note == on doubles */
  }
  return(result);
}

bool dym_weakly_dominates(dym *dx,dym *dy)
{
  dym *diff = mk_dym_subtract(dx,dy);
  bool result = ((dym_rows(dx)*dym_cols(dx))==0) ? TRUE : dym_min(diff) >= 0.0;
  free_dym(diff);
  return(result);
}

/* Sensible if args are NULL. False if different size */
bool dym_equal(dym *x1,dym *x2)
{
  bool result = TRUE;

  if ( EQ_PTR(x1,x2) )
    result = TRUE;
  else if ( x1 == NULL || x2 == NULL )
    result = FALSE;
  else if ( dym_rows(x1) != dym_rows(x2) )
    result = FALSE;
  else if ( dym_cols(x1) != dym_cols(x2) )
    result = FALSE;
  else
  {
    int i,j;
    for ( i = 0 ; result && i < dym_rows(x1) ; i++ ) 
      for ( j = 0 ; result && j < dym_cols(x1) ; j++ ) 
        result = result && dym_ref(x1,i,j) == dym_ref(x2,i,j); /* Note == on doubles */
  }
  return(result);
}

/* Returns d in which all components are divided by the dyv's
   magnitude. If d is the zero vector, so is the result
*/
dyv *mk_normalize_dyv(dyv *d)
{
  dyv *result = mk_dyv(dyv_size(d));
  int i;
  double mag = dyv_magnitude(d);
  if ( mag > 1e-20 )
  {
    for ( i = 0 ; i < dyv_size(d) ; i++ )
      dyv_set(result,i,dyv_ref(d,i)/mag);
  }
  return(result);
}

/* It's okay for the args to be the same dyv */
void normalize_dyv(dyv *src,dyv *dest)
{
  dyv *norm = mk_normalize_dyv(src);
  copy_dyv(norm,dest);
  free_dyv(norm);
}

dyv *mk_random_unit_dyv(int size)
{
  dyv *result = mk_dyv(size);
  random_unit_dyv(result);
  return(result);
}

void dyv_clip_in_unit(dyv *dv)
/*
   Updates dv so that (writing dv' as constents after and dv as contents
   before)

    forall i    d'[i] = max( 0.0 , min( dv[i] , 1.0 ) )
*/
{
  int i;
  for ( i = 0 ; i < dyv_size(dv) ; i++ )
    dyv_set(dv,i,real_max(0.0,real_min(1.0,dyv_ref(dv,i))));
}

dyv *mk_dyv_clip_in_unit(dyv *dv)
{
  dyv *result = mk_copy_dyv(dv);
  dyv_clip_in_unit(result);
  return(result);
}

  
void random_rotation(dym *d)
/*
   Let X(k) be a matrix with k rows and n columns in which each row
   is a unit length vector and each row is orthogonal.

   Let X(k) have been randomly generated from the space of all such matrices.

   Let p be a random unit vector in N space.

   Let y = (I - X^T X) p

   Let x(k+1) = y / sqrt(y^T y)
   
   Let X(k+1) be the matrix constructed by adding a row at the bottom of X(k)
   with bottom row x(k+1)

   Then X(k) is a matrix with k+1 rows and n columns in which each row
   is a unit length vector and each row is orthogonal.

   This algorithm invented by a smug Andrew. But presumably actually well
   known.
*/
{
  int k;
  int n = dym_rows(d);
  dym *x = mk_dym(0,n);
  dym *id = mk_identity_dym(n);

  for ( k = 0 ; k < n ; k++ )
  {
    dyv *y = mk_zero_dyv(n);
    double ymag = 0.0;
    dym *xt = mk_dym_transpose(x);
    dym *xtx = mk_dym_mult(xt,x);
    dym *i_minus_xtx = mk_dym_subtract(id,xtx);
    
    while ( ymag < 1e-4 )
    {
      dyv *p = mk_random_unit_dyv(n);
      dym_times_dyv(i_minus_xtx,p,y);
      ymag = dyv_magnitude(y);
      free_dyv(p);
    }

    dyv_scalar_mult(y,1/ymag,y);

    add_row(x);
    copy_dyv_to_dym_row(y,x,k);

    free_dyv(y);
    free_dym(xt);
    free_dym(xtx);
    free_dym(i_minus_xtx);
  }

  copy_dym(x,d);
  free_dym(x);
  free_dym(id);
}

void fprintf_oneline_dyv(FILE *s,const char *m1, const dyv *d, const char *m2)
{
  int i;
  fprintf(s,"%s ",m1);
  for ( i = 0 ; i < dyv_size(d) ; i++ )
    fprintf(s,"%8g%s",dyv_ref(d,i),(i==dyv_size(d)-1) ? "" : " ");
  fprintf(s,"%s",m2);
}

/* Makes a random subset of iv of size "k" */
ivec *mk_random_ivec_subset(ivec *iv,int k)
{
  ivec *temp = mk_copy_ivec(iv);
  ivec *subset = mk_ivec(k);
  int i;
  shuffle_ivec(temp);
  for ( i = 0 ; i < k ; i++ )
    ivec_set(subset,i,ivec_ref(temp,i));
  free_ivec(temp);
  return subset;
}

ivec *mk_ivec_from_dyv(dyv *d)
{
  ivec *iv = mk_ivec(dyv_size(d));
  int i;
  for ( i = 0 ; i < dyv_size(d) ; i++ )
    ivec_set(iv,i,(int) floor(0.5 + dyv_ref(d,i)));

  return(iv);
}

dyv *mk_dyv_from_ivec(ivec *iv)
{
  dyv *d = mk_dyv(ivec_size(iv));
  int i;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    dyv_set(d,i,(double) ivec_ref(iv,i));

  return(d);
}

ivec *mk_ivec_from_args(char *key,int argc,char *argv[],ivec *def)
{
  dyv *d1 = mk_dyv_from_ivec(def);
  dyv *d2 = mk_dyv_from_args(key,argc,argv,d1);
  ivec *iv = mk_ivec_from_dyv(d2);
  free_dyv(d1);
  free_dyv(d2);
  return(iv);
}

void copy_dym_row(dym *source,int source_row,dym *dest,int dest_row)
{
  int i;
  for ( i = 0 ; i < dym_cols(source) ; i++ )
    dym_set(dest,dest_row,i,dym_ref(source,source_row,i));
}

void copy_dym_col(dym *source,int source_col,dym *dest,int dest_col)
{
  int i;
  for ( i = 0 ; i < dym_rows(source) ; i++ )
    dym_set(dest,i,dest_col,dym_ref(source,i,source_col));
}

void swap_dym_rows(dym *dm,int i,int j)
/* Swaps rows i and j. Doesn't screw up if same row */
{
  if ( i != j )
  {
    dyv *row_i = mk_dyv_from_dym_row(dm,i);
    dyv *row_j = mk_dyv_from_dym_row(dm,j);
    copy_dyv_to_dym_row(row_i,dm,j);
    copy_dyv_to_dym_row(row_j,dm,i);
    free_dyv(row_i);
    free_dyv(row_j);
  }
}

void swap_dym_cols(dym *dm,int i,int j)
/* Swaps rows i and j. Doesn't screw up if same row */
{
  if ( i != j )
  {
    dyv *col_i = mk_dyv_from_dym_col(dm,i);
    dyv *col_j = mk_dyv_from_dym_col(dm,j);
    copy_dyv_to_dym_col(col_i,dm,j);
    copy_dyv_to_dym_col(col_j,dm,i);
    free_dyv(col_i);
    free_dyv(col_j);
  }
}

/* rows may be NULL dneoting "use all rows" */
dyv *mk_sum_of_dym_rows(dym *dm,ivec *rows)
{
  int num_rows = (rows == NULL) ? dym_rows(dm) : ivec_size(rows);
  dyv *sum = mk_zero_dyv(dym_cols(dm));
  int j;
  for ( j = 0 ; j < num_rows ; j++ )
  {
    int row = (rows==NULL) ? j : ivec_ref(rows,j);
    dyv *this_row = mk_dyv_from_dym_row(dm,row);
    dyv_plus(sum,this_row,sum);
    free_dyv(this_row);
  }
  return(sum);
}

/* rows may be NULL dneoting "use all rows" */
dyv *mk_sum_of_dyv_array(const dyv_array *da)
{
  const int num_rows = dyv_array_size(da);
  dyv *sum = mk_zero_dyv(dyv_size(dyv_array_ref(da, 0)));
  int j;
  for ( j = 0 ; j < num_rows ; j++ )
  {
    const dyv *this_row = dyv_array_ref(da, j);
    dyv_plus(sum, this_row, sum);
  }
  return(sum);
}

/* rows may be NULL dneoting "use all rows" */
dyv *mk_mean_of_dym_rows(dym *dm,ivec *rows)
{
  int num_rows = (rows == NULL) ? dym_rows(dm) : ivec_size(rows);
  dyv *mean = mk_sum_of_dym_rows(dm,rows);
  dyv_scalar_mult(mean,1.0 / int_max(1,num_rows),mean);
  return(mean);
}

/* rows may be NULL dneoting "use all rows" */
dyv *mk_sdev_of_dym_rows(dym *dm,ivec *rows)
{
  int num_rows = (rows == NULL) ? dym_rows(dm) : ivec_size(rows);
  int size = dym_cols(dm);
  dyv *mean = mk_mean_of_dym_rows(dm,rows);
  dyv *result = mk_zero_dyv(size);

  int i,k;
  for ( i = 0 ; i < size ; i++ )
    for ( k = 0 ; k < num_rows ; k++ )
    {
      int row = (rows==NULL) ? k : ivec_ref(rows,k);
      double d = dym_ref(dm,row,i) - dyv_ref(mean,i);
      dyv_increment(result,i,d * d);
    }
  for ( i = 0 ; i < size ; i++ )
    dyv_set(result,i,sqrt(dyv_ref(result,i)/int_max(1,dym_rows(dm)-1)));

  free_dyv(mean);
  return result;
}

dyv *mk_sum_of_dym_cols(dym *dm)
{
  dyv *sum = mk_zero_dyv(dym_rows(dm));
  int i;
  for ( i = 0 ; i < dym_cols(dm) ; i++ )
  {
    dyv *col = mk_dyv_from_dym_col(dm,i);
    dyv_plus(sum,col,sum);
    free_dyv(col);
  }
  return(sum);
}

dyv *mk_mean_of_dym_cols(dym *dm)
{
  dyv *mean = mk_sum_of_dym_cols(dm);
  dyv_scalar_mult(mean,1.0 / int_max(1,dym_cols(dm)),mean);
  return(mean);
}

ivec *mk_ivec_from_dym_col(dym *a,int col)
{
  dyv *x = mk_dyv_from_dym_col(a,col);
  ivec *result = mk_ivec_from_dyv(x);
  free_dyv(x);
  return result;
}

ivec *mk_ivec_from_dym_row(dym *a,int row)
{
  dyv *x = mk_dyv_from_dym_row(a,row);
  ivec *result = mk_ivec_from_dyv(x);
  free_dyv(x);
  return result;
}

/***** Now we'll play with dyv_arrays which are adjustable length
       arrays of dyvs
*****/

#define INITIAL_DYV_ARRAY_SIZE 10

dyv_array *mk_empty_dyv_array()
{
  dyv_array *da = AM_MALLOC(dyv_array);
  da -> size = 0;
  da -> array_size = INITIAL_DYV_ARRAY_SIZE;
  da -> array = AM_MALLOC_ARRAY(dyv_ptr,da->array_size);
  return(da);
}

void add_to_dyv_array(dyv_array *da, const dyv *dv)
/*
     Assume dyv_array is previously of size n. After this it is of size
   n+1, and the n+1'th element is a COPY of dv.
*/
{
  if ( da -> size == da -> array_size )
  {
    int new_size = 2 + 2 * da->array_size;
    dyv **new_array = AM_MALLOC_ARRAY(dyv_ptr,new_size);
    int i;
    for ( i = 0 ; i < da -> array_size ; i++ )
      new_array[i] = da->array[i];
    AM_FREE_ARRAY(da->array,dyv_ptr,da->array_size);
    da -> array = new_array;
    da -> array_size = new_size;
  }
  da->array[da->size] = (dv==NULL) ? NULL : mk_copy_dyv(dv);
  da->size += 1;
}

int dyv_array_size(const dyv_array *da)
{
  return(da->size);
}

dyv *safe_dyv_array_ref(const dyv_array *da,int index)
/*
     Returns a pointer (not a copy) to the index'th element stored in
   the dyv_array. Error if index < 0 or index >= size
*/
{
  dyv *result;
  if ( index < 0 || index >= dyv_array_size(da) )
  {
    result = NULL;
    my_error("dyv_array_ref");
  }
  else
    result = da->array[index];
  return(result);
}
  
void dyv_array_set(dyv_array *dva, int index, const dyv *dv)
{
  if ((index < 0) || (dva == NULL) || (index >= dva->size))
        my_error("dyv_array_set: called with incompatible arguments");
  if (dva->array[index] != NULL)
        free_dyv(dva->array[index]);
  dva->array[index] = (dv == NULL) ? NULL : mk_copy_dyv(dv);
}

void fprintf_dyv_array(FILE *s,char *m1,dyv_array *da,char *m2)
{
  if ( dyv_array_size(da) == 0 )
    fprintf(s,"%s = <dyv_array with zero entries>%s",m1,m2);
  else
  {
    int i;
    for ( i = 0 ; i < dyv_array_size(da) ; i++ )
    {
      char buff[100];
      sprintf(buff,"%s[%2d]",m1,i);
      fprintf_dyv(s,buff,dyv_array_ref(da,i),m2);
    }
  }
}

void free_dyv_array(dyv_array *da)
{
  int i;
  for ( i = 0 ; i < dyv_array_size(da) ; i++ )
    if ( da->array[i] != NULL )
      free_dyv(da->array[i]);
  AM_FREE_ARRAY(da->array,dyv_ptr,da->array_size);
  AM_FREE(da,dyv_array);
}

dyv_array *mk_copy_dyv_array(const dyv_array *da)
{
  dyv_array *new_ar = mk_empty_dyv_array();
  int i;

  for ( i = 0 ; i < dyv_array_size(da) ; i++ )
    add_to_dyv_array(new_ar,dyv_array_ref(da,i));

  return(new_ar);
}

void dyv_array_remove(dyv_array *dva,int index)
{
  int i;
  dyv *dv = dyv_array_ref(dva,index);
  if ( dv != NULL ) free_dyv(dv);
  for ( i = index ; i < dva->size-1 ; i++ )
    dva->array[i] = dva->array[i+1];
  dva->array[dva->size-1] = NULL;
  dva->size -= 1;
}

dyv_array *mk_array_of_zero_length_dyvs(int size)
{
  dyv_array *dva = mk_empty_dyv_array();
  dyv *dv = mk_dyv(0);
  int i;

  for (i = 0; i < size; i++)
        add_to_dyv_array(dva, dv);
  free_dyv(dv);
  return(dva);
}

/***** Now we'll play with dym_arrays which are adjustable length
       arrays of dyms
*****/

#define INITIAL_DYM_ARRAY_SIZE 10

dym_array *mk_empty_dym_array()
{
  dym_array *da = AM_MALLOC(dym_array);
  da -> size = 0;
  da -> array_size = INITIAL_DYM_ARRAY_SIZE;
  da -> array = AM_MALLOC_ARRAY(dym_ptr,da->array_size);
  return(da);
}

void add_to_dym_array(dym_array *da,dym *dm)
/*
     Assume dym_array is previously of size n. After this it is of size
   n+1, and the n+1'th element is a COPY of dm.
*/
{
  if ( da -> size == da -> array_size )
  {
    int new_size = 2 + 2 * da->array_size;
    dym **new_array = AM_MALLOC_ARRAY(dym_ptr,new_size);
    int i;
    for ( i = 0 ; i < da -> array_size ; i++ )
      new_array[i] = da->array[i];
    AM_FREE_ARRAY(da->array,dym_ptr,da->array_size);
    da -> array = new_array;
    da -> array_size = new_size;
  }
  da->array[da->size] = mk_copy_dym(dm);
  da->size += 1;
}

int dym_array_size(dym_array *da)
{
  return(da->size);
}

dym *dym_array_ref(dym_array *da,int index)
/*
     Returns a pointer (not a copy) to the index'th element stored in
   the dym_array. Error if index < 0 or index >= size
*/
{
  dym *result;
  if ( index < 0 || index >= dym_array_size(da) )
  {
    result = NULL;
    my_error("dym_array_ref");
  }
  else
    result = da->array[index];
  return(result);
}
  
void fprintf_dym_array(FILE *s,char *m1,dym_array *da,char *m2)
{
  if ( dym_array_size(da) == 0 )
    fprintf(s,"%s = <dym_array with zero entries>%s",m1,m2);
  else
  {
    int i;
    for ( i = 0 ; i < dym_array_size(da) ; i++ )
    {
      char buff[100];
      sprintf(buff,"%s[%2d]",m1,i);
      fprintf_dym(s,buff,dym_array_ref(da,i),m2);
    }
  }
}

void free_dym_array(dym_array *da)
{
  int i;
  for ( i = 0 ; i < dym_array_size(da) ; i++ )
    free_dym(da->array[i]);
  AM_FREE_ARRAY(da->array,dym_ptr,da->array_size);
  AM_FREE(da,dym_array);
}

dym_array *mk_copy_dym_array(dym_array *da)
{
  dym_array *new_ar = mk_empty_dym_array();
  int i;

  for ( i = 0 ; i < dym_array_size(da) ; i++ )
    add_to_dym_array(new_ar,dym_array_ref(da,i));

  return(new_ar);
}

/* 
   This function tries to create and return a dyv containing the specified
   numbers in the string_array sa.  In the simple case, format == NULL and
   we return a dyv with the same number of elements as sa, the ith element
   being the numeric value of the ith string in sa.  If format is non-NULL,
   then it should have the same length as sa, and contain only the characters
   'i' and '-'.  An 'i' signifies that the corresponding string in sa 
   should be included in the result dyv; a '-' specifies that it should be
   skipped.  Thus for a non-NULL format, the jth element of the result
   dyv is the numeric value of the jth non-ignored string in sa.

   If there is an error (because a substring didn't represent a number,
   or the format's length didn't equal the size of sa), this function
   returns NULL and sets err_mess to a newly created string describing the
   problem....  If the function succeeds, it will always return a
   non-NULL dyv (though the dyv could, in theory, have size 0) and will
   set err_mess to NULL.

   Added by AWM March 16th 1997: If the string consists of the single token "empty"
   "zero" or "none" in any case, a zero-length dyv is returned without complaint.
*/
dyv *mk_dyv_from_string_array_with_error_message(string_array *sa,
                                                 char *format,
                                                 char **err_mess)
{
  bool ok = TRUE;
  int sa_size = string_array_size(sa);
  dyv *result = NULL;
  bool empty;  
  *err_mess = NULL;

  if ( sa_size == 0 )
    empty = TRUE;
  else if ( sa_size == 1 )
  {
    char *f = string_array_ref(sa,0);
    empty = caseless_eq_string(f,"empty") || caseless_eq_string(f,"zero") || caseless_eq_string(f,"none");
  }
  else
    empty = FALSE;

  if ( format != NULL && !empty && (((int) strlen(format)) != sa_size))
  {
    char buff[200];

    ok = FALSE;
    sprintf(buff, "Expected %d fields, not %d.",
            (int) strlen(format),
            sa_size);
    *err_mess = mk_copy_string(buff);
  }
  else if ( empty )
    result = mk_dyv(0);
  else
  {
    int i;     /* Indexes string_array */
    int j = 0; /* indexes result */
    int result_size;

    if ( format == NULL )
      result_size = sa_size;
    else
      result_size = num_of_char_in_string(format, 'i');

    result = mk_dyv(result_size);

    for ( i = 0 ; ok && i < sa_size ; i++ )
    {
      bool use_me = (format == NULL || format[i] == 'i');
      if ( use_me )
      {
        char *s = string_array_ref(sa,i);
        if ( !is_a_number(s) )
        {
          char buff[200];

          ok = FALSE;
          sprintf(buff, "Expected a number in column %d, not '%s'.",
                  i + 1,
                  s);
          *err_mess = mk_copy_string(buff);
        }
        else
        {
          dyv_set(result,j,atof(s));
          j += 1;
        }
      }
    }

    if ( ok && j != result_size )
      my_error("Assert j must be result_size");
  }

  if ( !ok && result != NULL )
  {
    free_dyv(result);
    result = NULL;
  }

  return(result);
}

/* This function is identical to mk_dyv_from_string_array_with_error_message
   except that it doesn't create an error message.  */
dyv *mk_dyv_from_string_array(string_array *sa, char *format)
{
  char *tem_str;
  dyv *res = mk_dyv_from_string_array_with_error_message(sa, format, &tem_str);

  if (tem_str != NULL) free_string(tem_str);
  return(res);
}

/* 
   This function tries to create and return a dyv by reading the specified
   numbers from string.  Each comma or section of whitespace in the string
   is interpreted as a separator between two fields of the string.  In the
   simple case, format == NULL and this function assumes that string 
   simply contains a series of numbers (separated by commas or whitespace),
   and returns a dyv whose ith element is the ith such number.  
   
   If format is non-NULL, then it should have the same length as there are
   fields in string, and in addition format should contain only the characters
   'i' and '-'.  An 'i' signifies that the corresponding field in string
   should be included in the result dyv; a '-' specifies that it should be
   skipped.  Thus for a non-NULL format, the jth element of the result
   dyv is the numeric value of the jth non-ignored field in string.  (In
   theory, the ignored fields of string could contain non-numbers.)

   If there is an error (because a non-number is found in string where a 
   number was expected, or because the format's length didn't equal the 
   number of fields in the string), this function returns NULL and sets
   err_mess to a newly created string describing the problem.  If the 
   function succeeds, it will always return a non-NULL dyv (though the dyv 
   could, in theory, have size 0), and will set err_mess to NULL.

   Added by AWM March 16th 1997: If the string consists of the single token "empty"
   "zero" or "none" in any case, a zero-length dyv is returned without complaint.
*/
dyv *mk_dyv_from_string_with_error_message(char *string,
                                           char *format,
                                           char **err_mess)
{
  string_array *sa = mk_broken_string_using_seppers(string,",");
  dyv *res = mk_dyv_from_string_array_with_error_message(sa, format, err_mess);

  free_string_array(sa);
  return(res);
}

/* This function is identical to mk_dyv_from_string_with_error_message
   except that it doesn't create an error message.  */
dyv *mk_dyv_from_string(char *string,char *format)
{
  char *tem_str;
  dyv *res = mk_dyv_from_string_with_error_message(string, format, &tem_str);

  if (tem_str != NULL) free_string(tem_str);
  return(res);
}

/* Makes two dyvs from reading numbers out of a string.
   If format is NULL, it assumes the string consists of numbers
   separated by spaces and/or commas (which are ignored, except as separators)
   
   If format is non-null, it must be zero-terminated string in which
   each character is an 'i' or an 'o' or a -. The string is broken, again by
   replacing commas with spaces then breaking into substrings. If the
   j'th format character is 'i' includes the numeric value of j'th
   substring at right end of input_dyv. The length of result dyv is the number
   of 'i's in format. If the j'th format character is 'o' includes
   at the right of the output_dyv.

   If format is NULL, uses all the fields, and puts all but the last in
   input_dyv, and the last in output_dyv.

   If there's an error (becasue a substring didn't represent a number,
   of the number of substrings wasn't sufficient to satisfy all the
   i's in the format string) returns NULL
*/
void mk_io_dyvs_from_string(
    char *string,
    char *format,
    dyv **r_in_dyv,
    dyv **r_out_dyv
  )
{
  dyv *in;
  dyv *out;
  int i;

  if ( format == NULL )
  {
    dyv *both = mk_dyv_from_string(string,NULL);
    if ( both == NULL )
    {
      in = NULL;
      out = NULL;
    }
    else
    {
      int in_size = dyv_size(both)-1;
      in = mk_dyv(in_size);
      out = mk_dyv_1(dyv_ref(both,in_size));
      for ( i = 0 ; i < in_size ; i++ )
        dyv_set(in,i,dyv_ref(both,i));
      free_dyv(both);
    }
  }
  else
  {
    char *in_format = mk_copy_string(format);
    char *out_format = mk_copy_string(format);
    int i;
    for ( i = 0 ; format[i] != '\0' ; i++ )
    {
      if ( format[i] == 'o' )
      {
        in_format[i] = '-';
        out_format[i] = 'i';
      }
      else if ( format[i] == 'i' )
        out_format[i] = '-';
    }
    in = mk_dyv_from_string(string,in_format);
    out = mk_dyv_from_string(string,out_format);
    free_string(in_format);
    free_string(out_format);
  }

  if ( in == NULL && out != NULL )
  {
    free_dyv(out);
    out= NULL;
  }
  
  if ( out == NULL && in != NULL )
  {
    free_dyv(in);
    in= NULL;
  }
  
  *r_in_dyv = in;
  *r_out_dyv = out;
}
/* Changed by AWM March 16th 1997 to turn an empty
   dyv into a 1-element string array: "empty" */
string_array *mk_string_array_from_dyv(dyv *d)
{
  string_array *sa = mk_string_array(dyv_size(d));
  if ( dyv_size(d) == 0 )
    add_to_string_array(sa,"empty");
  else
  {
    int i;
    for ( i = 0 ; i < dyv_size(d) ; i++ )
    {
      char buff[100];
      sprintf(buff,"%g",dyv_ref(d,i));
      string_array_set(sa,i,buff);
    }
  }
  return(sa);
}

/* Makes a string of numbers, each separated by whitespace.
   No quotes or anything. Just numbers. */
/* Changed by AWM March 16th 1997 to turn an empty
   dyv into a 1-element string: "empty" */
char *mk_string_from_dyv(dyv *d)
{
  string_array *sa = mk_string_array_from_dyv(d);
  char *s = mk_string_from_string_array(sa);
  free_string_array(sa);
  return(s);
}

/* Returns a dyv_array of size "num_dyvs" in which the i'th
   element (forall 0 <= i < num_dyvs) is a vector of
   size "dyv_size" in which each element is 0.0 */
dyv_array *mk_dyv_array_of_zeroed_dyvs(int num_dyvs,int dyv_size)
{
  dyv *zero = mk_zero_dyv(dyv_size);
  dyv_array *result = mk_empty_dyv_array();
  int i;
  for ( i = 0 ; i < num_dyvs ; i++ )
    add_to_dyv_array(result,zero);
  free_dyv(zero);
  return result;
}

/* Returns a dym_array of size "num_dyms" in which the i'th
   element (forall 0 <= i < num_dyms) is a matrix of
   size "dym_size x dym_size" in which each element is 0.0 */
dym_array *mk_dym_array_of_zeroed_dyms(int num_dyms,int dym_size)
{
  dym *zero = mk_zero_dym(dym_size,dym_size);
  dym_array *result = mk_empty_dym_array();
  int i;
  for ( i = 0 ; i < num_dyms ; i++ )
    add_to_dym_array(result,zero);
  free_dym(zero);
  return result;
}

/* Returns a dym_array of size "num_dyms" in which the i'th
   element (forall 0 <= i < num_dyms) is a matrix of
   size "dym_size_r x dym_size_c" in which each element is 0.0 */
dym_array *mk_dym_array_of_zeroed_nonrect_dyms(int num_dyms, int dym_size_r,
											   int dym_size_c)
{
  dym *zero = mk_zero_dym(dym_size_r,dym_size_c);
  dym_array *result = mk_empty_dym_array();
  int i;
  for ( i = 0 ; i < num_dyms ; i++ )
    add_to_dym_array(result,zero);
  free_dym(zero);
  return result;
}

dyv *mk_midpoint_dyv(dyv *a,dyv *b)
{
  dyv *sum = mk_dyv_plus(a,b);
  dyv_scalar_mult(sum,0.5,sum);
  return(sum);
}

double median_of_three(double x,double y,double z)
{
  double result;
  if ( x <= y && x <= z )
    result = real_min(y,z);
  else if ( y <= x && y <= z )
    result = real_min(x,z);
  else
    result = real_min(x,y);
  return(result);
}

int find_index_of_kth_smallest(const dyv *x,int k)
{
  int size = dyv_size(x);
  ivec *indices = mk_identity_ivec(size);
  int result = -1;

  while ( result < 0 )
  {
    /* Invariant: we must find which member of indices has the
       kth smallest value and return it.

       size == size of indices */
    double pivot = median_of_three(dyv_ref(x,ivec_ref(indices,0)),
				   dyv_ref(x,ivec_ref(indices,size/2)),
				   dyv_ref(x,ivec_ref(indices,size-1)));
    ivec *lower = mk_ivec(0);
    ivec *equal = mk_ivec(0);
    ivec *higher = mk_ivec(0);
    int i;

    for ( i = 0 ; i < size ; i++ )
    {
      int index = ivec_ref(indices,i);
      double v = dyv_ref(x,index);
      if ( v < pivot )
        add_to_ivec(lower,index);
      else if ( v > pivot )
        add_to_ivec(higher,index);
      else
        add_to_ivec(equal,index);
    }

    if ( k < ivec_size(lower) )
    {
      free_ivec(equal);
      free_ivec(higher);
      free_ivec(indices);
      indices = lower;
    }
    else if ( k < ivec_size(lower) + ivec_size(equal) )
    {
      result = ivec_ref(equal,0);
      free_ivec(lower);
      free_ivec(higher);
      free_ivec(indices);
      indices = equal;
    }
    else
    {
      k -= (ivec_size(lower) + ivec_size(equal));
      free_ivec(lower);
      free_ivec(equal);
      free_ivec(indices);
      indices = higher;
    }

    size = ivec_size(indices);
    if ( size == 0 )
      my_error("Can't happen");
    else if ( size == 1 )
      result = ivec_ref(indices,0);
  }

  free_ivec(indices);
  return(result);
}

double dyv_kth_smallest(const dyv *d,int k)
{
  int index = find_index_of_kth_smallest(d,k);
  return(dyv_ref(d,index));
}

double dyv_median(const dyv *d)
{
  int size = dyv_size(d);
  double result;

  if ( size == 0 )
  {
    my_error("dyv_median: zero size is not allowed");
    result = -77.8;
  }
  else if ( size % 2 == 1 )
  {
    int k = (size - 1)/2;
    result = dyv_kth_smallest(d,k);
  }
  else
  {
    int k = size/2;
    result = 0.5 * ( dyv_kth_smallest(d,k-1) + dyv_kth_smallest(d,k) );
  }
 
  return(result);
}



/*************************************************************************/


/* Saves m to the file in a format that can later be read back with
   mk_dym_from_file. */
void save_dym_to_file(FILE *s, dym *m)
{
  int i, j;
  int rows = dym_rows(m);
  int cols = dym_cols(m);

  fprintf(s, "DynamicMatrix: %d rows, %d columns\n", rows, cols);
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
      fprintf(s, "%g ", dym_ref(m, i, j));
    if( cols )
      fprintf(s, "\n");
  }

  return;
}

/* Reads a dym from a file.  The dym must have been saved in the
   format produced by save_dym_to_file.  The file pointer must be
   at the start of the dym on entry, and will be at the end of the
   dym on exit.  

   If this succeeds, r_errmess is set to NULL.  Otherwise it is
   set to a message explaining the problem, and mk_dym_from_file
   returns NULL. */
dym *mk_dym_from_file(FILE *s, char **r_errmess)
{
  dym *m = NULL;
  string_array *sa;

  *r_errmess = NULL;
  sa = mk_string_array_from_line(s);  
  if (sa == NULL)
    *r_errmess = mk_copy_string("Failed when trying to read a matrix (started at end of file)");
  else if (string_array_size(sa) < 5 || 
           (strcmp(string_array_ref(sa, 0), "DynamicMatrix:") != 0) ||
           (strcmp(string_array_ref(sa, 2), "rows,") != 0) ||
           (strcmp(string_array_ref(sa, 4), "columns") != 0))
  {
    *r_errmess = mk_copy_string("Failed when reading matrix. "
                                "First line should be: DynamicMatrix: <n> rows, <m> columns");
  }
  if (!(*r_errmess))
  {
    int num_rows, num_cols=0;
    bool okayp;

    num_rows = int_from_string(string_array_ref(sa, 1), &okayp);
    if (okayp)
      num_cols = int_from_string(string_array_ref(sa, 3), &okayp);
    if (!okayp)
      *r_errmess = mk_copy_string("Failed when trying to read #rows/#columns of a matrix");
    else
    {
      int cur_row = 0;
      
      m = mk_dym(num_rows, num_cols);
      if( num_rows*num_cols )
	{
	  while (okayp && (cur_row < num_rows))
	    {
	      string_array *curline = mk_string_array_from_line(s);
	      
	      if (curline == NULL)
		{
		  okayp = FALSE;
		  *r_errmess = mk_copy_string("Failed when trying to read a matrix (reached end of file)");
		}
	      else if (string_array_size(curline) != num_cols)
		{
		  char buff[200];
		  okayp = FALSE;
		  
		  sprintf(buff, "Failed when trying to read a matrix: row %d has wrong number of columns", cur_row);
		  *r_errmess = mk_copy_string(buff);
		}
	      else
		{
		  int j;
		  
		  for (j = 0; j < num_cols; j++)
		    {         
		      dym_set(m, cur_row, j, 
			      double_from_string(string_array_ref(curline, j), &okayp));
		      if (!okayp)
			{
			  char buff[200];
			  
			  sprintf(buff, "Failed when trying to read element at row %d, column %d of matrix", 
				  cur_row, j);
			  *r_errmess = mk_copy_string(buff);
			}
		    }
		  cur_row++;
		}
	      if (curline != NULL) free_string_array(curline);
	    }
	}
    }
  }
  if (sa != NULL) free_string_array(sa);
  if ((*r_errmess != NULL) && (m != NULL))
    {
      free_dym(m);
      m = NULL;
    }
  return(m);
}

/* Saves dyv to the file in a format that can later be read back with
   mk_dyv_from_file. */
void save_dyv_to_file(FILE *s, dyv *v)
{
  int i;
  int size = dyv_size(v);

  fprintf(s, "DynamicVector: %d components\n", size );
  for (i = 0; i < size; i++)
    fprintf(s, "%g ", dyv_ref(v, i));
  if( size )
    fprintf(s, "\n");

  return;
}

/* Reads a dyv from a file.  The dyv must have been saved in the
   format produced by save_dyv_to_file.  The file pointer must be
   at the start of the dyv on entry, and will be at the end of the
   dyv on exit.  

   If this succeeds, r_errmess is set to NULL.  Otherwise it is
   set to a message explaining the problem, and mk_dyv_from_file
   returns NULL. */
dyv *mk_dyv_from_file(FILE *s, char **r_errmess)
{
  dyv *v = NULL;
  string_array *sa;

  *r_errmess = NULL;
  sa = mk_string_array_from_line(s);  
  if (sa == NULL)
    *r_errmess = mk_copy_string("Failed when trying to read a vector (started at end of file)");
  else if (string_array_size(sa) < 3 || 
           (strcmp(string_array_ref(sa, 0), "DynamicVector:") != 0) ||
           (strcmp(string_array_ref(sa, 2), "components") != 0))
  {
    *r_errmess = mk_copy_string("Failed when trying to read a vector. "
                                "First line should have the format: "
                                "DynamicVector: <n> components");
  }
  if (!(*r_errmess))
  {
    int size;
    bool okayp;

    size = int_from_string(string_array_ref(sa, 1), &okayp);
    if (!okayp)
      *r_errmess = mk_copy_string("Failed when trying to read #components of a vector");
    else
    {
      v = mk_dyv( size );

      if( size )
	{
	  string_array *curline = mk_string_array_from_line(s);
	  
	  if (curline == NULL)
	    {
	      okayp = FALSE;
	      *r_errmess = mk_copy_string("Failed when trying to read a vector (reached end of file)");
	    }
	  else if (string_array_size(curline) != size)
	    {
	      char buff[200];
	      okayp = FALSE;
	      
	      sprintf(buff, "Failed when trying to read a vector: wrong number of components");
	      *r_errmess = mk_copy_string(buff);
	    }
	  else
	    {
	      int j;
	      
	      for (j = 0; j < size; j++)
		{         
		  dyv_set( v, j, 
			   double_from_string(string_array_ref(curline, j), &okayp));
		  if (!okayp)
		    {
		      char buff[200];
		      
		      sprintf(buff, "Failed when trying to read %d-th element of a vector", j);
		      *r_errmess = mk_copy_string(buff);
		    }
		}
	    }
       
	  if (curline != NULL) free_string_array(curline);
	}
    }
  }

  if (sa != NULL) free_string_array(sa);
  if ((*r_errmess != NULL) && (v != NULL))
  {
    free_dyv(v);
    v = NULL;
  }
  return(v);
}


/* Saves ivec to the file in a format that can later be read back with
   mk_ivec_from_file. */
void save_ivec_to_file(FILE *s, ivec *v)
{
  int i;
  int size = ivec_size(v);

  fprintf(s, "IntegerVector: %d components\n", size );
  for (i = 0; i < size; i++)
    fprintf(s, "%d ", ivec_ref(v, i));
  if( size )
    fprintf(s, "\n");

  return;
}

/* Reads a ivec from a file.  The ivec must have been saved in the
   format produced by save_ivec_to_file.  The file pointer must be
   at the start of the ivec on entry, and will be at the end of the
   ivec on exit.  

   If this succeeds, r_errmess is set to NULL.  Otherwise it is
   set to a message explaining the problem, and mk_ivec_from_file
   returns NULL. */
ivec *mk_ivec_from_file(FILE *s, char **r_errmess)
{
  ivec *v = NULL;
  string_array *sa;

  *r_errmess = NULL;
  sa = mk_string_array_from_line(s);  
  if (sa == NULL)
    *r_errmess = mk_copy_string("Failed when trying to read a vector (started at end of file)");
  else if (string_array_size(sa) < 3 || 
           (strcmp(string_array_ref(sa, 0), "IntegerVector:") != 0) ||
           (strcmp(string_array_ref(sa, 2), "components") != 0))
  {
    *r_errmess = mk_copy_string("Failed when trying to read a vector. "
                                "First line should have the format: "
                                "IntegerVector: <n> components");
  }
  if (!(*r_errmess))
  {
    int size;
    bool okayp;

    size = int_from_string(string_array_ref(sa, 1), &okayp);
    if (!okayp)
      *r_errmess = mk_copy_string("Failed when trying to read #components of a vector");
    else
    {
      v = mk_ivec( size );

      if( size )
	{
	  string_array *curline = mk_string_array_from_line(s);
	  
	  if (curline == NULL)
	    {
	      okayp = FALSE;
	      *r_errmess = mk_copy_string("Failed when trying to read a vector (reached end of file)");
	    }
	  else if (string_array_size(curline) != size)
	    {
	      char buff[200];
	      okayp = FALSE;
	      
	      sprintf(buff, "Failed when trying to read a vector: wrong number of components");
	      *r_errmess = mk_copy_string(buff);
	    }
	  else
	    {
	      int j;
	      
	      for (j = 0; j < size; j++)
		{         
		  ivec_set( v, j, 
			    int_from_string(string_array_ref(curline, j), &okayp));
		  if (!okayp)
		    {
		      char buff[200];
		      
		      sprintf(buff, "Failed when trying to read %d-th element of a vector", j);
		      *r_errmess = mk_copy_string(buff);
		    }
		}
	    }
	
	  if (curline != NULL) free_string_array(curline);
	}
    }
  }
  
  if (sa != NULL) free_string_array(sa);
  if ((*r_errmess != NULL) && (v != NULL))
    {
      free_ivec(v);
      v = NULL;
    }
  return(v);
}



/* Saves ma to the file in a format that can later be read back with
   mk_dym_array_from_file. */
void save_dym_array_to_file(FILE *s, dym_array *ma)
{
  int i;
  int size = dym_array_size(ma);
  dym *d;

  fprintf(s, "DynamicMatrixArray: %d components\n", size );
  for (i = 0; i < size; i++)
  {
    d = dym_array_ref( ma, i );
    save_dym_to_file( s, d );
  }

  return;
}


/* Reads a dym_array from a file.  The dym_array must have been saved in the
   format produced by save_dym_array_to_file.  The file pointer must be
   at the start of the dym_array on entry, and will be at the end of the
   dym_array on exit.  

   If this succeeds, r_errmess is set to NULL.  Otherwise it is
   set to a message explaining the problem, and mk_dym_array_from_file
   returns NULL. */
dym_array *mk_dym_array_from_file(FILE *s, char **r_errmess)
{
  dym_array *ma = NULL;
  string_array *sa;
  dym *d;
  int j;
  char *temp_string;

  *r_errmess = NULL;
  sa = mk_string_array_from_line(s);  
  if (sa == NULL)
    *r_errmess = mk_copy_string("Failed when trying to read a dynamic matrix array (started at end of file)");
  else if (string_array_size(sa) < 3 || 
           (strcmp(string_array_ref(sa, 0), "DynamicMatrixArray:") != 0) ||
           (strcmp(string_array_ref(sa, 2), "components") != 0))
    {
      *r_errmess = mk_copy_string("Failed when trying to read a dynamic matrix array. "
                                  "First line should have the format: "
                                  "DynamicMatrixArray: <n> components");
    }
  if (!(*r_errmess))
    {
      int size;
      bool okayp;
      
      size = int_from_string(string_array_ref(sa, 1), &okayp);
      if (!okayp)
	*r_errmess = mk_copy_string("Failed when trying to read #components of a dynamic matrix array");
      else
	{
	  ma = mk_empty_dym_array();
	  
	  for (j = 0; j < size; j++)
	    {         
	      if ( (d = mk_dym_from_file( s, &temp_string)) != NULL )
		add_to_dym_array( ma, d );
	      else
		{
		  *r_errmess = mk_printf("Failed when trying to read %d-th element "
                                         "of a dynamic matrix array: %s", j,temp_string);
		  free_string(temp_string);
		}
	      free_dym( d );
	    }
	}
    }
  
  if (sa != NULL) free_string_array(sa);
  if ((*r_errmess != NULL) && (ma != NULL))
    {
      free_dym_array(ma);
      ma = NULL;
    }
  return(ma);
}



string_array *mk_string_array_from_argc_argv(int argc,char *argv[])
{
  string_array *sa = mk_string_array(argc);
  int i;
  for ( i = 0 ; i < argc ; i++ )
    string_array_set(sa,i,argv[i]);
  return(sa);
}

string_array *mk_string_array_from_stream_tokens(FILE *s)
{
  string_array *sa = mk_string_array(0);
  bool finished = FALSE;
  while ( !finished )
  {
    string_array *this_line = mk_string_array_from_line(s);

    if ( this_line == NULL )
      finished = TRUE;
    else if ( string_array_size(this_line) > 0 &&
              string_array_ref(this_line,0)[0] != '#' )
    {
      int j;
      for ( j = 0 ; j < string_array_size(this_line) ; j++ )
        add_to_string_array(sa,string_array_ref(this_line,j));
    }
    if ( this_line != NULL ) free_string_array(this_line);
  }
  return(sa);
}

string_array *mk_string_array_from_file_tokens(char *filename)
{
  FILE *s = safe_fopen(filename,"r");
  string_array *sa = mk_string_array_from_stream_tokens(s);
  fclose(s);
  return(sa);
}

void make_argc_argv_from_string_array(string_array *sa,int *actual_argc,
                                      char ***actual_argv)
{
  int size = string_array_size(sa);
  int i;
  *actual_argv = AM_MALLOC_ARRAY(char *,size);
  for ( i = 0 ; i < size ; i++ )
    (*actual_argv)[i] = mk_copy_string(string_array_ref(sa,i));
  *actual_argc = size;  
}


/* This function AM_MALLOCS and RETURNS (in actual_argc and
   actual_argv) a new argc and argv. These are usually the same
   as argc and argv. But if argfile <filename> appears on the
   commandline, reads in all the tokens in <filename> and strings them
   together to make a longer argc argv. These new elements are
   appended onto the end of copies of the original argc argv.

   When you are done with these you should call

     free_loaded_argc_argv(actual_argc,actual_argv)

   Example:

     If the contents of file plop are:
        # This is a comment
        nsamples 35
        name Andrew

   ...and if someone runs the program foo with

       foo height 35 argfile plop noisy t

    Then after load_actual_argc_argv is called it will be as though
    the following command line was used:

       foo height 35 noisy t nsamples 35 name Andrew
*/

void load_actual_argc_argv(int argc,char *argv[],
                           int *actual_argc,char ***actual_argv)
{
  string_array *sa = mk_string_array_from_argc_argv(argc,argv);
  int i = 0;
 
  while ( i < string_array_size(sa) && string_array_size(sa) < 5000 )
  {
    if ( caseless_eq_string(string_array_ref(sa,i),"argfile") )
    {
      if ( i == string_array_size(sa)-1 )
        my_error("argfile must be followed by a filename on the commandline");
      else
      {
        char *filename = string_array_ref(sa,i+1);
        string_array *file_tokens = mk_string_array_from_file_tokens(filename);
        int j;
        string_array_remove(sa,i);
        string_array_remove(sa,i);

        fprintf(stdout,"Loaded the following command line "
                "extension from file %s:\n",filename);
  
        for ( j = 0 ; j < string_array_size(file_tokens) ; j++ )
	{
          fprintf(stdout,"%s ",string_array_ref(file_tokens,j));
          add_to_string_array(sa,string_array_ref(file_tokens,j));
	}
        fprintf(stdout,"\n");
        free_string_array(file_tokens);
      }
    }
    i += 1;
  }

  make_argc_argv_from_string_array(sa,actual_argc,actual_argv);
  free_string_array(sa);
}

void free_loaded_argc_argv(int argc,char **argv)
{
  int i;
  for ( i = 0 ; i < argc ; i++ )
    free_string(argv[i]);
  AM_FREE_ARRAY(argv,char *,argc);
}

/* rows may be NULL denotinf "use all datapoints" */
void make_vector_limits(dym *x,ivec *rows,dyv **xlo,dyv **xhi)
{
  int k;
  int num_points = (rows==NULL) ? dym_rows(x) : ivec_size(rows);

  *xlo = NULL;
  *xhi = NULL;

  for ( k = 0 ; k < num_points ; k++ )
  {
    int row = (rows==NULL) ? k : ivec_ref(rows,k);

    if ( *xlo == NULL )
    {
      *xlo = mk_dyv_from_dym_row(x,row);
      *xhi = mk_dyv_from_dym_row(x,row);
    }
    else
    {
      int j;
      for ( j = 0 ; j < dym_cols(x) ; j++ )
      {
        dyv_set(*xlo,j,real_min(dyv_ref(*xlo,j),dym_ref(x,row,j)));
        dyv_set(*xhi,j,real_max(dyv_ref(*xhi,j),dym_ref(x,row,j)));
      }
    }
  }
}

void make_vector_limits_from_dyv_array(const dyv_array *da, dyv **xlo,dyv **xhi)
{
  int k;
  int num_points = dyv_array_size(da);

  *xlo = NULL;
  *xhi = NULL;

  for ( k = 0 ; k < num_points ; k++ )
  {
    const dyv *row = dyv_array_ref(da, k);

    if ( *xlo == NULL )
    {
      *xlo = mk_copy_dyv(row);
      *xhi = mk_copy_dyv(row);
    }
    else
    {
      int j;
      for ( j = 0 ; j < dyv_size(row) ; j++ )
      {
        dyv_set(*xlo,j,real_min(dyv_ref(*xlo,j),dyv_ref(row,j)));
        dyv_set(*xhi,j,real_max(dyv_ref(*xhi,j),dyv_ref(row,j)));
      }
    }
  }
}
    
void set_vector_limits_sensible(dyv *xlo,dyv *xhi)
{
  int j;
  for ( j = 0 ; j < dyv_size(xlo) ; j++ )
  {
    double lo,hi,delta;
    sensible_limits(dyv_ref(xlo,j),dyv_ref(xhi,j),&lo,&hi,&delta);
    dyv_set(xlo,j,lo);
    dyv_set(xhi,j,hi);
  }
}

/* rows may be NULL denoting "use all datapoints" */
void make_sensible_vector_limits(dym *x,ivec *rows,dyv **xlo,dyv **xhi)
{
  make_vector_limits(x,rows,xlo,xhi);
  set_vector_limits_sensible(*xlo,*xhi);
}

/*** Some basic file utilities. Should really be in amiv.ch ****/

/* A line is interesting if its not all white space and
the leftmost non whitespace character isnt # */
bool line_string_is_interesting(char *line_string)
{
  int i;
  char first_non_whitespace = ' ';
  char second_non_whitespace = ' ';
  bool result;

  for ( i = 0 ; first_non_whitespace == ' ' && line_string[i] != '\0' ; i++ )
  {
    if ( line_string[i] != ' ' && line_string[i] != '\t' && 
         line_string[i] != '\r' )
      first_non_whitespace = line_string[i];
    if (first_non_whitespace != '\0') second_non_whitespace = line_string[i+1];
  }
  result = ( first_non_whitespace != ' ' );
  /* we allow the special sequence '##' to be a "machine readable comment */
  if ((first_non_whitespace == '#') && (second_non_whitespace != '#'))
    result = FALSE;

  return(result);
}

/* Searches the file for the next line that isn't all whitespace and
   that doesn't have # as its first non-whitespace character. 

   If no-such line before file-end, returns NULL */
char *mk_next_interesting_line_string(FILE *s,int *line_number)
{
  char *line_string = NULL;
  bool finished = FALSE;

  while ( !finished )
  {
    line_string = mk_string_from_line(s);
    *line_number += 1;
    if ( line_string == NULL )
      finished = TRUE;
    else
      finished = line_string_is_interesting(line_string);

    if ( !finished && line_string != NULL )
    {
      free_string(line_string);
      line_string = NULL;
    }
  }

  return(line_string);
}

/* As above excepts breaks resulting line into a string array of tokens... */
string_array *mk_next_interesting_line(FILE *s,int *line_number)
{
  char *str = mk_next_interesting_line_string(s,line_number);
  string_array *sa = (str == NULL) ? NULL : mk_broken_string(str);
  if ( str != NULL ) free_string(str);
  return(sa);
}

void save_dyv_to_file_plain(char *fname,dyv *x,char **r_errmess)
{
  FILE *s = fopen(fname,"w");
  *r_errmess = NULL;
  if ( s == NULL )
    *r_errmess = mk_printf("Can't open %s for writing",fname);
  else
  {
    int i;
    for ( i = 0 ; i < dyv_size(x) ; i++ )
      fprintf(s,"%g ",dyv_ref(x,i));
    fclose(s);
  }
}

bool file_exists(char *fname)
{
  FILE *s = fopen(fname,"r");
  bool result = s != NULL;
  if ( s != NULL ) fclose(s);
  return result;
}

void remove_file(char *fname,char **r_errmess)
{
  *r_errmess = NULL;
  unlink(fname);
  if ( file_exists(fname) )
    *r_errmess = mk_printf("Couldn't remove %s",fname);
}

void execute_command(char *exec_string,char **r_errmess)
{
  *r_errmess = NULL;
  printf("Will execute command: %s\n",exec_string);
  system(exec_string);
  printf("Finished executing command: %s\n",exec_string);
}

dyv *mk_dyv_from_file_plain(char *fname,int size,char **r_errmess)
{
  FILE *s = fopen(fname,"r");
  dyv *x = NULL;
  *r_errmess = NULL;

  if ( s == NULL )
    *r_errmess = mk_printf("Can't open %s for reading",fname);
  else
  {
    string_array *sa = mk_string_array_from_stream_tokens(s);  
    x = mk_dyv_from_string_array_with_error_message(sa,NULL,r_errmess);
    if ( x != NULL && size >= 0 && dyv_size(x) != size )
    {
      *r_errmess = mk_printf("I wanted to load a vector of size %d from "
			     "file %s, but the one I got was of size %d",
			     size,fname,dyv_size(x));
      free_dyv(x);
      x = NULL;
    }
    free_string_array(sa);
  }
  return x;
}

/***********************************************************************/

/******* FILE UTILITIES *********/

bool contains_a_number(string_array *sa)
{
  bool result = FALSE;
  int i;
  for ( i = 0 ; !result && i < string_array_size(sa) ; i++ )
    result = is_a_number(string_array_ref(sa,i));
  return(result);
}

bool exload_white_space(char c)
{
  return(c == ' ' || c == '\t' || c == '\n' || c == '\r');
}

bool format_comma_separated(int line_format)
{
  return(line_format == COMMA_FORMAT || line_format == AUTO_FORMAT);
}

bool format_whitespace_separated(int line_format)
{
  return(line_format == WHITESPACE_FORMAT || line_format == AUTO_FORMAT);
}

/************* NEW LINE PARSING CODE *************/

/* If line_format is WHITESPACE then the line is read SPACE_STYLE
   if line_format is COMMA      then the line is read COMMA_STYLE
   if lineformat is  ANY        then
        if there's an unquoted , anywhere on the line then use COMMA_STYLE
                                                      else use SPACE_STYLE

   The line parser runs through a finite state machine. On
   each character it looks at the character type:

     S Space       - The character is the ' ' char
     C Comma       - The character is the ',' char
     A SingleQuote - The character is the '\'' char
     Q DoubleQuote - The character is the '\"' char
     T Token       - The character is something other than the above
     
   The line parser is building up an array of tokens. It begins with
   an empty array of tokens. It has a current token being built. It begins
   with the current token empty. After each character is read, it performs
   one of the following actions:

     ADD   Add the curent token to the array. Set the current token to empty
     PUSH  Put the current character at the end of the current token
     NIL   Do nothing
     DASH  Put a dash character at the end of the current token
     DP    Put a dash, then the current character at end of token
     UNKN  Add the UNKNOWN_STRING to the array. Clear current token


  COMMA_STYLE parsing:

       All whitespace to immediate left and right of commas is removed.
       All other contiguous blocks of whitespace are replaced with - symbols
         (outside quotes, N contiguous spaces are replaced with one -.
          inside quotes, N contiguous spaces are replaced with N -'s)
       The resulting tokens between commas are used.
       Empty string between commas is turned into UNKNOWN STRING
  
  SPACE_STYLE parsing:

       All whitespace inside quotes are turned to dashes
       All other CONTIGUOUS blocks of whitespace are collapsed to one space
       Then the resulting tokens between whitespaces are used.
*/
typedef struct parse_state
{
  int id;
  int s_action;  int s_next;
  int c_action;  int c_next;
  int a_action;  int a_next;
  int q_action;  int q_next;
  int t_action;  int t_next;
  int end_action;
} parse_state;

#define CGO  0
#define C1   1
#define C2   2
#define CQS  3
#define CQ   4
#define CAS  5
#define CA   6
#define SGO  7
#define S1   8
#define SQST 9
#define SQ   10
#define SAST 11
#define SA   12

#define ADD  0
#define PUSH 1
#define NIL  2
#define DASH 3
#define DP   4
#define UNKN 5

parse_state Parse_array[] =
/* State   S(act,next)  C(act,next)  A(act,next)  Q(act,next)  T(act,next) END(act)*/
{{ CGO ,   NIL ,CGO  ,  UNKN,CGO  ,  NIL ,CAS  ,  NIL ,CQS  ,  PUSH,C1   , UNKN},
 { C1  ,   NIL ,C2   ,  ADD ,CGO  ,  PUSH,C1   ,  PUSH,C1   ,  PUSH,C1   , ADD },
 { C2  ,   NIL ,C2   ,  ADD ,CGO  ,  DP  ,C1   ,  DP  ,C1   ,  DP  ,C1   , ADD },
 { CQS ,   DASH,CQ   ,  PUSH,CQ   ,  PUSH,CQ   ,  NIL ,CGO  ,  PUSH,CQ   , UNKN},
 { CQ  ,   DASH,CQ   ,  PUSH,CQ   ,  PUSH,CQ   ,  NIL ,C1   ,  PUSH,CQ   , ADD },
 { CAS ,   DASH,CA   ,  PUSH,CA   ,  NIL ,CGO  ,  PUSH,CA   ,  PUSH,CA   , UNKN},
 { CA  ,   DASH,CA   ,  PUSH,CA   ,  NIL ,C1   ,  PUSH,CA   ,  PUSH,CA   , ADD },
 { SGO ,   NIL ,SGO  ,  PUSH,S1   ,  NIL ,SAST ,  NIL ,SQST ,  PUSH,S1   , NIL },
 { S1  ,   ADD ,SGO  ,  PUSH,S1   ,  PUSH,S1   ,  PUSH,S1   ,  PUSH,S1   , ADD },
 { SQST,   DASH,SQ   ,  PUSH,SQ   ,  PUSH,SQ   ,  UNKN,SGO  ,  PUSH,SQ   , UNKN},
 { SQ  ,   DASH,SQ   ,  PUSH,SQ   ,  PUSH,SQ   ,  ADD ,SGO  ,  PUSH,SQ   , ADD },
 { SAST,   DASH,SA   ,  PUSH,SA   ,  UNKN,SA   ,  PUSH,SA   ,  PUSH,SA   , UNKN},
 { SA  ,   DASH,SA   ,  PUSH,SA   ,  ADD ,SA   ,  PUSH,SA   ,  PUSH,SA   , ADD }};

string_array *mk_string_array_from_parse(char *s,bool comma_style)
{
  int state = (comma_style) ? CGO : SGO;
  bool finished = FALSE;
  int s_ptr = 0;
  string_array *tokarray = mk_string_array(0);
  int currtok_size = strlen(s) + 1;
  char *currtok = AM_MALLOC_ARRAY(char,currtok_size);
  int currtok_ptr = 0;
  bool is_a_number=TRUE, neg_number=FALSE;
  int num_periods = 0, e_pos = 0;

  while ( !finished )
  {
    parse_state *ps = &(Parse_array[state]);
    char c = s[s_ptr];
    int action;
    int next;

    if ( state != ps->id ) my_error("Parse_array misspecified");

    if ( c == '\0' )
    {
      finished = TRUE;
      next = -1;
      action = ps->end_action;
    }
    else if ( c == ' '  ) { action = ps->s_action ; next = ps->s_next; }
    else if ( c == ','  ) { action = ps->c_action ; next = ps->c_next; }
    else if ( c == '\'' ) { action = ps->a_action ; next = ps->a_next; }
    else if ( c == '\"' ) { action = ps->q_action ; next = ps->q_next; }
    else                  { action = ps->t_action ; next = ps->t_next; }

    /*This section added by DTH on 9/28/99 to eliminate confusion with minus signs*/
    if(c=='-' && currtok_ptr!=0 && currtok_ptr!=e_pos)
      c = '_';
    if(is_a_number && action!=ADD && (action!=NIL || currtok_ptr!=0) 
        && !((c>='0' && c<='9')  || c=='-'  || c=='+' || (c=='.' && num_periods==0) 
          || ((c=='e' || c=='E') && !e_pos)) ) 
      is_a_number = FALSE;
    if(c=='.') num_periods++;
    if((c=='e' || c=='E') && !e_pos) 
      e_pos = currtok_ptr+1;

    switch ( action )
    {
      case ADD :
        currtok[currtok_ptr] = '\0';
        /*This if added by DTH on 9/28/99 to eliminate confusion with minus signs*/
        if(!is_a_number){
          if(currtok[0]=='-') currtok[0] = '_';
          if(currtok[e_pos]=='-') currtok[e_pos] = '_';
        }
        add_to_string_array(tokarray,currtok);
        currtok_ptr = 0;
        is_a_number = TRUE;
        num_periods = 0;
        e_pos = 0;
      break;
      case PUSH:
        currtok[currtok_ptr] = c;
        currtok_ptr += 1;
      break;
      case NIL :
        /* skip */
      break;
      case DASH:
        currtok[currtok_ptr] = '_';
        currtok_ptr += 1;
      break;
      case DP  :
        currtok[currtok_ptr] = '_';
        currtok_ptr += 1;
        currtok[currtok_ptr] = c;
        currtok_ptr += 1;
      break;
      case UNKN:
        add_to_string_array(tokarray,"?");
        currtok_ptr = 0;
        neg_number = FALSE;
        is_a_number = TRUE;
        num_periods = 0;
      break;
      default: my_error("ljdnlkjs"); break;
    }

    state = next;
    s_ptr += 1;
  }

  AM_FREE_ARRAY(currtok,char,currtok_size);
  return(tokarray);
}

#define UQ_START    0
#define UQ_MIDDLE   1
#define UQ_INSIDE_Q 2
#define UQ_INSIDE_A 3
#define UQ_STOP_YES 4
#define UQ_STOP_NO  5

bool line_has_unquoted_comma(char *string)
{
  int state = UQ_START;
  int i = 0;
  while ( state != UQ_STOP_YES && state != UQ_STOP_NO )
  {
    char c = string[i];
    if ( c == '\0' )
      state = UQ_STOP_NO;
    else
    {
      switch ( state )
      {
        case UQ_START:
          if ( c == ' ' ) state = UQ_START;
          else if ( c == '\"' ) state = UQ_INSIDE_Q;
          else if ( c == '\'' ) state = UQ_INSIDE_A;
          else if ( c == ',' ) state = UQ_STOP_YES;
          else state = UQ_MIDDLE;
        break;
        case UQ_MIDDLE:
          if ( c == ' ' ) state = UQ_START;
          else if ( c == ',' ) state = UQ_STOP_YES;
          else state = UQ_MIDDLE;
        break;
        case UQ_INSIDE_A:
          if ( c == '\'' ) state = UQ_START;
          else state = UQ_INSIDE_A;
        break;
        case UQ_INSIDE_Q:
          if ( c == '\"' ) state = UQ_START;
          else state = UQ_INSIDE_Q;
        break;
        default: my_error("wiudbiuwb"); break;
      }
    }
    i += 1;
  }
  return(state == UQ_STOP_YES);
}


string_array *mk_parse_data_line(char *string,int line_format)
{
  bool comma_style = (line_format == COMMA_FORMAT) ||
                     (line_format == AUTO_FORMAT && line_has_unquoted_comma(string));
  string_array *sa = mk_string_array_from_parse(string,comma_style);
  return(sa);
}

string_array *mk_next_tokens(FILE *s,int *line_number,int line_format)
{
  char *line_string = mk_next_interesting_line_string(s,line_number);
  string_array *tokens;
  
  if ( line_string == NULL )
    tokens = NULL;
  else
    tokens = mk_parse_data_line(line_string,line_format);

  if ( line_string != NULL )
    free_string(line_string);

  return(tokens);
}

string_array *mk_default_attribute_names(int num_atts)
{
  string_array *sa = mk_string_array(num_atts);
  int i;
  for ( i = 0 ; i < num_atts ; i++ )
  {
    char buff[100];
    sprintf(buff,"x%d",i+1);
    string_array_set(sa,i,buff);
  }
  return(sa);
} 

bool all_numeric(string_array *sa)
{
  int i;
  bool result = TRUE;
  for ( i = 0 ; result && i < string_array_size(sa) ; i++ )
    result = is_a_number(string_array_ref(sa,i));
  return result;
}

void add_to_x_from_string_array(dym *x,string_array *sa,char **r_errmess)
{
  dyv *newrow = mk_dyv_from_string_array_with_error_message(sa,NULL,r_errmess);
  if ( newrow != NULL )
  {
    if ( dyv_size(newrow) != dym_cols(x) )
      *r_errmess = mk_printf("I expected %d items in the row of datapoint "
			     "values, but in fact I found %d",
                             dym_cols(x),dyv_size(newrow));
    else
    {
      add_row(x);
      copy_dyv_to_dym_row(newrow,x,dym_rows(x)-1);
    }
  }
  if ( newrow != NULL ) free_dyv(newrow);
}

void make_attnames_and_dym_from_filename(char *filename,int argc,char *argv[],
					 string_array **r_attnames,
                                         dym **r_x,char **r_errmess)
{
  FILE *s = fopen(filename,"r");
  bool finished = FALSE;
  int linenum = 0;
  *r_errmess = NULL;
  *r_attnames = NULL;
  *r_x = NULL;
  
  if ( s == NULL )
    *r_errmess = mk_printf("Can't open %s for reading",filename);
  
  while ( !finished && *r_errmess == NULL )
  {
    string_array *sa = mk_next_tokens(s,&linenum,AUTO_FORMAT);
    if ( sa == NULL )
      finished = TRUE;
    else if ( *r_attnames == NULL && all_numeric(sa) )
    {
      int num_atts = string_array_size(sa);
      *r_attnames = mk_default_attribute_names(num_atts);
      *r_x = mk_dym(0,num_atts);
      add_to_x_from_string_array(*r_x,sa,r_errmess);
    }
    else if ( *r_attnames == NULL )
    {
      int num_atts = string_array_size(sa);
      *r_attnames = mk_copy_string_array(sa);
      *r_x = mk_dym(0,num_atts);
    }
    else
      add_to_x_from_string_array(*r_x,sa,r_errmess);

    if ( sa != NULL ) free_string_array(sa);
  }

  if ( *r_x == NULL && *r_errmess == NULL )
    *r_errmess = mk_printf("No data, not even attribute names!");

  if ( *r_errmess != NULL )
  {
    if ( *r_x != NULL ) free_dym(*r_x);
    if ( *r_attnames != NULL ) free_string_array(*r_attnames);
    *r_x = NULL;
    *r_attnames = NULL;
    if ( s != NULL )
    {
      char buff[1000];
      sprintf(buff,"Error on line %d of file %s: %s",
	      linenum,filename,*r_errmess);
      free_string(*r_errmess);
      *r_errmess = mk_copy_string(buff);
    }
  }

  if ( s != NULL ) fclose(s);
}

dym *mk_dym_from_filename(char *filename,char **r_errmess)
{
  string_array *attnames;
  dym *x;
  make_attnames_and_dym_from_filename(filename,0,NULL,&attnames,&x,r_errmess);
  if ( attnames != NULL ) free_string_array(attnames);
  return x;
}

dym *mk_dym_from_filename_simple(char *filename)
{
  char *errmess = NULL;
  dym *x = mk_dym_from_filename(filename,&errmess);
  if ( x == NULL ) my_error(errmess);
  return x;
}

void save_attnames_and_dym(FILE *s,string_array *attnames,dym *x)
{
  int i,j;
  int rows = dym_rows(x);
  int cols = dym_cols(x);
  for ( j = 0 ; j < cols ; j++ )
    fprintf(s,"%s%s",string_array_ref(attnames,j),(j==cols-1)?"\n":",");
  for ( i = 0 ; i < rows ; i++ )
    for ( j = 0 ; j < cols ; j++ )
      fprintf(s,"%g%s",dym_ref(x,i,j),(j==cols-1)?"\n":",");
}

void save_dym_to_filename(char *filename,dym *x)
{
  FILE *s = safe_fopen(filename,"w");
  int i,j;
  int rows = dym_rows(x);
  int cols = dym_cols(x);

  for ( i = 0 ; i < rows ; i++ )
    for ( j = 0 ; j < cols ; j++ )
      fprintf(s,"%g%s",dym_ref(x,i,j),(j==cols-1)?"\n":",");
  fclose(s);
}

/* Find least index i such that value = dyv_ref(iv,i).
  If not found, returns -1
*/
int find_index_in_dyv(dyv *dv,double value,double error)
{
  int result = -1;
  int i;

  for ( i = 0 ; i < dyv_size(dv) && result < 0 ; i++ )
    if (fabs(value-dyv_ref(dv,i))<=error) 
      result = i;
  return(result);
}

bool is_in_dyv(dyv *dv,double value,double error)
{
  return(find_index_in_dyv(dv,value,error) >= 0);
}


ivec_array *mk_array_of_empty_ivecs(int size)
{
  ivec_array *iva = mk_empty_ivec_array();
  int i;
  ivec *empty_iv = mk_ivec(0);
  for ( i = 0 ; i < size ; i++ )
    add_to_ivec_array(iva,empty_iv);
  free_ivec(empty_iv);
  return(iva);
}

int index_of_longest_ivec(ivec_array *iva)
{
  int result = -1;
  int size = ivec_array_size(iva);
  if ( size == 0 )
    my_error("index_of_longest_ivec zero length");
  else
  {
    int i;
    result = 0;
    for ( i = 1 ; i < size ; i++ )
    {
      if ( ivec_size(ivec_array_ref(iva,result)) < 
           ivec_size(ivec_array_ref(iva,i)) )
        result = i;
    }
  }
  return(result);
}


/***************** SIVEC ***************/

/* A sivec is a regular old ivec, except it is treated as a set of
   integers.

   An ivec is a legal sivec if it is sorted in increasing order with no
   duplicates.

   The following set of functions consititute a reasonable simple
   package of set-theory operations.

   Note that everything is as efficient as possible for a set package 
   except for adding a single element and deleting a single element, 
   which (because of our representation by means of sorted ivecs) could
   take time linear in set size. */

bool is_sivec(const ivec *iv)
{
  int size = ivec_size(iv);
  int i;
  bool result = TRUE;
  for ( i = 0 ; result && i < size-1 ; i++ )
    if ( ivec_ref(iv,i) >= ivec_ref(iv,i+1) )
      result = FALSE;
  return result;
}

#ifdef AMFAST
#define debugger_assert_is_sivec(siv)
#else
void debugger_assert_is_sivec(const ivec *siv)
{
  if ( !is_sivec(siv) )
  {
    fprintf_ivec(stdout,"siv",siv,"\n");
    printf("The above ivec was not a sivec (SetIVEC) because it was not\n"
           "sorted. It was passed to a sivec operation in amdmex.c\n"
           "Use the debugger to find the ofending call\n");
  }
}
#endif

/* Makes { 0 , 1 , ... size-1 } */
ivec *mk_identity_sivec(int size)
{
  return mk_identity_ivec(size);
}

/* Returns number of elements in sivec */
int sivec_size(const ivec *siv)
{
  return ivec_size(siv);
}

/* Returns the minimum value in sivec. 
   Time cost: Constant */
int sivec_min(const ivec *siv)
{
#ifndef AMFAST
  if ( ivec_size(siv) <= 0 )
    my_error("sivec_min : empty set");
#endif
  return ivec_ref(siv,0);
}

/* Returns the minimum value in sivec. 
   Time cost: Constant */
int sivec_max(const ivec *siv)
{
#ifndef AMFAST
  if ( ivec_size(siv) <= 0 )
    my_error("sivec_max : empty set");
#endif
  return ivec_ref(siv,ivec_size(siv)-1);
}

/* If siv has 0 elements returns 0
   If value > ivec_max(siv) returns size
   If value <= ivec_min(siv) returns 0
   Else returns index such that
      value <= ivec_ref(siv,index) 
      ivec_ref(siv,index-1) < value
      
   It returns the value such that ivec_insert(iv,index,value)
   would represent the set with value added to iv (assuming value
   wasn't already in iv). */
int find_sivec_insert_index(ivec *siv,int value)
{
  int size = ivec_size(siv);
  int result;

  debugger_assert_is_sivec(siv);

  if ( size == 0 )
    result = 0;
  else
  {
    int loval = ivec_ref(siv,0);
    if ( value <= loval ) 
      result = 0;
    else
    {
      int lo = 0;
      int hi = size-1;
      int hival = ivec_ref(siv,hi);
      if ( value > hival )
	      result = size;
      else
      {
        while ( hi > lo + 1 )
        {
          int mid = (lo + hi) / 2;
          int midval = ivec_ref(siv,mid);
          if ( midval < value )
          {
            lo = mid;
            loval = midval;
          }
          else
          {
            hi = mid;
            hival = midval;
          }
        }
        if ( loval == value )
          result = lo;
        else
          result = hi;
      }
    }
  }
#ifndef AMFAST
  {
    bool ok;
    if ( size == 0 )
      ok = result == 0;
    else if ( value > sivec_max(siv) )
      ok = result == size;
    else if ( value <= sivec_min(siv) )
      ok = result == 0;
    else
      ok = value <= ivec_ref(siv,result) &&
           ivec_ref(siv,result-1) < value;

    if ( !ok )
      my_error("find_sivec_insert_index() : bug\n");
  }
#endif

  return result;
}

/* Adds the element while maintaining legal siveckiness.
   (If element already there, no change)
   Time cost: O(size) */
void add_to_sivec(ivec *siv,int value)
{
  int insert_index = find_sivec_insert_index(siv,value);
  if ( insert_index >= ivec_size(siv) )
    add_to_ivec(siv,value);
  else if ( value != ivec_ref(siv,insert_index) )
    ivec_insert(siv,insert_index,value);
}

ivec *mk_add_to_sivec(ivec *siv,int value)
{
  ivec *result = mk_copy_ivec(siv);
  add_to_sivec(result,value);
  return result;
}

/* Returns -1 if the value does not exist in siv.
   Else returns index such that
      value == ivec_ref(siv,value) 
  Time cost: O(log(size)) */
int index_in_sivec(ivec *siv,int value)
{
  int index = find_sivec_insert_index(siv,value);
  if ( index >= ivec_size(siv) || 
       value != ivec_ref(siv,index) )
    index = -1;
  return index;
}

/* Returns true iff siv contains value
   Time cost: O(log(size)) */
bool is_in_sivec(ivec *siv,int value)
{
  return index_in_sivec(siv,value) >= 0;
}

void sivec_remove_at_index(ivec *siv,int index)
{
  ivec_remove(siv,index);
}

/* Does nothing if value is not in siv.
   If value is in siv, the sivec is updated to
   represent siv \ { value } */
void sivec_remove_value(ivec *siv,int value)
{
  int index = find_sivec_insert_index(siv,value);
  if ( index < ivec_size(siv) && value == ivec_ref(siv,index) )
    ivec_remove(siv,index);
}
  
/* Returns answer to A subset-of B?
   Returns true if and only if the set of integers in a is
   a subset of the set of integers in b */
bool sivec_subset(const ivec *siva,const ivec *sivb)
{
  bool result = ivec_size(siva) <= ivec_size(sivb) ;
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sivec(siva);
  debugger_assert_is_sivec(sivb);

  while ( result && ai < ivec_size(siva) )
  {
    int aval = ivec_ref(siva,ai);
    while ( ivec_ref(sivb,bi) < aval && bi < ivec_size(sivb) )
      bi += 1;
    if ( bi >= ivec_size(sivb) || ivec_ref(sivb,bi) > aval )
      result = FALSE;
    else
      ai += 1;
  }
  return result;
}


bool sivec_equal(const ivec *siva,const ivec *sivb)
{
  return ivec_equal(siva,sivb);
}

/* Returns TRUE iff A is a subset of B and A != B */
bool sivec_strict_subset(const ivec *siva,const ivec *sivb)
{
  return sivec_subset(siva,sivb) && !sivec_equal(siva,sivb);
}

ivec *mk_sivec_union(const ivec *siva,const ivec *sivb)
{
  ivec *result = mk_ivec(0);
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sivec(siva);
  debugger_assert_is_sivec(sivb);
  while ( ai < ivec_size(siva) && bi < ivec_size(sivb) )
  {
    int aval = ivec_ref(siva,ai);
    int bval = ivec_ref(sivb,bi);
    if ( aval == bval )
    {
      add_to_ivec(result,aval);
      ai += 1;
      bi += 1;
    }
    else if ( aval < bval )
    {
      add_to_ivec(result,aval);
      ai += 1;
    }
    else
    {
      add_to_ivec(result,bval);
      bi += 1;
    }
  }

  while ( ai < ivec_size(siva) )
  {
    add_to_ivec(result,ivec_ref(siva,ai));
    ai += 1;
  }

  while ( bi < ivec_size(sivb) )
  {
    add_to_ivec(result,ivec_ref(sivb,bi));
    bi += 1;
  }

  debugger_assert_is_sivec(result);
  return result;
}

/* Returns A \ B.
   This is { x : x in A and x not in B } */
ivec *mk_sivec_difference(const ivec *siva,const ivec *sivb)
{
  ivec *result = mk_ivec(0);
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sivec(siva);
  debugger_assert_is_sivec(sivb);
  while ( ai < ivec_size(siva) && bi < ivec_size(sivb) )
  {
    while ( ai < ivec_size(siva) && ivec_ref(siva,ai) < ivec_ref(sivb,bi) )
    {
      add_to_ivec(result,ivec_ref(siva,ai));
      ai += 1;
    }
    if ( ai < ivec_size(siva) )
    {
      while ( bi < ivec_size(sivb) && ivec_ref(sivb,bi) < ivec_ref(siva,ai) )
        bi += 1;
      while ( ai < ivec_size(siva) && 
	      bi < ivec_size(sivb) && 
	      ivec_ref(siva,ai) == ivec_ref(sivb,bi) )
      {
	ai += 1;
	bi += 1;
      }
    }
  }
  while ( ai < ivec_size(siva) )
  {
    add_to_ivec(result,ivec_ref(siva,ai));
    ai += 1;
  }
  debugger_assert_is_sivec(result);
  return result;
}

ivec *mk_sivec_intersection(const ivec *siva,const ivec *sivb)
{
  ivec *result = mk_ivec(0);
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sivec(siva);
  debugger_assert_is_sivec(sivb);
  while ( ai < ivec_size(siva) && bi < ivec_size(sivb) )
  {
    while ( ai < ivec_size(siva) && ivec_ref(siva,ai) < ivec_ref(sivb,bi) )
      ai += 1;
    if ( ai < ivec_size(siva) )
    {
      while ( bi < ivec_size(sivb) && ivec_ref(sivb,bi) < ivec_ref(siva,ai) )
        bi += 1;
      if ( bi < ivec_size(sivb) && ivec_ref(siva,ai) == ivec_ref(sivb,bi) )
      {
	add_to_ivec(result,ivec_ref(siva,ai));
	ai += 1;
	bi += 1;
      }
    }
  }
  debugger_assert_is_sivec(result);
  return result;
}

/* Returns TRUE iff A intersect B is empty. O(size) time */
bool sivec_disjoint(ivec *siva,ivec *sivb)
{
  bool result = TRUE;
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sivec(siva);
  debugger_assert_is_sivec(sivb);
  while ( result && ai < ivec_size(siva) && bi < ivec_size(sivb) )
  {
    while ( ai < ivec_size(siva) && ivec_ref(siva,ai) < ivec_ref(sivb,bi) )
      ai += 1;
    if ( ai < ivec_size(siva) )
    {
      while ( bi < ivec_size(sivb) && ivec_ref(sivb,bi) < ivec_ref(siva,ai) )
        bi += 1;
      if ( bi < ivec_size(sivb) && ivec_ref(siva,ai) == ivec_ref(sivb,bi) )
	result = FALSE;
    }
  }
  return result;
}

ivec *mk_sivec_from_ivec(ivec *iv)
{
  ivec *siv = mk_ivec(0);
  if ( ivec_size(iv) > 0 )
  {
    ivec *sorted = mk_ivec_sort(iv);
    int i;
    int prev = ivec_ref(sorted,0);
    add_to_ivec(siv,prev);
    for ( i = 1 ; i < ivec_size(sorted) ; i++ )
    {
      int value = ivec_ref(sorted,i);
      if ( value < prev ) 
	my_error("No way");
      else if ( value != prev )
	add_to_ivec(siv,value);

      prev = value;
    }
    free_ivec(sorted);
  }
  return siv;
}

ivec *mk_ivec_from_string(char *s)
{
  dyv *dv = mk_dyv_from_string(s,NULL);
  ivec *iv = mk_ivec_from_dyv(dv);
  free_dyv(dv);
  return iv;
}

/* Turns a space separated string into a sivec.
   Example: "3 1 4 1 5 9" ====> { 1 , 3 , 4 , 5 , 9 } */
ivec *mk_sivec_from_string(char *s)
{
  ivec *iv = mk_ivec_from_string(s);
  ivec *siv = mk_sivec_from_ivec(iv);
  free_ivec(iv);
  return siv;
}

void paf(char *s,ivec *siv)
{
  fprintf_ivec(stdout,s,siv,"\n");
  free_ivec(siv);
}
    
void pb(char *s,bool v)
{
  printf("%s = %s\n",s,(v)?"True":"False");
}

void sivec_main(int argc,char *argv[])
{
  if ( argc < 4 )
    printf("%s %s \"<integers>\" \"integers\"\n",argv[0],argv[1]);
  else
  {
    ivec *siva = mk_sivec_from_string(argv[2]);
    ivec *sivb = mk_sivec_from_string(argv[3]);
    fprintf_ivec(stdout,"siva",siva,"\n");
    fprintf_ivec(stdout,"sivb",sivb,"\n");
    paf("union",mk_sivec_union(siva,sivb));
    paf("difference",mk_sivec_difference(siva,sivb));
    paf("intersection",mk_sivec_intersection(siva,sivb));
    
    pb("subset",sivec_subset(siva,sivb));
    pb("equal",sivec_equal(siva,sivb));
    pb("strict_subset",sivec_strict_subset(siva,sivb));
    pb("disjoint",sivec_disjoint(siva,sivb));

    add_to_sivec(siva,5);
    sivec_remove_value(sivb,5);

    fprintf_ivec(stdout,"siva U { 5 }",siva,"\n");
    fprintf_ivec(stdout,"sivb \\ { 5 }",sivb,"\n");
    free_ivec(siva);
    free_ivec(sivb);
  }
}

/********** Utilities added by Artur. But is "mk_ranks_fast" just
            the same as "mk_indices_of_sorted_dyv"??? (-AWM) ***/

//----------------------- returns a random permutation of v
void shuffle_dyv(dyv *v)
{
  int size = dyv_size(v);
  int i;
  for ( i = 0 ; i < size ; i++ )
  {
    int j = int_random(size-i);
    if ( j > 0 )
    {
      double swap_me_1 = dyv_ref(v,i);
      double swap_me_2 = dyv_ref(v,j+i);
      dyv_set(v,i,swap_me_2);
      dyv_set(v,j+i,swap_me_1);
    }
  }
  return;
}

/*************
This does the same as above, but may be faster in case of large datasets
with very many ties in the processed attributes.
**************/
dyv *mk_ranks_fast( dyv *source )
{
  int n = dyv_size(source);

  if( n > 0 ) // precaution against zero-dimensional source dyvs
  {
    dyv *ranks = mk_dyv( n );
    ivec *ind = mk_indices_of_sorted_dyv( source );
    ivec *ties = mk_zero_ivec( n );
    int i,j;
    double rank=1.0;
    int itsatie = 0;

    for(j=n-1;j>0;j--)
    {
     if( dyv_ref(source,ivec_ref(ind,j)) == dyv_ref(source,ivec_ref(ind,j-1)) )
     {
       itsatie++;
       ivec_set( ties, j, itsatie );
     }
     else
     {
       if( itsatie )
       {
          ivec_set( ties, j, itsatie+1 );
          itsatie=0;
       }
     }       
    }
    if( itsatie )
      ivec_set( ties, 0, itsatie+1 );
       
    j=0;
    rank = 1.0;
    while( j<n )
    {
      if( (itsatie = ivec_ref(ties,j)) )
      {
        double medrank = rank+0.5*(itsatie-1);
        for(i=j;i<j+itsatie;i++)
        {
          dyv_set(ranks,ivec_ref(ind,i),medrank);
        }
        rank+=(double)itsatie;
        j+=itsatie;
     }
      else
      {
       dyv_set(ranks,ivec_ref(ind,j),rank);
        rank+=1.0;
        j++;
      }
    }

    free_ivec( ind );
    free_ivec( ties );
    return( ranks );
  }
  else
    return(NULL);
}


/***************** sosarray ***************/

bool string_less(char *s1,char *s2)
{
  return strcmp(s1,s2) < 0;
}

bool string_greater(char *s1,char *s2)
{
  return string_less(s2,s1);
}

bool string_leq(char *s1,char *s2)
{
  return !string_greater(s1,s2);
}

bool string_geq(char *s1,char *s2)
{
  return !string_less(s1,s2);
}

/* A sosarray is a regular old string_array, except it is treated as a set of
   integers.

   An string_array is a legal sosarray if it is sorted in increasing order with no
   duplicates.

   The following set of functions consititute a reasonable simple
   package of set-theory operations.

   Note that everything is as efficient as possible for a set package 
   except for adding a single element and deleting a single element, 
   which (because of our representation by means of sorted string_arrays) could
   take time linear in set size. */

bool is_sosarray(string_array *sa)
{
  int size = string_array_size(sa);
  int i;
  bool result = TRUE;
  for ( i = 0 ; result && i < size-1 ; i++ )
    if ( string_geq(string_array_ref(sa,i),string_array_ref(sa,i+1)) )
      result = FALSE;
  return result;
}

#ifdef AMFAST
#define debugger_assert_is_sosarray(sosarr)
#else
void debugger_assert_is_sosarray(string_array *sosarr)
{
  if ( !is_sosarray(sosarr) )
  {
    fprintf_string_array(stdout,"sosarr",sosarr,"\n");
    printf("The above string_array was not a sosarray (sorted string array) because it was not\n"
           "sorted. It was passed to a sosarray operation in amdmex.c\n"
           "Use the debugger to find the ofending call\n");
  }
}
#endif

/* Returns number of elements in sosarray */
int sosarray_size(string_array *sosarr)
{
  return string_array_size(sosarr);
}

char *sosarray_first(string_array *sosarr)
{
  return string_array_ref(sosarr,0);
}

char *sosarray_last(string_array *sosarr)
{
  return string_array_ref(sosarr,string_array_size(sosarr)-1);
}

/* If sosarr has 0 elements returns 0
   If value > string_array_max(sosarr) returns size
   If value <= string_array_min(sosarr) returns 0
   Else returns index such that
      value <= string_array_ref(sosarr,index) 
      string_array_ref(sosarr,index-1) < value
      
   It returns the value such that string_array_insert(sa,index,value)
   would represent the set with value added to sa (assuming value
   wasn't already in sa). */
int find_sosarray_insert_index(string_array *sosarr,char *string)
{
  int size = string_array_size(sosarr);
  int result;

  debugger_assert_is_sosarray(sosarr);

  if ( size == 0 )
    result = 0;
  else
  {
    char *loval = string_array_ref(sosarr,0);
    if ( string_leq(string,loval) ) 
      result = 0;
    else
    {
      int lo = 0;
      int hi = size-1;
      char *hival = string_array_ref(sosarr,hi);
      if ( string_greater(string,hival) )
	      result = size;
      else
      {
        while ( hi > lo + 1 )
        {
          int mid = (lo + hi) / 2;
          char *midval = string_array_ref(sosarr,mid);
          if ( string_less(midval,string) )
          {
            lo = mid;
            loval = midval;
          }
          else
          {
            hi = mid;
            hival = midval;
          }
        }
        if ( eq_string(loval,string) )
          result = lo;
        else
          result = hi;
      }
    }
  }
#ifndef AMFAST
  {
    bool ok;
    if ( size == 0 )
      ok = result == 0;
    else if ( string_greater(string,sosarray_last(sosarr)) )
      ok = result == size;
    else if ( string_leq(string,sosarray_first(sosarr)) )
      ok = result == 0;
    else
      ok = string_leq(string,string_array_ref(sosarr,result)) &&
           string_less(string_array_ref(sosarr,result-1),string);

    if ( !ok )
      my_error("find_sosarray_insert_index() : bug\n");
  }
#endif

  return result;
}

/* Adds the element while maintaining legal sosarraykiness.
   (If element already there, no change)
   Time cost: O(size) */
void add_to_sosarray(string_array *sosarr,char *string)
{
  int insert_index = find_sosarray_insert_index(sosarr,string);
  if ( insert_index >= string_array_size(sosarr) )
    add_to_string_array(sosarr,string);
  else if ( !eq_string(string,string_array_ref(sosarr,insert_index)) )
    insert_in_string_array(sosarr,insert_index,string);
}

/* Returns -1 if the string does not exist in sosarr.
   Else returns index such that
      string == string_array_ref(sosarr,string) 
  Time cost: O(log(size)) */
int index_in_sosarray(string_array *sosarr,char *string)
{
  int index = find_sosarray_insert_index(sosarr,string);
  if ( index >= string_array_size(sosarr) || 
       !eq_string(string,string_array_ref(sosarr,index)) )
    index = -1;
  return index;
}

/* Returns true iff sosarr contains string
   Time cost: O(log(size)) */
bool is_in_sosarray(string_array *sosarr,char *string)
{
  return index_in_sosarray(sosarr,string) >= 0;
}

void sosarray_remove_at_index(string_array *sosarr,int index)
{
  string_array_remove(sosarr,index);
}

/* Does nothing if string is not in sosarr.
   If string is in sosarr, the sosarray is updated to
   represent sosarr \ { string } */
void sosarray_remove_string(string_array *sosarr,char *string)
{
  int index = find_sosarray_insert_index(sosarr,string);
  if ( index < string_array_size(sosarr) && eq_string(string,string_array_ref(sosarr,index)) )
    string_array_remove(sosarr,index);
}
  
/* Returns answer to A subset-of B?
   Returns true if and only if the set of integers in a is
   a subset of the set of integers in b */
bool sosarray_subset(string_array *sosarra,string_array *sosarrb)
{
  bool result = string_array_size(sosarra) <= string_array_size(sosarrb) ;
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sosarray(sosarra);
  debugger_assert_is_sosarray(sosarrb);

  while ( result && ai < string_array_size(sosarra) )
  {
    char *aval = string_array_ref(sosarra,ai);
    while ( string_array_ref(sosarrb,bi) < aval && bi < string_array_size(sosarrb) )
      bi += 1;
    if ( bi >= string_array_size(sosarrb) || 
	     string_greater(string_array_ref(sosarrb,bi),aval) )
      result = FALSE;
    else
      ai += 1;
  }
  return result;
}

bool equal_string_array(string_array *sa1,string_array *sa2)
{
  int size = string_array_size(sa1);
  bool result = size = string_array_size(sa2);
  int i;
  for ( i = 0 ; result && i < size ; i++ )
    result = eq_string(string_array_ref(sa1,i),string_array_ref(sa2,i));
  return result;
}

bool sosarray_equal(string_array *sosarra,string_array *sosarrb)
{
  return equal_string_array(sosarra,sosarrb);
}

/* Returns TRUE iff A is a subset of B and A != B */
bool sosarray_strict_subset(string_array *sosarra,string_array *sosarrb)
{
  return sosarray_subset(sosarra,sosarrb) && !sosarray_equal(sosarra,sosarrb);
}

string_array *mk_sosarray_union(string_array *sosarra,string_array *sosarrb)
{
  string_array *result = mk_string_array(0);
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sosarray(sosarra);
  debugger_assert_is_sosarray(sosarrb);
  while ( ai < string_array_size(sosarra) && bi < string_array_size(sosarrb) )
  {
    char *aval = string_array_ref(sosarra,ai);
    char *bval = string_array_ref(sosarrb,bi);
    if ( eq_string(aval,bval) )
    {
      add_to_string_array(result,aval);
      ai += 1;
      bi += 1;
    }
    else if ( string_less(aval,bval) )
    {
      add_to_string_array(result,aval);
      ai += 1;
    }
    else
    {
      add_to_string_array(result,bval);
      bi += 1;
    }
  }

  while ( ai < string_array_size(sosarra) )
  {
    add_to_string_array(result,string_array_ref(sosarra,ai));
    ai += 1;
  }

  while ( bi < string_array_size(sosarrb) )
  {
    add_to_string_array(result,string_array_ref(sosarrb,bi));
    bi += 1;
  }

  debugger_assert_is_sosarray(result);
  return result;
}

/* Returns A \ B.
   This is { x : x in A and x not in B } */
string_array *mk_sosarray_difference(string_array *sosarra,string_array *sosarrb)
{
  string_array *result = mk_string_array(0);
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sosarray(sosarra);
  debugger_assert_is_sosarray(sosarrb);
  while ( ai < string_array_size(sosarra) && bi < string_array_size(sosarrb) )
  {
    while ( ai < string_array_size(sosarra) && 
	        string_less(string_array_ref(sosarra,ai),string_array_ref(sosarrb,bi)) )
    {
      add_to_string_array(result,string_array_ref(sosarra,ai));
      ai += 1;
    }
    if ( ai < string_array_size(sosarra) )
    {
      while ( bi < string_array_size(sosarrb) && 
	          string_less(string_array_ref(sosarrb,bi),string_array_ref(sosarra,ai)) )
        bi += 1;
      while ( ai < string_array_size(sosarra) && 
	      bi < string_array_size(sosarrb) && 
	      eq_string(string_array_ref(sosarra,ai),string_array_ref(sosarrb,bi)) )
      {
	    ai += 1;
	    bi += 1;
      }
    }
  }
  while ( ai < string_array_size(sosarra) )
  {
    add_to_string_array(result,string_array_ref(sosarra,ai));
    ai += 1;
  }
  debugger_assert_is_sosarray(result);
  return result;
}

string_array *mk_sosarray_intersection(string_array *sosarra,string_array *sosarrb)
{
  string_array *result = mk_string_array(0);
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sosarray(sosarra);
  debugger_assert_is_sosarray(sosarrb);
  while ( ai < string_array_size(sosarra) && bi < string_array_size(sosarrb) )
  {
    while ( ai < string_array_size(sosarra) && 
	        string_less(string_array_ref(sosarra,ai),string_array_ref(sosarrb,bi)) )
      ai += 1;
    if ( ai < string_array_size(sosarra) )
    {
      while ( bi < string_array_size(sosarrb) && 
	          string_less(string_array_ref(sosarrb,bi),string_array_ref(sosarra,ai)) )
        bi += 1;
      if ( bi < string_array_size(sosarrb) && 
	       eq_string(string_array_ref(sosarra,ai),string_array_ref(sosarrb,bi)) )
      {
	add_to_string_array(result,string_array_ref(sosarra,ai));
	ai += 1;
	bi += 1;
      }
    }
  }
  debugger_assert_is_sosarray(result);
  return result;
}

/* Returns TRUE iff A intersect B is empty. O(size) time */
bool sosarray_disjoint(string_array *sosarra,string_array *sosarrb)
{
  bool result = TRUE;
  int ai = 0;
  int bi = 0;
  debugger_assert_is_sosarray(sosarra);
  debugger_assert_is_sosarray(sosarrb);
  while ( result && ai < string_array_size(sosarra) && bi < string_array_size(sosarrb) )
  {
    while ( ai < string_array_size(sosarra) && 
	        string_less(string_array_ref(sosarra,ai),string_array_ref(sosarrb,bi)) )
      ai += 1;
    if ( ai < string_array_size(sosarra) )
    {
      while ( bi < string_array_size(sosarrb) && 
	          string_less(string_array_ref(sosarrb,bi),string_array_ref(sosarra,ai)) )
        bi += 1;
      if ( bi < string_array_size(sosarrb) && 
	       eq_string(string_array_ref(sosarra,ai),string_array_ref(sosarrb,bi)) )
	result = FALSE;
    }
  }
  return result;
}

string_array *mk_sosarray_from_string_array(string_array *sa)
{
  string_array *sosarr = mk_string_array(0);
  if ( string_array_size(sa) > 0 )
  {
    string_array *sorted = mk_sort_string_array(sa);
    int i;
    char *prev = string_array_ref(sorted,0);
    add_to_string_array(sosarr,prev);
    for ( i = 1 ; i < string_array_size(sorted) ; i++ )
    {
      char *string = string_array_ref(sorted,i);
      if ( string_less(string,prev) ) 
    	my_error("No way");
      else if ( !eq_string(string,prev) )
	    add_to_string_array(sosarr,string);

      prev = string;
    }
    free_string_array(sorted);
  }
  return sosarr;
}

string_array *mk_string_array_from_string(char *s)
{
  dyv *dv = mk_dyv_from_string(s,NULL);
  string_array *sa = mk_string_array_from_dyv(dv);
  free_dyv(dv);
  return sa;
}

/* Turns a space separated string into a sosarray.
   Example: "3 1 4 1 5 9" ====> { 1 , 3 , 4 , 5 , 9 } */
string_array *mk_sosarray_from_string(char *s)
{
  string_array *sa = mk_string_array_from_string(s);
  string_array *sosarr = mk_sosarray_from_string_array(sa);
  free_string_array(sa);
  return sosarr;
}

/*************** The following stuff is for displaying 
                 string arrays in tabular form **********/

/* Makes an 'LS' style string matrix. That means it takes the
   string array and puts each element into cells of a string
   matrix. The string matrix has "cols" columns. The order in
   which string_array elements are placed is

      sa[0]      sa[r+0]    ....    sa[(cols-1)r+0]
      sa[1]      sa[r+1]    ....    sa[(cols-1)r+1]
        :                                  :
        :                                  :
        :                                  :
      sa[r-1]    sa[2r-1]   ....    sa[(cols-1)r-1]

    where r is the least r such that r*cols >= string_array_size

   ...and some of the rightmost column might be filled with
      empty cells.
*/
string_matrix *mk_ls_style_string_matrix_given_cols(string_array *names,int cols)
{
  int row;
  int size = string_array_size(names);
  int rows = (size + cols - 1) / cols;
  string_matrix *sm = mk_string_matrix(rows,cols);
  for ( row = 0 ; row < rows ; row++ )
  {
    int col;
    for ( col = 0 ; col < cols ; col++ )
    {
      int index = row + col * rows;
      char *name = (index < size) ?
                   string_array_ref(names,index) : "";
      string_matrix_set(sm,row,col,name);
    }
  }
  return sm;
}

/* Returns the max string length in sa */
int string_array_max_length(string_array *sa)
{
  int result = 0;
  int i;
  for ( i = 0 ; i < string_array_size(sa) ; i++ )
    result = int_max(result,(int)strlen(string_array_ref(sa,i)));
  return result;
}

/* Makes an 'LS' style string matrix cleverly designed so that when
   printed it uses less than "max_chars" characters per line. (it auto-chooses
   and sizes the columns) */
string_matrix *mk_ls_style_string_matrix(string_array *names,int max_chars)
{
  bool ok = FALSE;
  string_matrix *prev_sm = NULL;
  int cols = 1;
  while ( !ok )
  {
    string_matrix *sm = mk_ls_style_string_matrix_given_cols(names,cols);
    string_array *tabular = mk_tabular_string_array(sm,NULL);
    int chars = string_array_max_length(tabular);
    if ( chars > max_chars )
      ok = TRUE;
    if ( cols == 1 || !ok )
    {
      if ( prev_sm != NULL ) free_string_matrix(prev_sm);
      prev_sm = sm;
      cols += 1;
    }
    else
      free_string_matrix(sm);
    free_string_array(tabular);
  }
  return prev_sm;
}

/* Prints the contents of string_array cleverly much in the same way
   that "short ls" in unix displays filenames. It's cleverly designed so that when
   printed it uses less than "max_chars" characters per line. (it auto-chooses
   and sizes the columns) */
void display_names(FILE *s,string_array *names)
{
  string_matrix *sm = mk_ls_style_string_matrix(names,78);
  render_string_matrix(s,"",sm);
  free_string_matrix(sm);
}

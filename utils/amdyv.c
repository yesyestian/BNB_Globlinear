
/*
   File:        amdyv.c
   Author:      Andrew W. Moore
   Created:     Thu Sep 15 21:01:13 EDT 1994
   Updated:     amdm was split into amdyv, amdym and svd by Frank Dellaert, Aug 14 1997
   Description: Header for Dynamically allocated and deallocated vectors

   Copyright 1996, Schenley Park Research
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "amdyv.h"     /* Dynamically allocated and deallocated vectors */
#include "amma.h"      /* Fast, non-fragmenting, memory management */
#include "amar.h"      /* Obvious operations on 1-d arrays */

#ifndef AMFAST

void check_dyv_code(const dyv *d, const char *name)
{
  if ( d == NULL )
  {
    fprintf(stderr,"NULL dyv passed in operation %s\n",name);
    my_error("dyv data structure");
  }
  if ( d->dyv_code != DYV_CODE )
  {
    fprintf(stderr,"Attempt to access a non-allocated DYnamic Vector\n");
    fprintf(stderr,"This is in the operation %s\n",name);
    my_error("dyv data structure error");
  }
  if ( d->array_size < 0 || d->size < 0 || d->size > d->array_size )
    my_error("d->array_size or d->size muddled. Found in check_dyv_code.");
}

#endif /* #ifdef AMFAST */


/*
* check_dyv_access is only called in safe_dyv_xxx
* It used to be non-functional with AMFAST,
* but Frank Dellaert reinstated it June 30 1997
*/

void check_dyv_access(const dyv *d,int i, char *name)
{
  check_dyv_code(d,name); /* non-functional if AMFAST */

  if ( i < 0 || i >= d->size )
  {
    fprintf(stderr,"In operation \"%s\"\n",name);
    fprintf(stderr,"the dyv (dynamic vector) has size = %d\n",d->size);
    fprintf(stderr,"You tried to use index i=%d\n",i);
    fprintf(stderr,"Here is the dyv that was involved:\n");
    fprintf_dyv(stderr,"dv",d,"\n");
    my_error("check_dyv_access");
  }
}


#ifndef AMFAST

void assert_dyv_shape(dyv *d,int size,char *name)
{
  check_dyv_code(d,name);

  if ( size != d->size )
  {
    fprintf(stderr,"In operation \"%s\"\n",name);
    fprintf(stderr,"the dyv (dynamic vector) has size = %d\n", d->size);
    fprintf(stderr,"But should have been predefined with the shape:\n");
    fprintf(stderr,"size = %d\n",size);
    my_error("assert_dyv_shape");
  }
}

#endif /* #ifdef AMFAST */

int Dyvs_mallocked = 0;
int Dyvs_freed = 0;

dyv *mk_dyv(int size)
{
  dyv *result = AM_MALLOC(dyv);
  result -> dyv_code = DYV_CODE;
  result -> array_size = size;
  result -> size = size;
  result -> farr = am_malloc_realnums(size);
  Dyvs_mallocked += 1;
  return(result);
}

double dyv_sum(const dyv *dv)
{
  double result = 0.0;
  int i;
  for ( i = 0 ; i < dyv_size(dv) ; i++ )
    result += dyv_ref(dv,i);
  return(result);
}

bool dyv_is_ill_defined(dyv *x)
{
  return is_ill_defined(dyv_sum(x));
}

void free_dyv(dyv *d)
{
  check_dyv_code(d,"free_dyv");
#ifndef AMFAST
  if ( dyv_is_ill_defined(d) ) my_error("free_dyv: dyv contained NaN's");
#endif

  d -> dyv_code = 7777;

  am_free_realnums(d->farr,d->array_size);
  am_free((char *)d,sizeof(dyv));

  Dyvs_freed += 1;
}

void dyv_malloc_report(void)
{
  if ( Dyvs_mallocked > 0 )
  {
    fprintf(stdout,"# Dynamic Vectors  (datatype dyv) currently allocated:  %d\n",
           Dyvs_mallocked - Dyvs_freed
          );
    if ( Dyvs_mallocked - Dyvs_freed != 0 )
    {
      fprintf(stdout,"#       Number of dyv allocations since program start:  %d\n",
             Dyvs_mallocked
            );
      fprintf(stdout,"#       Number of dyv frees       since program start:  %d\n#\n",
             Dyvs_freed
            );
    }
  }
}

double safe_dyv_ref(const dyv *d, int i)
{
  check_dyv_access(d,i,"dyv_ref");
  return(d->farr[i]);
}

void safe_dyv_set(dyv *d,int i,double value)
{
  check_dyv_access(d,i,"dyv_set");
  d->farr[i] = value;
}

void safe_dyv_increment(dyv *d,int i,double value)
{
  check_dyv_access(d,i,"dyv_increment");
  d->farr[i] += value;
}

void dyv_increase_length(dyv *d,int extra_size)
{
  for ( ; extra_size > 0 ; extra_size-- )
    add_to_dyv(d,0.0);
}

void add_to_dyv(dyv *d,double new_val)
{
  if ( d->array_size < 0 || d->size > d->array_size )
    my_error("dyv size or array size has got muddles. Talk to AWM");

  if ( d->size == d->array_size )
  {
    int new_array_size = 2 * d->size + 2;
    double *farr_new = AM_MALLOC_ARRAY(double,new_array_size);
    int i;
    for ( i = 0 ; i < d->size ; i++ )
      farr_new[i] = d->farr[i];
    AM_FREE_ARRAY(d->farr,double,d->size);
    d -> farr = farr_new;
    d -> array_size = new_array_size;
  }
  d->farr[d->size] = new_val;
  d -> size += 1;
}

void copy_dyv_to_farr(dyv *d, double *farr)
{
  copy_realnums(d->farr,farr,d->size);
}
  
double *mk_farr_from_dyv(dyv *d)
{
  double *result;
  check_dyv_code(d,"make_copy_farr");
  result = am_malloc_realnums(d->size);
  copy_dyv_to_farr(d,result);
  return(result);
}

void copy_farr_to_dyv(double *farr,int size,dyv *r_d)
{
  assert_dyv_shape(r_d,size,"copy_farr_to_dyv");
  copy_realnums(farr,r_d->farr,size);
}

dyv *mk_dyv_from_farr(double *farr,int size)
{
  dyv *result = mk_dyv(size);
  copy_farr_to_dyv(farr,size,result);
  return(result);
}

/***** Copying dyvs to and from rows and columns of tdarrs. And
       Making dyvs from rows and columns of tdarrs too. *********/

void copy_dyv_to_tdarr_row(dyv *dv,double **tdarr,int row)
{
  int i;
  for ( i = 0 ; i < dv->size ; i++ )
    tdarr[row][i] = dv->farr[i];
}

void copy_dyv_to_tdarr_col(dyv *dv,double **tdarr,int col)
{
  int i;
  for ( i = 0 ; i < dv->size ; i++ )
    tdarr[i][col] = dv->farr[i];
}

void copy_tdarr_row_to_dyv(double **tdarr,dyv *dv,int row)
{
  int i;
  for ( i = 0 ; i < dv->size ; i++ )
    dv->farr[i] = tdarr[row][i];
}

dyv *mk_dyv_from_tdarr_row(double **tdarr,int row,int tdarr_cols)
{
  dyv *result = mk_dyv(tdarr_cols);
  copy_tdarr_row_to_dyv(tdarr,result,row);
  return(result);
}

void copy_tdarr_col_to_dyv(double **tdarr,dyv *dv,int col)
{
  int i;
  for ( i = 0 ; i < dv->size ; i++ )
    dv->farr[i] = tdarr[i][col];
}

dyv *mk_dyv_from_tdarr_col(double **tdarr,int col,int tdarr_rows)
{
  dyv *result = mk_dyv(tdarr_rows);
  copy_tdarr_col_to_dyv(tdarr,result,col);
  return(result);
}

void constant_dyv(dyv *r_d,double v)
{
  check_dyv_code(r_d,"constant_dyv");
  set_realnums_constant(r_d->farr,r_d->size,v);
}

void zero_dyv(dyv *r_d)
{
  check_dyv_code(r_d,"zero_dyv");
  constant_dyv(r_d,0.0);
}

/* Returns TRUE iff d is a non-NULL dyv all of whose elements are zero.
   Added by Mary on 11 Dec 96.
*/
bool zero_dyvp(dyv *d)
{
  bool res = TRUE;

  if (d == NULL)
    res = FALSE;
  else
  {
    int i;

    for (i = 0; res && i < dyv_size(d); i++)
    {
      if (dyv_ref(d, i) != 0.0) res = FALSE;
    }
  }
  return (res);
}

dyv *mk_constant_dyv(int size,double v)
{
  dyv *result = mk_dyv(size);
  constant_dyv(result,v);
  return(result);
}

dyv *mk_zero_dyv(int size)
{
  dyv *result = mk_dyv(size);
  zero_dyv(result);
  return(result);
}

/********* Standard operations of dyvs *********/

void dyv_scalar_mult(const dyv *d, double alpha, dyv *r_d)
{
  int i, n = d->size;
  assert_dyv_shape(r_d,n,"dyv_scalar_mult");
  for ( i = 0 ; i < n ; i++ )
    r_d -> farr[i] = d->farr[i] * alpha;
}

dyv *mk_dyv_scalar_mult(const dyv *d, double alpha)
{
  dyv *result;
  check_dyv_code(d,"mk_dyv_scalar_mult");
  result = mk_dyv(d->size);
  dyv_scalar_mult(d,alpha,result);
  return(result);
}

void dyv_scalar_add(dyv *d, double alpha, dyv *r_d)
{
  int i;
  assert_dyv_shape(r_d,d->size,"dyv_scalar_add");
  for ( i = 0 ; i < r_d -> size ; i++ )
    r_d -> farr[i] = d->farr[i] + alpha;
}

dyv *mk_dyv_scalar_add(dyv *d,double alpha)
{
  dyv *result;
  check_dyv_code(d,"mk_dyv_scalar_add");
  result = mk_dyv(d->size);
  dyv_scalar_add(d,alpha,result);
  return(result);
}

void copy_dyv(const dyv *d, dyv *r_d)
{
  assert_dyv_shape(r_d,d->size,"copy_dyv");
  dyv_scalar_mult(d,1.0,r_d);
}
    
dyv *mk_copy_dyv(const dyv *d)
{
  check_dyv_code(d,"mk_copy_dyv");
  return(mk_dyv_scalar_mult(d,1.0));
}
    
void dyv_plus(const dyv *d_1, const dyv *d_2, dyv *r_d)
{
  int i;
  check_dyv_code(d_1,"dyv_plus (1st arg)");
  check_dyv_code(d_2,"dyv_plus (2nd arg)");

  if ( d_1 -> size != d_2 -> size )
  {
    fprintf_dyv(stderr,"d_1",d_1,"\n");
    fprintf_dyv(stderr,"d_2",d_2,"\n");
    my_error("dyv_plus: dyvs (DYnamic Vectors) different shape");
  }

  assert_dyv_shape(r_d,d_1->size,"dyv_plus");
  for ( i = 0 ; i < r_d -> size ; i++ )
    r_d -> farr[i] = d_1->farr[i] + d_2 -> farr[i];
}

dyv *mk_dyv_plus(const dyv *a,const dyv *b)
{
  dyv *result = mk_dyv(a->size);
  dyv_plus(a,b,result);
  return(result);
}

void dyv_subtract(const dyv *d_1,const dyv *d_2,dyv *r_d)
{
  int i;
  check_dyv_code(d_1,"dyv_subtract (1st arg)");
  check_dyv_code(d_2,"dyv_subtract (2nd arg)");

  if ( d_1 -> size != d_2 -> size )
  {
    fprintf_dyv(stderr,"d_1",d_1,"\n");
    fprintf_dyv(stderr,"d_2",d_2,"\n");
    my_error("dyv_subtract: dyvs (DYnamic Vectors) different shape");
  }

  assert_dyv_shape(r_d,d_1->size,"dyv_subtract");
  for ( i = 0 ; i < r_d -> size ; i++ )
    r_d -> farr[i] = d_1->farr[i] - d_2 -> farr[i];
}

dyv *mk_dyv_subtract(const dyv *a,const dyv *b)
{
  dyv *result; 
  check_dyv_code(a,"mk_dyv_subtract (1st arg)");
  check_dyv_code(b,"mk_dyv_subtract (2nd arg)");
  result = mk_dyv(a->size);
  dyv_subtract(a,b,result);
  return(result);
}

/***** More complex operations ******/

double dyv_scalar_product(const dyv *a,const dyv *b)
{
  int i,n= a->size;
  double result = 0.0;
  if ( b -> size != n)
  {
    fprintf(stderr,"dyv_scalar_product: sizes differ\n");
    wait_for_key();
    fprintf_dyv(stderr,"a",a,"\n");
    fprintf_dyv(stderr,"b",b,"\n");
    my_error("dyv_scalar_product: sizes differ");
  }

  for ( i = 0 ; i < n ; i++ )
    result += a->farr[i] * b->farr[i];
  return(result);
}

/* FD rewrote this on aug 17: does not call mk_dyv_subtract anymore */
double dyv_dsqd(const dyv *a,const dyv *b)
{
  int i,n= a->size;
  double result = 0.0;
  if ( b -> size != n)
  {
    fprintf(stderr,"dyv_dsqd: sizes differ\n");
    wait_for_key();
    fprintf_dyv(stderr,"a",a,"\n");
    fprintf_dyv(stderr,"b",b,"\n");
    my_error("dyv_dsqd: sizes differ");
  }

  for ( i = 0 ; i < n ; i++ )
    {
      double dip = a->farr[i] - b->farr[i];
      result += dip*dip;
    }
  return(result);
}

double dyv_magnitude(const dyv *a)
{
  double result = sqrt(dyv_scalar_product(a,a));
  return(result);
}

double dyv_mean(const dyv *dv)
{
  check_dyv_code(dv,"dyv_mean");
  return(doubles_mean(dv->farr,dv->size));
}

double dyv_sdev(const dyv *dv)
{
  double sum_sq = 0.0;
  int i;
  double mean = dyv_mean(dv);
  double result;

  for ( i = 0 ; i < dv->size ; i++ )
    sum_sq += real_square(dv->farr[i] - mean);
    
  result = sqrt(sum_sq / int_max(dv->size-1,1));

  return(result);
}

double dyv_min(const dyv *dv)
{
  if ( dyv_size(dv) < 1 )
    my_error("dyv_min: empty dyv");
  return(doubles_min(dv->farr,dv->size));
}

double dyv_max(const dyv *dv)
{
  if ( dyv_size(dv) < 1 )
    my_error("dyv_max: empty dyv");
  return(doubles_max(dv->farr,dv->size));
}

int dyv_argmin(const dyv *dv)
{
  if ( dyv_size(dv) < 1 )
    my_error("dyv_argmin: empty dyv");
  return(doubles_argmin(dv->farr,dv->size));
}

int dyv_argmax(const dyv *dv)
{
  if ( dyv_size(dv) < 1 )
    my_error("dyv_argmax: empty dyv");
  return(doubles_argmax(dv->farr,dv->size));
}

dyv *mk_dyv_1(double x0)
{
  dyv *result = mk_dyv(1);
  dyv_set(result,0,x0);
  return(result);
}

dyv *mk_dyv_2(double x0,double x1)
{
  dyv *result = mk_dyv(2);
  dyv_set(result,0,x0);
  dyv_set(result,1,x1);
  return(result);
}

dyv *mk_dyv_3(double x0,double x1,double x2)
{
  dyv *result = mk_dyv(3);
  dyv_set(result,0,x0);
  dyv_set(result,1,x1);
  dyv_set(result,2,x2);
  return(result);
}

dyv *mk_dyv_4(double x0,double x1,double x2,double x3)
{
  dyv *result = mk_dyv(4);
  dyv_set(result,0,x0);
  dyv_set(result,1,x1);
  dyv_set(result,2,x2);
  dyv_set(result,3,x3);
  return(result);
}

dyv *mk_dyv_5(double x0,double x1,double x2,double x3,double x4)
{
  dyv *result = mk_dyv(5);
  dyv_set(result,0,x0);
  dyv_set(result,1,x1);
  dyv_set(result,2,x2);
  dyv_set(result,3,x3);
  dyv_set(result,4,x4);
  return(result);
}

dyv *mk_dyv_6(double x0,double x1,double x2,double x3,double x4,double x5)
{
  dyv *result = mk_dyv(6);
  dyv_set(result,0,x0);
  dyv_set(result,1,x1);
  dyv_set(result,2,x2);
  dyv_set(result,3,x3);
  dyv_set(result,4,x4);
  dyv_set(result,5,x5);
  return(result);
}

dyv *mk_user_input_dyv(char *message,int dims)
{
  dyv *res = mk_dyv(dims);
  int i = 0;
  char buff[100];
  for ( i = 0 ; i < dims ; i++ )
  {
    sprintf(buff,"%s (Dyv Component %d)> ",message,i);
    dyv_set(res,i,input_realnum(buff));
  }

  return(res);
}

dyv *mk_basic_dyv_from_args(char *name,int argc,char *argv[],int size)
{
  int index = index_of_arg(name,argc,argv);
  dyv *result;

  if ( index < 0 )
    result = NULL;
  else
  {
    int i;
    bool ok = TRUE;
    result = mk_dyv(size);
    for ( i = 0 ; i < size && ok ; i++ )
    {
      int j = index + i + 1;
      if ( j >= argc || !is_a_number(argv[j]) )
        ok = FALSE;
      else
        dyv_set(result,i,atof(argv[j]));
    }

    if ( !ok )
    {
      free_dyv(result);
      result = NULL;
    }
  }

  return(result);
}

dyv *mk_dyv_from_args(char *name,int argc,char *argv[],dyv *deflt)
/* COPIES in deflt (if so required) */
{
  bool name_there = index_of_arg(name,argc,argv) >= 0;
  dyv *result;
  if ( !name_there )
    result = mk_copy_dyv(deflt);
  else
  {
    result = mk_basic_dyv_from_args(name,argc,argv,dyv_size(deflt));
    if ( result == NULL )
    {
      fprintf(stderr,"COMMAND LINE USER ERROR (it's YOUR fault)\n");
      fprintf(stderr,"...when attempting to read a dyv identified by\n");
      fprintf(stderr,"the name \"%s\". Perhaps a non-number, or the\n",name);
      fprintf(stderr,"command line finished before all args found?\n");
      fprintf_dyv(stderr,"deflt_dyv",deflt,"\n");
      my_error("mk_dyv_from_args()");
    }
  }

  return(result);
}

/* as above except the dyv is expected to occupy a single argument rather
   than a sequence of them
 */
dyv *mk_dyv_from_args1(char *name,int argc,char *argv[],dyv *deflt)
{
  int index = index_of_arg(name,argc,argv)+1;
  dyv *result;
  int i, nfound = 0, nsize = 100;
  double *farr,*tmp;
  char *tok, *valstr;
  int valstrlen;

  if (index < 1)
  {
    if (!deflt) return NULL;
    else        return (mk_copy_dyv(deflt));
  }
  /* make copy of our value-string, which strtok will modify */
  valstr = mk_copy_string(argv[index]);
  valstrlen = strlen(valstr);
  
  farr = (double *)am_malloc(nsize * sizeof(double));
  tok = strtok(valstr," ");
  if (!tok)
    result=NULL;
  else {
    sscanf(tok,"%lf",&farr[nfound++]);
    while((tok = strtok(NULL," ")))
      {
	sscanf(tok,"%lf",&farr[nfound++]);
	if (nfound == (nsize-1))
	  {
	    tmp = (double *)am_malloc(nsize*2*sizeof(double));
	    for(i=0;i<nfound;i++) tmp[i] = farr[i];
	    am_free((char *)farr,nsize*sizeof(double));
	    nsize *= 2;
	    farr = tmp;
	  }
      }
    result = mk_dyv_from_farr(farr,nfound);
    am_free((char *)farr,nsize*sizeof(double));
  }
  am_free((char *) valstr, 1+valstrlen);
  return result;
}

int index_in_sorted_dyv(dyv *d,double t){
  int i1 = dyv_size(d), i0 = 0;
  int i = (i1-i0)>>1;
  double n2 = dyv_ref(d,i), n1 = (i)? dyv_ref(d,i-1):-1;
  while((i1-i0)>1 && !(n2>=t && (!i || n1<t))){
    if(n2>t){
      i1 = i;
      i = i0+((i-i0)>>1);
    } else {
      i0 = i+1;
      i = i+((i1-i)>>1);
    }
  }
  return i;
}

void dyv_sort(dyv *dv,dyv *r_dv)
{
  int size;
  double *farr;

  check_dyv_code(dv,"dyv_sort (1st arg)");
  assert_dyv_shape(r_dv,dv->size,"dyv_sort");

  size = dyv_size(dv);
  farr = mk_farr_from_dyv(dv);
  sort_realnums(farr,size,farr);
  copy_farr_to_dyv(farr,size,r_dv);
  am_free_realnums(farr,size);
}

dyv *mk_dyv_sort(dyv *dv)
{
  dyv *result;
  check_dyv_code(dv,"mk_dyv_sort");
  result = mk_dyv(dyv_size(dv));
  dyv_sort(dv,result);
  return(result);
}


/**** Removal functions on dyvs ****/

/* Reduces the size of d by one.
   dyv_ref(d,index) disappears.
   Everything to the right of dyv_ref(d,index) is copied one to the left.

   Formally: Let dold be the dyv value beore calling this function
             Let dnew be the dyv value after calling this function.

PRE: dyv_size(dold) > 0
     0 <= index < dyv_size(dold)

POST: dyv_size(dnew) = dyv_size(dold)-1
      for j = 0 , 1, 2 ... index-1  : 
         dyv_ref(dnew,j) == dyv_ref(dold,j)

      for j = index , index+1 , ... dyv_size(dnew)-1:
         dyv_ref(dnew,j) == dyv_ref(dold,j+1)
*/
void dyv_remove(dyv *d,int index)
{
  int i;
  int dsize = dyv_size(d);

#ifndef AMFAST
  if ( dsize <= 0 ) my_error("dyv_remove: empty dyv");
  if ( index < 0 || index >= dsize ) my_error("dyv_remove: bad index");
#endif /* #ifndef AMFAST */

  for ( i = index ; i < dsize - 1 ; i++ )
    dyv_set(d,i,dyv_ref(d,i+1));
  d -> size -= 1;
}

/* Shrinks d by one element by removing the rightmost element. 
   Example:

  Before: d == ( 3 1 4 1 5 )
    dyv_remove_last_element(d)
  After d == ( 3 1 4 1 )
*/
void dyv_remove_last_element(dyv *d)
{
#ifndef AMFAST
  if (dyv_size(d) <= 0) my_error("dyv_remove_last_elt: empty dyv");
#endif /* #ifndef AMFAST */
  d->size -= 1;
}

/* Increases dv in length by 1 and shifts all elements
   with original index greater or equal to index one to the
   right and inserts val at index. */
void dyv_insert(dyv *dv,int index,double val)
{
  int i;
  add_to_dyv(dv,0);
  for ( i = dyv_size(dv)-1 ; i > index ; i-- )
    dyv_set(dv,i,dyv_ref(dv,i-1));
  dyv_set(dv,index,val);
}



double *mk_nrecipes_vector_from_dyv(dyv *d)
{
  double *farr = am_malloc_realnums(d->size+1);
  int i;
  for ( i = 0 ; i < d->size ; i++ )
    farr[i+1] = d->farr[i];
  return(farr);
}

void copy_nrecipes_vector_to_dyv(double *nrfarr,dyv *d)
{
  int i;
  for ( i = 0 ; i < d->size ; i++ )
    d->farr[i] = nrfarr[i+1];
}

void free_nrecipes_vector(double *nrfarr,dyv *d)
{
  am_free_realnums(nrfarr,d->size+1);
}

/*
 * following 2 dyv functions are DEFINED in amdym.c
 * since they have code common with fprintf_dym
 
void fprintf_dyv(FILE *s,char *m1,dyv *d,char *m2)
void pdyv(dyv *d)

 */

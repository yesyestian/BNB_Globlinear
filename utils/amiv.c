/* *
   File:        amiv.c
   Author:      Andrew W. Moore
   Created:     Sat Apr  8 18:48:26 EDT 1995
   Updated:     29 March 97
   Description: Integer and Boolean Dynamic vectors

   Copyright 1997, Schenley Park Research
*/

#include <stdio.h>
#include <math.h>
#include "amiv.h"
#include <limits.h>
#include "adgui.h"

#define IVEC_CODE 20541

int Ivecs_mallocked = 0;
int Ivecs_freed = 0;

#ifdef AMFAST

#define NOTHING_TO_DO

#define check_ivec_code(iv,name) NOTHING_TO_DO

#else

void check_ivec_code(const ivec *iv,char *name)
{
  if ( iv == NULL )
  {
    fprintf(stderr,"NULL ivec passed in operation %s\n",name);
    my_error("ivec data structure");
  }
  if ( iv->ivec_code != IVEC_CODE )
  {
    fprintf(stderr,"Attempt to access a non-allocated Integer Vector\n");
    fprintf(stderr,"This is in the operation %s\n",name);
    my_error("ivec data structure error");
  }
  if ( iv->array_size < iv->size )
    my_error("check_ivec_code: array_size and size muddled");
}

#endif

#ifdef AMFAST

#define check_ivec_access(iv, i, name) NOTHING_TO_DO

#else

void check_ivec_access(const ivec *iv, int i, char *name)
{
  check_ivec_code(iv, name);

  if ( i < 0 || i >= iv->size )
  {
    fprintf(stderr,"In operation \"%s\"\n",name);
    fprintf(stderr,"the ivec (int integer vector) has size = %d\n",iv->size);
    fprintf(stderr,"You tried to use index i=%d\n",i);
    fprintf(stderr,"Here is the ivec that was involved:\n");
    fprintf_ivec(stderr,"ivv",iv,"\n");
    my_error("check_ivec_access");
  }
  if ( iv->array_size < iv->size )
    my_error("check_ivec_code: array_size and size muddled");
}

#endif

void privec(ivec *iv)
{
  fprintf_ivec(stdout,"iv",iv,"\n");
}

void assert_ivec_shape(ivec *iv,int size,char *name)
{
  check_ivec_code(iv,name);

  if ( size != iv->size )
  {
    fprintf(stderr,"In operation \"%s\"\n",name);
    fprintf(stderr,"the ivec (int dynamic vector) has size = %d\n", iv->size);
    fprintf(stderr,"But should have been predefined with the shape:\n");
    fprintf(stderr,"size = %d\n",size);
    my_error("assert_ivec_shape");
  }
}

ivec *mk_ivec(int size)
{
  ivec *result = AM_MALLOC(ivec);
  if ( size < 0 ) my_error("mk_ivec : size < 0 illegal");
  result -> ivec_code = IVEC_CODE;
  result -> array_size = size;
  result -> size = size;
  result -> iarr = am_malloc_ints(size);
  Ivecs_mallocked += 1;
  return(result);
}

void free_ivec(ivec *iv)
{
  check_ivec_code(iv,"free_ivec");
  iv -> ivec_code = 7777;

  am_free_ints(iv->iarr,iv->array_size);
  am_free((char *)iv,sizeof(ivec));

  Ivecs_freed += 1;
}

void fprintf_ivec(FILE *s,char *m1, const ivec *iv,char *m2)
{
  if ( iv == NULL )
    fprintf(s,"%s = (ivec *)NULL%s",m1,m2);
  else
  {
    char buff[100];
    check_ivec_code(iv,"fprintf_ivec");
    sprintf(buff,"%s = ",m1);
    fprintf_ints(s,buff,iv->iarr,iv->size,m2);
  }
}

int safe_ivec_ref(const ivec *iv, int i)
{
  check_ivec_access(iv,i,"ivec_ref");
  return(iv->iarr[i]);
}

void safe_ivec_set(ivec *iv,int i,int value)
{
  check_ivec_access(iv,i,"ivec_set");
  iv->iarr[i] = value;
}

void safe_ivec_increment(ivec *iv,int i,int value)
{
  check_ivec_access(iv,i,"ivec_increment");
  iv->iarr[i] += value;
}

void ivec_plus(ivec *d_1, ivec *d_2, ivec *r_d)
{
  int i;

  if ( d_1 -> size != d_2 -> size )
  {
    fprintf_ivec(stderr,"d_1",d_1,"\n");
    fprintf_ivec(stderr,"d_2",d_2,"\n");
    my_error("ivec_plus: ivecs have different shape");
  }

  for ( i = 0 ; i < r_d -> size ; i++ )
    r_d -> iarr[i] = d_1->iarr[i] + d_2 -> iarr[i];
}

ivec *mk_ivec_plus(ivec *a,ivec *b)
{
  ivec *result = mk_ivec(a->size);
  ivec_plus(a,b,result);
  return(result);
}

void ivec_subtract(ivec *d_1,ivec *d_2,ivec *r_d)
{
  int i;
  ivec *a = mk_ivec(d_2->size);
  for (i=0;i<a->size;i++) ivec_set(a,i,-ivec_ref(d_2,i));

  if ( d_1 -> size != d_2 -> size )
  {
    fprintf_ivec(stderr,"d_1",d_1,"\n");
    fprintf_ivec(stderr,"d_2",d_2,"\n");
    my_error("ivec_subtract: ivecs have different shape");
  }

  ivec_plus(d_1,a,r_d);
  free_ivec(a);
}

ivec *mk_ivec_subtract(ivec *a,ivec *b)
{
  ivec *result; 
  result = mk_ivec(a->size);
  ivec_subtract(a,b,result);
  return(result);
}

int safe_ivec_size(const ivec *iv)
{
  check_ivec_code(iv, "ivec_size");
  return(iv->size);
}

void copy_iarr_to_ivec(int *iarr,int size,ivec *r_iv)
{
  assert_ivec_shape(r_iv,size,"copy_iarr_to_ivec");
  copy_ints(iarr,r_iv->iarr,size);
}

ivec *mk_ivec_from_iarr(int *iarr,int size)
{
  ivec *result = mk_ivec(size);
  copy_iarr_to_ivec(iarr,size,result);
  return(result);
}

void copy_ivec_to_iarr(ivec *iv, int *iarr)
{
  check_ivec_code(iv,"copy_ivec_to_iarr");
  copy_ints(iv->iarr,iarr,iv->size);
}
  
int *mk_iarr_from_ivec(ivec *iv)
{
  int *result;
  check_ivec_code(iv,"make_copy_iarr");
  result = am_malloc_ints(iv->size);
  copy_ivec_to_iarr(iv,result);
  return(result);
}

/* Makes an ivec of size end - start:
   { start , start+1 , .... end-2 , end-1 } */
ivec *mk_sequence_ivec(int start_value,int end_value)
{
  int size = end_value - start_value;
  ivec *iv = mk_ivec(size);
  int i;
  for ( i = 0 ; i < size ; i++ )
    ivec_set(iv,i,start_value+i);
  return iv;
}

/*
   Allocates and returns an ivec of size size in which ivec[i] = i
   ivec[0] = 0
   ivec[1] = 1
    .
    .
   ivec[size-1] = size-1
*/
ivec *mk_identity_ivec(int size)
{
  return mk_sequence_ivec(0,size);
}

void shuffle_ivec(ivec *iv)
/*
   A random permutation of iv is returned
*/
{
  int size = ivec_size(iv);
  int i;
  for ( i = 0 ; i < size-1 ; i++ )
  {
    int j = int_random(size - i);
    if ( j > 0 )
    {
      int swap_me_1 = ivec_ref(iv,i);
      int swap_me_2 = ivec_ref(iv,i+j);
      ivec_set(iv,i,swap_me_2);
      ivec_set(iv,i+j,swap_me_1);
    }
  }
}

void constant_ivec(ivec *iv,int value)
{
  int i;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    ivec_set(iv,i,value);
}

ivec *mk_constant_ivec(int size,int value)
{
  ivec *iv = mk_ivec(size);
  constant_ivec(iv,value);
  return(iv);
}

void zero_ivec(ivec *iv)
{
  constant_ivec(iv,0);
}

ivec *mk_zero_ivec(int size)
{
  return(mk_constant_ivec(size,0));
}

void ivec_scalar_mult(ivec *iv,int scale,ivec *riv)
{
  int i;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    ivec_set(riv,i,ivec_ref(iv,i) * scale);
}

ivec *mk_ivec_scalar_mult(ivec *iv,int scale)
{
  ivec *result = mk_ivec(ivec_size(iv));
  ivec_scalar_mult(iv,scale,result);
  return(result);
}

void copy_ivec(ivec *iv,ivec *r_iv)
{
  ivec_scalar_mult(iv,1,r_iv);
}

ivec *mk_copy_ivec(ivec *iv)
{
  ivec *result = mk_ivec(ivec_size(iv));
  copy_ivec(iv,result);
  return(result);
}

void copy_ivec_subset(ivec *iv, int start, int size, ivec *r_iv)
{
  int i;
  for (i = start; i < start + size; i++)
    ivec_set(r_iv, i - start, ivec_ref(iv, i));
}

ivec *mk_copy_ivec_subset(ivec *iv, int start, int size)
{
  ivec *result = mk_ivec(size);
  copy_ivec_subset(iv, start, size, result);
  return (result);
}

int num_of_given_value(ivec *iv,int value)
{
  int i;
  int result = 0;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    if ( ivec_ref(iv,i) == value ) result += 1;

  return(result);
}

int num_zero_entries(ivec *iv)
{
  return(num_of_given_value(iv,0));
}

int num_nonzero_entries(ivec *iv)
{
  return(ivec_size(iv) - num_zero_entries(iv));
}

int ivec_extreme(const ivec *iv,bool lowest)
{
  int result;
  int i;

  if ( ivec_size(iv) <= 0 )
    my_error("Can't find min or max of empty ivec");
  
  result = ivec_ref(iv,0);

  for ( i = 1 ; i < ivec_size(iv) ; i++ )
  {
    int v = ivec_ref(iv,i);
    if ( lowest )
    {
      if ( v < result ) result = v;
    }
    else if ( v > result ) result = v;
  }

  return(result);
}

int ivec_min(const ivec *iv)
{
  return(ivec_extreme(iv,TRUE));
}

int ivec_max(const ivec *iv)
{
  return(ivec_extreme(iv,FALSE));
}

/* Sensible if args are NULL. False if different size */
bool ivec_equal(const ivec *x1,const ivec *x2)
{
  bool result = TRUE;

  if ( EQ_PTR(x1,x2) )
    result = TRUE;
  else if ( x1 == NULL || x2 == NULL )
    result = FALSE;
  else if ( ivec_size(x1) != ivec_size(x2) )
    result = FALSE;
  else
  {
    int i;
    for ( i = 0 ; result && i < ivec_size(x1) ; i++ ) 
      result = result && ivec_ref(x1,i) == ivec_ref(x2,i);
  }
  return(result);
}

ivec *mk_ivec_1(int x0)
{
  ivec *result = mk_ivec(1);
  ivec_set(result,0,x0);
  return(result);
}

ivec *mk_ivec_2(int x0,int x1)
{
  ivec *result = mk_ivec(2);
  ivec_set(result,0,x0);
  ivec_set(result,1,x1);
  return(result);
}

ivec *mk_ivec_3(int x0,int x1,int x2)
{
  ivec *result = mk_ivec(3);
  ivec_set(result,0,x0);
  ivec_set(result,1,x1);
  ivec_set(result,2,x2);
  return(result);
}

ivec *mk_ivec_4(int x0,int x1,int x2,int x3)
{
  ivec *result = mk_ivec(4);
  ivec_set(result,0,x0);
  ivec_set(result,1,x1);
  ivec_set(result,2,x2);
  ivec_set(result,3,x3);
  return(result);
}

ivec *mk_ivec_5(int x0,int x1,int x2,int x3,int x4)
{
  ivec *result = mk_ivec(5);
  ivec_set(result,0,x0);
  ivec_set(result,1,x1);
  ivec_set(result,2,x2);
  ivec_set(result,3,x3);
  ivec_set(result,4,x4);
  return(result);
}

ivec *mk_ivec_6(int x0,int x1,int x2,int x3,int x4,int x5)
{
  ivec *result = mk_ivec(6);
  ivec_set(result,0,x0);
  ivec_set(result,1,x1);
  ivec_set(result,2,x2);
  ivec_set(result,3,x3);
  ivec_set(result,4,x4);
  ivec_set(result,5,x5);
  return(result);
}

/**** Removal functions on ivecs ****/

/* Reduces the size of iv by one.
   ivec_ref(iv,index) disappears.
   Everything to the right of ivec_ref(iv,index) is copied one to the left.

   Formally: Let ivold be the ivec value beore calling this function
             Let ivnew be the ivec value after calling this function.

PRE: ivec_size(ivold) > 0
     0 <= index < ivec_size(ivold)

POST: ivec_size(ivnew) = ivec_size(ivold)-1
      for j = 0 , 1, 2 ... index-1  : 
         ivec_ref(ivnew,j) == ivec_ref(ivold,j)

      for j = index , index+1 , ... ivec_size(ivnew)-1:
         ivec_ref(ivnew,j) == ivec_ref(ivold,j+1)
*/
void ivec_remove(ivec *iv,int index)
{
  int i;
  int isize = ivec_size(iv);

#ifndef AMFAST
  if ( isize <= 0 ) my_error("ivec_remove: empty ivec");
  if ( index < 0 || index >= isize )
    my_error("ivec_remove: bad index");
#endif /* #ifndef AMFAST */

  for ( i = index ; i < isize - 1 ; i++ )
    ivec_set(iv,i,ivec_ref(iv,i+1));
  iv -> size -= 1;
}

/* Shrinks d by one element by removing the rightmost element. 
   Example:

  Before: d == ( 3 1 4 1 5 )
    ivec_remove_last_element(d)
  After d == ( 3 1 4 1 )
*/
void ivec_remove_last_element(ivec *iv)
{
#ifndef AMFAST
  if (ivec_size(iv) <= 0) my_error("ivec_remove_last_elt: empty ivec");
#endif /* #ifndef AMFAST */

  iv->size -= 1;
}

/* Remove all occurences of val in iv. Does nothing
   if val doesn't appear in ivec. Keeps remaining
   iv elements in the same order. */
void ivec_remove_value(ivec *iv,int val)
{
  int i = 0;
  while ( i < ivec_size(iv) )
  {
    if ( ivec_ref(iv,i) == val )
      ivec_remove(iv,i);
    else
      i += 1;
  }
}

int ivec_sum(ivec *iv)
{
  int result = 0;
  int i;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    result += ivec_ref(iv,i);
  return(result);
}

bool equal_ivecs(ivec *iv1,ivec *iv2)
{
  bool result = TRUE;
  int i;

  if ( ivec_size(iv1) != ivec_size(iv2) )
    my_error("qual_ivecs wrong sizes");

  for ( i = 0 ; i < ivec_size(iv1) && result ; i++ )
    result = ivec_ref(iv1,i) == ivec_ref(iv2,i);
  return(result);
}

void ivec_malloc_report()
{
  if ( Ivecs_mallocked > 0 )
  {
    fprintf(stdout,"# Integer Vectors  (datatype ivec) currently allocated: %d\n",
           Ivecs_mallocked - Ivecs_freed
          );
    if ( Ivecs_mallocked - Ivecs_freed != 0 )
    {
      fprintf(stdout,"#       Number of ivec allocations since program start:  %d\n",
             Ivecs_mallocked
            );
      fprintf(stdout,"#       Number of ivec frees       since program start:  %d\n#\n",
            Ivecs_freed
            );
    }
  }
}


/* Changed by Paul Komarek 19 May 1999 */
/* Increases length of ivec by 1 and puts val into
   new rightmost entry
*/
void add_to_ivec(ivec *iv,int new_val)
{
  int size;

  if ( iv->array_size < 0 || iv->size > iv->array_size )
    my_error("ivec size or array size has got muddles. Talk to AWM");

  size = iv->size;

  if ( size == iv->array_size )
  {
    /* Paul Komarek added the more conservative new_array_size on
       19 May 1999. */
    /* int new_array_size = 2 * iv->size + 2; */
    int s3 = 3 * size / 10;
    int new_array_size =
      2 + 11 * size / 10  +
      (size<10) * s3  +
      (size<100) * s3 +
      (size<500) * s3;

    int *iarr_new = AM_MALLOC_ARRAY(int,new_array_size);
    int i;
    for ( i = 0 ; i < size ; i++ )
      iarr_new[i] = iv->iarr[i];
    AM_FREE_ARRAY(iv->iarr,int,size);
    iv -> iarr = iarr_new;
    iv -> array_size = new_array_size;
  }
  iv->iarr[size] = new_val;
  iv -> size += 1;
}

     /* Increases iv in length by 1 and shifts all elements
        with original index greater or equal to index one to the
        right and inserts val at index. */
void ivec_insert(ivec *iv,int index,int val)
{
  int i;
  add_to_ivec(iv,0);
  for ( i = ivec_size(iv)-1 ; i > index ; i-- )
    ivec_set(iv,i,ivec_ref(iv,i-1));
  ivec_set(iv,index,val);
}


/* Find least index i such that value = ivec_ref(iv,i).
  If not found, returns -1
*/
int find_index_in_ivec(const ivec *iv, int value)
{
  int result = -1;
  int i;

  for ( i = 0 ; i < ivec_size(iv) && result < 0 ; i++ )
    if (value == ivec_ref(iv,i)) 
      result = i;
  return(result);
}

int find_index_in_ivec_hint(const ivec *iv, int value, int hint)
{
  int i;
  int sign = -1;
  int sz = ivec_size(iv);

  for (i=1; i<=sz; i++) {
	int disp = i/2;
	int idx = (hint + sign*disp) % sz;
    if (value == ivec_ref(iv, idx)) 
      return idx;
	sign *= -1;
  }
  return -1;
}

/* Finds leftmost index such that siv[index] > val. If none
   such, returns ivec_size(siv) */
int index_in_sorted_ivec(ivec *siv,int val)
{
  int result = -1;
  int i;
  for ( i = 0 ; i < ivec_size(siv) && result < 0 ; i++ )
    if ( ivec_ref(siv,i) > val )
      result = i;
  if ( result < 0 ) result = ivec_size(siv);
  return(result);
}

void add_to_sorted_ivec(ivec *siv, int val)
{
  int i;
  add_to_ivec(siv, 0);
  for (i = ivec_size(siv) - 2; i >= 0 && val < ivec_ref(siv, i); i--)
    ivec_set(siv, i + 1, ivec_ref(siv, i));
  ivec_set(siv, i + 1, val);
}

bool is_in_ivec(ivec *iv,int value)
{
  return(find_index_in_ivec(iv,value) >= 0);
}

void ivec_sort(ivec *dv,ivec *r_dv)
{
  check_ivec_code(dv,"ivec_sort (1st arg)");
  assert_ivec_shape(r_dv,dv->size,"ivec_sort");
  copy_ivec(dv,r_dv);
  sort_ints(r_dv->iarr,ivec_size(r_dv),r_dv->iarr);
}

ivec *mk_ivec_sort(ivec *iv)
{
  ivec *result = mk_ivec(ivec_size(iv));
  ivec_sort(iv,result);
  return result;
}

/*Could be made faster*/
ivec *mk_ivec_union(ivec *v1, ivec *v2)
{
  ivec *result = NULL;
  int i;
  if(v1 == NULL || v2 == NULL)
    my_error("Null arg to 'mk_ivec_union()'");

  result = mk_copy_ivec(v1);
  for(i = 0; i< ivec_size(v2); i++)
  {
    int x = ivec_ref(v2, i);

    if(!is_in_ivec(result, x))
      add_to_ivec(result, x);
  }

  return (result);
}

/*Pre: ivecs are ordered
*/

ivec *mk_ivec_union_ordered(ivec *v1, ivec *v2)
{
  ivec *result = NULL;
  int i1 = 0, i2 = 0;
  int num1 = 0, num2 = 0;

  if(v1 == NULL || v2 == NULL)
    my_error("Null arg to 'mk_ivec_union()'");

  result = mk_ivec(0);

  while(1)
  {
    if(i1 >  ivec_size(v1) - 1)
    {
      if(i2 > ivec_size(v2) -1)
	break;
      else
      {
	num1 = INT_MAX;
	num2 = ivec_ref(v2, i2);
      }
    }
    else if(i2 > ivec_size(v2)  - 1)
    {
      num2 = INT_MAX;
      num1 = ivec_ref(v1, i1);
    }
    else 
    {
      num1 = ivec_ref(v1, i1);
      num2 = ivec_ref(v2, i2);
    }
    
    if(num1 < num2)
    {
      add_to_ivec(result, num1);
      i1++;
    }
    else if(num2 < num1)
    {
      add_to_ivec(result, num2);
      i2++;
    }
    else
    {
      add_to_ivec(result, num1);
      i1++;
      i2++;
    }
  }
  return(result);
}
  

/*Pre ivecs are ordered*/
ivec *mk_ivec_diff_ordered(ivec *v1, ivec *v2)
{
  int i = 0;
  int j = 0;
  ivec *result = NULL;
  int	num1, num2 = INT_MAX;

  if(v1 == NULL)
    my_error("Null arg #1 to 'set_diff()'");

  if(v2 == NULL)
    result = mk_copy_ivec(v1);
  else
  {
    result = mk_ivec(0);
    while( i < ivec_size(v1) )
    {
		
      num1 = ivec_ref(v1, i);
      if( j < ivec_size(v2))
	num2 = ivec_ref(v2, j);
      else
	num2 = ivec_ref(v1, ivec_size(v1) -1 ) + 1;
	
      if(num1 > num2)
      {
	j++;
      }
      if(num1 < num2)
      {
	i++;
	add_to_ivec(result, num1);
      }
      if(num1 == num2)
      {
	i++;
	j++;
      }
    }	
  }
  return(result);
}	

ivec *mk_append_ivecs(ivec *a,ivec *b)
{
  int size_a = ivec_size(a);
  int size_b = ivec_size(b);
  ivec *c = mk_ivec(size_a + size_b);
  int i;
  for ( i = 0 ; i < size_a ; i++ )
    ivec_set(c,i,ivec_ref(a,i));
  for ( i = 0 ; i < size_b ; i++ )
    ivec_set(c,size_a + i,ivec_ref(b,i));
  return c;
}


/* And now, the same for arrays of strings. Note that strings internally
   are always copied around and dynamically alolocated and deallocated.
*/


#define STRING_ARRAY_CODE 20542

void check_string_array_code(string_array *sar,char *name)
{
  if ( sar == NULL )
  {
    fprintf(stderr,"NULL string_array passed in operation %s\n",name);
    my_error("string_array data structure");
  }
  if ( sar->string_array_code != STRING_ARRAY_CODE )
  {
    fprintf(stderr,"Attempt to access a non-allocated String Array\n");
    fprintf(stderr,"This is in the operation %s\n",name);
    my_error("string_array data structure error");
  }
}

void check_string_array_access(string_array *sar,int i, char *name)
{
  check_string_array_code(sar,name);

  if ( i < 0 || i >= sar->size )
  {
    fprintf(stderr,"In operation \"%s\"\n",name);
    fprintf(stderr,"the string_array has size = %d\n",sar->size);
    fprintf(stderr,"You tried to use index i=%d\n",i);
    fprintf(stderr,"Here is the string_array that was involved:\n");
    fprintf_string_array(stderr,"sarv",sar,"\n");
    my_error("check_string_array_access");
  }
}

void assert_string_array_shape(string_array *sar,int size,char *name)
{
  check_string_array_code(sar,name);

  if ( size != sar->size )
  {
    fprintf(stderr,"In operation \"%s\"\n",name);
    fprintf(stderr,"the string_array has size = %d\n", sar->size);
    fprintf(stderr,"But should have been predefined with the shape:\n");
    fprintf(stderr,"size = %d\n",size);
    my_error("assert_string_array_shape");
  }
}

int String_Arrays_mallocked = 0;
int String_Arrays_freed = 0;

string_array *mk_string_array(int size)
{
  string_array *result = AM_MALLOC(string_array);
  int i;

  result -> string_array_code = STRING_ARRAY_CODE;
  result -> size = size;
  result -> sarr_size = size;
  result -> sarr = AM_MALLOC_ARRAY(char_ptr,size);

  for ( i = 0 ; i < size ; i++ )
    result->sarr[i] = mk_copy_string("<UnDeFiNeD>");
  String_Arrays_mallocked += 1;
  return(result);
}

void free_string_array(string_array *sar)
{
  int i;
  check_string_array_code(sar,"free_string_array");
  sar -> string_array_code = 7777;

  for ( i = 0 ; i < sar->size ; i++ )
    if ( sar->sarr[i] != NULL )
      am_free_string(sar->sarr[i]);

  AM_FREE_ARRAY(sar->sarr,char_ptr,sar->sarr_size);
  AM_FREE(sar,string_array);

  String_Arrays_freed += 1;
}

string_array *mk_copy_string_array(string_array *sa)
{
  string_array *nsa = mk_string_array(string_array_size(sa));
  int i;
  for ( i = 0 ; i < string_array_size(sa) ; i++ )
    string_array_set(nsa,i,string_array_ref(sa,i));
  return(nsa);
}

string_array *mk_string_array_1(char *x)
{
  string_array *sa = mk_string_array(1);
  string_array_set(sa,0,x);
  return(sa);
}

string_array *mk_string_array_2(char *x1,char *x2)
{
  string_array *sa = mk_string_array(2);
  string_array_set(sa,0,x1);
  string_array_set(sa,1,x2);
  return(sa);
}

string_array *mk_string_array_3(char *x1,char *x2,char *x3)
{
  string_array *sa = mk_string_array(3);
  string_array_set(sa,0,x1);
  string_array_set(sa,1,x2);
  string_array_set(sa,2,x3);
  return(sa);
}

string_array *mk_string_array_4(char *x1,char *x2,char *x3,char *x4)
{
  string_array *sa = mk_string_array(4);
  string_array_set(sa,0,x1);
  string_array_set(sa,1,x2);
  string_array_set(sa,2,x3);
  string_array_set(sa,3,x4);
  return(sa);
}

string_array *mk_string_array_5(char *x1,char *x2,char *x3,char *x4,char *x5)
{
  string_array *sa = mk_string_array(5);
  string_array_set(sa,0,x1);
  string_array_set(sa,1,x2);
  string_array_set(sa,2,x3);
  string_array_set(sa,3,x4);
  string_array_set(sa,4,x5);
  return(sa);
}

void swap_string_array(string_array *sa, int i, int j)
{
	char *tmp = mk_copy_string(string_array_ref(sa,i));
	string_array_set(sa,i,string_array_ref(sa,j));
	string_array_set(sa,j,tmp);
  free_string(tmp);
}

void qsort_string_array(string_array *sa, int left, int right)
{
	int last, i;
	if(left>=right) return; /*basis*/
	swap_string_array(sa, left, (left+right)/2);
	last = left;
	for(i=left+1;i<=right;i++)
	{
		if(strcmp(string_array_ref(sa,i),string_array_ref(sa,left)) < 0)
		{
			swap_string_array(sa,++last,i);
		}
	}
		swap_string_array(sa, left, last);
	qsort_string_array(sa, left, last-1);
	qsort_string_array(sa,last+1,right);
}

/* Sort into lexicographic order */
void sort_string_array(string_array *sa)
{
  qsort_string_array(sa,0, string_array_size(sa)-1);
  
  /*the built in qsort didn't work for some reason.*/
  /*qsort((char *)sa->sarr,sa->size,sizeof(char *),strcmp);*/
}

string_array *mk_sort_string_array(string_array *sa)
{
  string_array *result = mk_copy_string_array(sa);
  sort_string_array(result);
  return result;
}

void fprintf_string_array(FILE *s,char *m1,string_array *sar,char *m2)
{
  int i;

  if (sar == NULL)
    fprintf(s,"%s = <NULL string array>%s",m1,m2);
  else
  {
    check_string_array_code(sar,"fprintf_string_array");
    if ( string_array_size(sar) == 0 )
      fprintf(s,"%s = <empty string array>%s",m1,m2);
    else
    {
      for ( i = 0 ; i < string_array_size(sar) ; i++ )
	fprintf(s,"%s[%3d] = %s%s",
		m1,i,
		(string_array_ref(sar,i)==NULL) ? "NULL" : string_array_ref(sar,i),
		m2
		);
    }
  }
}

void pstring_array(string_array *sa)
{
  fprintf_string_array(stdout,"sa",sa,"\n");
}


void fprintf_string_array_contents_on_one_line(FILE *s,string_array *sar)
{
  int i;
  check_string_array_code(sar,"fprintf_string_array_contents");

  for ( i = 0 ; i < string_array_size(sar) ; i++ )
    fprintf(s,"%s ",
            (string_array_ref(sar,i)==NULL) ? "NULL" : string_array_ref(sar,i)
           );
}

void fprintf_string_array_contents(FILE *s,string_array *sar)
{
  int i;
  check_string_array_code(sar,"fprintf_string_array_contents");

  for ( i = 0 ; i < string_array_size(sar) ; i++ )
    fprintf(s,"%s\n",
            (string_array_ref(sar,i)==NULL) ? "NULL" : string_array_ref(sar,i)
           );
}

char *string_array_ref(string_array *sar, int i)
{
  check_string_array_access(sar,i,"string_array_ref");
  return(sar->sarr[i]);
}

void string_array_set(string_array *sar,int i,char *value)
/* value is COPIED in */
{
  check_string_array_access(sar,i,"string_array_set");
  if ( sar->sarr[i] != NULL )
    am_free_string(sar->sarr[i]);
  sar->sarr[i] = mk_copy_string(value);
}

int string_array_size(string_array *sar)
{
  check_string_array_code(sar,"string_array_size");
  return(sar->size);
}

string_array *mk_string_array_from_array(char **sarr,int size)
{
  int i;
  string_array *result = mk_string_array(size);
  for ( i = 0 ; i < size ; i++ )
    string_array_set(result,i,sarr[i]);
  return(result);
}

/* Extends the length of the string_array by one (efficiently, if
   amortized over time) and adds "string" as the size-1'th entry where
   size is the new size.

   COPIES in string. It's okay if string is NULL.
*/
void add_to_string_array(string_array *sa,char *string)
{
  if ( sa -> size == sa -> sarr_size )
  {
    int new_sarr_size = 2 + 2 * sa->sarr_size;
    char **new_sarr = AM_MALLOC_ARRAY(char_ptr,new_sarr_size);
    int i;
    for ( i = 0 ; i < sa->size ; i++ )
      new_sarr[i] = sa->sarr[i];
    AM_FREE_ARRAY(sa->sarr,char_ptr,sa->sarr_size);
    sa -> sarr_size = new_sarr_size;
    sa -> sarr = new_sarr;
  }

  sa -> size += 1;
  sa -> sarr[sa->size-1] = (string==NULL) ? NULL : mk_copy_string(string);
}

void insert_in_string_array(string_array *sa, int pos, char *string)
{
  int i;
  if (pos > sa->size)
    my_error("Illegal string_array index.");
  i = sa->size - 1;
  add_to_string_array(sa, NULL);
  for (; i >= pos; i--)
    sa->sarr[i + 1] = sa->sarr[i];
  sa->sarr[pos] = (string == NULL) ? NULL : mk_copy_string(string);
}

/* Synonym for add_to_string_array. Not the preferred name, retained
   for compatibility
*/
void string_array_add(string_array *sa,char *string)
{
  add_to_string_array(sa,string);
}

void string_array_malloc_report()
{
  if ( String_Arrays_mallocked > 0 )
  {
    fprintf(stdout,"#         Number of string_arrays currently allocated: %d\n",
           String_Arrays_mallocked - String_Arrays_freed
          );
    if ( String_Arrays_mallocked - String_Arrays_freed != 0 )
    {
      fprintf(stdout,"#       Number of string_array allocs since prog start:  %d\n",
             String_Arrays_mallocked
            );
      fprintf(stdout,"#       Number of string_array frees  since prog start:  %d\n#\n",
            String_Arrays_freed
            );
    }
  }
}

/* Returns true if c is whitespace or containined in the 0-terminated
    sepper string. If seppers is NULL it's ignored.
*/
bool is_sepper(char c,char *seppers)
{
  bool result = c == ' ' || c == '\n' || c == '\r' || c == '\t';
  if ( !result && seppers != NULL )
  {
    int i;
    for ( i = 0 ; !result && seppers[i] != '\0' ; i++ )
      if ( seppers[i] == c ) result = TRUE;
  }
  return(result);
}

/**** Removal functions on string_arrays ****/

/* Reduces the size of sa by one.
   string_array_ref(sa,index) disappears.
   Everything to the right of string_array_ref(sa,index) is copied one to the left.

   Formally: Let saold be the string_array value beore calling this function
             Let sanew be the string_array value after calling this function.

PRE: string_array_size(saold) > 0
     0 <= index < string_array_size(saold)

POST: string_array_size(sanew) = string_array_size(saold)-1
      for j = 0 , 1, 2 ... index-1  : 
         string_array_ref(sanew,j) == string_array_ref(saold,j)

      for j = index , index+1 , ... string_array_size(sanew)-1:
         string_array_ref(sanew,j) == string_array_ref(saold,j+1)
*/
void string_array_remove(string_array *sa,int index)
{
  int i;
  if ( string_array_size(sa) <= 0 ) my_error("string_array_remove: empty string_array");
  if ( index < 0 || index >= string_array_size(sa) )
    my_error("string_array_remove: bad index");

  /* WKW - New version shuffles the pointers over by one instead
     of using string_array_set, which copies the (i+1)th element into the
     ith spot and then frees the (i+1)th element.  This is faster and
     it avoids some memory bugs in the old version */
  if( sa->sarr[index] != NULL )
  {
	  am_free_string(sa->sarr[index]);
	  sa->sarr[index]=NULL;
  }
  for ( i = index ; i < (sa->size - 1); i++ )
  {
	  sa->sarr[i] = sa->sarr[i+1];
  }

  /* Old Code - Removed by WKW
  for ( i = index ; i < string_array_size(sa) - 1 ; i++ )
    string_array_set(sa,i,string_array_ref(sa,i+1));
  */

  /* Old Code - Removed by WKW
  if ( sa->sarr[sa->size-1] != NULL ) 
    free_string(sa->sarr[sa->size-1]);
  */
  sa->sarr[sa->size-1] = NULL; /* Not strictly necessary but neater */

  /* END OF NEW CODE */

  sa -> size -= 1;
}

/* Shrinks sa by one element by removing the rightmost element. 
   Example:

  Before: sa == ( "I" "hate" "to" "see" "you" "leave" )
    string_array_remove_last_element(sa)
  Before: sa == ( "I" "hate" "to" "see" "you" )
*/
void string_array_remove_last_element(string_array *sa)
{
  string_array_remove(sa,string_array_size(sa)-1);
}


/* Stop characters are any whitespace character, or any
   character in the seppers string.

   If seppers is NULL it's ignored.

   Any stop character is removed and used as a separator breaking string
   into smaller chunks.
   Each word (chunk) is placed as a separate entry into the string array.

   If string is NULL or "" or has only sepper characters, a string_array of
   length 0 is returned.

   Modified 1-23-96 JS: A single string may include whitespace including
   a newline character if delimited by single or double quotes.  Note the
   delimiting characters are still left in the string, as are all backslashes
   used to produce the characters '"\ in the string.
   Modified 1-31-96 JS: The single quote is no longer considered a delimiter
*/
string_array *mk_broken_string_using_seppers(char *string,char *seppers)
{
  string_array *sa = mk_string_array(0);
  int i = 0;
  while ( string != NULL && string[i] != '\0' )
  {
    int j = 0;
    while ( is_sepper(string[i],seppers) )
      i += 1;

    if ( string[i] != '\0' )
    {
      int part_size,backslashcount = 0;
      char *part,stopchar = ' ';
      int k;

      while ( string[i+j] != '\0' )
      {
	if (stopchar == ' ')
	{
	  if (is_sepper(string[i+j],seppers)) break;
	  else if ((string[i+j]=='\"') && !(backslashcount%2)) stopchar = '\"';
	}
	else if (stopchar == '\"')
	{
	  /* bug fix? this used to say stopchar = '\n' which made it put the
             whole rest of the line into one string once it had seen a double
             quote.  Now it only includes up to the next double quote.
             8/24/99 JGS
	     */
	  if ((string[i+j] == '\"') && !(backslashcount %2)) stopchar = ' ';
	}

	if (string[i+j] == '\\') backslashcount++;
	else                     backslashcount=0;
	
        j++;
      }

      part_size = j+1;
      part = AM_MALLOC_ARRAY(char,part_size);

      for ( k = 0 ; k < j ; k++ )
        part[k] = string[i+k];
      if ( k != part_size-1 ) my_error("oaibxwibxpwibx");
      part[k] = '\0';
      string_array_add(sa,part);
      AM_FREE_ARRAY(part,char,part_size);
    }
    i = i+j;
  }
  return(sa);
}

/* whitespace is removed, and each word is placed as a separate entry
   int the string array.
   If string is NULL or "" or has only whitespace, a string_array of
   length 0 is returned.
*/
string_array *mk_broken_string(char *string)
{
  string_array *result = mk_broken_string_using_seppers(string,NULL);
  return(result);
}

/*  As above except removes all double quote marks before doing the splitting.
*/
string_array *mk_broken_quoteless_string(char *string)
{
  char *quoteless = mk_quoteless_string(string);
  string_array *result = mk_broken_string_using_seppers(quoteless,NULL);
  free_string(quoteless);
  return(result);
}

/* Removes all occurences of c from string. If input was a string of
   length m containing n 'c' characters, then after it's a string of
   length m-n. */
char *mk_string_without_character(char *string,char c)
{
  int string_num_chars = (int)strlen(string);
  int string_num_bytes = string_num_chars + 1;
  char *temp = AM_MALLOC_ARRAY(char,string_num_bytes);
  int si = 0;
  int ti = 0;
  char *result;

  for ( si = 0 ; si < string_num_chars ; si++ )
  {
    if ( string[si] != c )
    {
      temp[ti] = string[si];
      ti += 1;
    }
  }

  temp[ti] = '\0';
  result = mk_copy_string(temp);
  AM_FREE_ARRAY(temp,char,string_num_bytes);
  return result;
}

char *mk_quoteless_string(char *s)
{
  return mk_string_without_character(s,'\"');
}

/* In below, entries are separated by whitespace. If n == 0 returns "",
   else returns rightmost "n" entries in sa, concatenated with " "
   char between them.

   PRE: n <= size
*/
char *mk_string_from_last_n_with_separator(string_array *sa,int n,char sep)
{
  char *result;
  if ( n < 0 )
  {
    my_error("odcbowdn");
    result = NULL;
  }
  else if ( n == 0 )
    result = mk_copy_string("");
  else
  {
    int i;
    int array_size = (sep)? n : 1;
    int k = 0;
    int start = string_array_size(sa) - n;

/* array_size = sum of string lengths + n-1 seperators + 1 \0 character.
   (note wouldnt be true for n==0 but thats dealt with above)
*/
    for ( i = 0 ; i < n ; i++ )
      array_size += strlen(string_array_ref(sa,start+i));

    result = AM_MALLOC_ARRAY(char,array_size);

    for ( i = 0 ; i < n ; i++ )
    {
      int j;
      char *s = string_array_ref(sa,start+i);
      for ( j = 0 ; s[j] != '\0' ; j++ )
      {
        result[k] = s[j];
        k += 1;
      }
      if(i==n-1)
      {
        result[k] = '\0';
        k += 1;
      } 
      else if(sep) 
      {
        result[k] = sep;
        k += 1;
      }
    }

    if ( k != array_size ) my_error("lkaxcnowsdcbolub");
  }
  return(result);
}

char *mk_string_from_last_n(string_array *sa,int n){
  return mk_string_from_last_n_with_separator(sa,n,' ');
}

char *mk_string_from_string_array_with_separator(string_array *sa,char sep){
  return mk_string_from_last_n_with_separator(sa,string_array_size(sa),sep);
}

char *mk_string_from_string_array(string_array *sa){
  return mk_string_from_last_n_with_separator(sa,string_array_size(sa),' ');
}

/* Find least index i such that eq_string(string,string_array_ref(sa,i).
  If not found, returns -1
*/
int find_index_in_string_array(string_array *sa,char *string)
{
  int result = -1;
  int i;
  for ( i = 0 ; i < string_array_size(sa) && result < 0 ; i++ )
    if ( eq_string(string,string_array_ref(sa,i)) )
      result = i;
  return(result);
}

/* Find least index i such that:
      caseless_eq_string(string,string_array_ref(sa,i).
   If not found, returns -1
*/
int caseless_find_index_in_string_array(string_array *sa,char *string)
{
  int result = -1;
  int i;
  for ( i = 0 ; i < string_array_size(sa) && result < 0 ; i++ )
    if ( caseless_eq_string(string,string_array_ref(sa,i)) )
      result = i;
  return(result);
}

/* Returns TRUE if and only if there exists i such that
   eq_string(string_array_ref(sa,i),string)
*/
bool string_array_member(string_array *sa,char *string)
{
  return(find_index_in_string_array(sa,string)>=0);
}

/* PRE: 0 <= start_index <= end_index <= size
   Returns a new string array of size end_index - start_index
   consisting of 
     { sa[start_index] , sa[start_index+1] ... sa[end_index-1] }
*/
string_array *mk_string_array_segment(string_array *sa,
                                      int start_index,int end_index)
{
  int seg_size = end_index - start_index;
  string_array *seg = mk_string_array(seg_size);
  int i;
  for ( i = 0 ; i < seg_size ; i++ )
    string_array_set(seg,i,string_array_ref(sa,start_index+i));
  return seg;
}

string_matrix *mk_string_matrix(int rows,int cols)
{
  string_matrix *sm = AM_MALLOC(string_matrix);
  int i;
  sm -> array_size = rows;
  sm -> rows = rows;
  sm -> cols = cols;
  sm -> sas = AM_MALLOC_ARRAY(string_array_ptr,rows);
  for ( i = 0 ; i < rows ; i++ )
  {
    int j;
    sm->sas[i] = mk_string_array(cols);
    for ( j = 0 ; j < cols ; j++ )
      string_array_set(sm->sas[i],j,"");
  }
  return(sm);
}

int string_matrix_rows(string_matrix *sm)
{
  return(sm->rows);
}

int string_matrix_cols(string_matrix *sm)
{
  return(sm->cols);
}

char *string_matrix_ref(string_matrix *sm,int row,int col)
{
  if ( row < 0 || col < 0 || row >= sm->rows || col >= sm->cols )
    my_error("string_matrix_ref bad row and/or col");
  return(string_array_ref(sm->sas[row],col));
}

/* Makes a COPY of char *value */
void string_matrix_set(string_matrix *sm,int row,int col,char *value)
{
  if ( row < 0 || col < 0 || row >= sm->rows || col >= sm->cols )
    my_error("string_matrix_set bad row and/or col");
  string_array_set(sm->sas[row],col,value);
}

void fprintf_string_matrix(FILE *s,char *m1,string_matrix *sm,char *m2)
{
  if ( sm == NULL )
    fprintf(s,"%s = (string_matrix *)NULL%s",m1,m2);
  else if ( sm -> rows == 0 )
    fprintf(s,"%s = <string_matrix with 0 rows>%s",m1,m2);
  else
  {
    int i;
    for ( i = 0 ; i < string_matrix_rows(sm) ; i++ )
    {
      char buff[100];
      sprintf(buff,"%s[%2d]",m1,i);
      fprintf_string_array(s,buff,sm->sas[i],"\n");
    }
  }
}

void string_matrix_add_row(string_matrix *sm)
{
  int j;

  if ( sm->array_size < sm->rows )
    my_error("oujbcowlbucv");
  if ( sm->array_size == sm->rows )
  {
    int new_size = int_max(100,(int) (2.5 * sm->array_size));
    string_array **new_ar = AM_MALLOC_ARRAY(string_array_ptr,new_size);
    int i;
    for ( i = 0 ; i < new_size ; i++ )
      new_ar[i] = NULL;
    for ( i = 0 ; i < sm->rows ; i++ )
      new_ar[i] = sm->sas[i];

    am_free((char *)sm->sas,sm->array_size * sizeof(string_array_ptr));
    sm->sas = new_ar;
    sm->array_size = new_size;
  }
  sm->sas[sm->rows] = mk_string_array(sm->cols);
  
  /* Initialize new row with empty string entries. */
  for ( j = 0 ; j < sm->cols ; j++ )
    string_array_set(sm->sas[sm->rows],j,"");

  sm->rows += 1;
}


void free_string_matrix(string_matrix *sm)
{
  int i;
  for ( i = 0 ; i < string_matrix_rows(sm) ; i++ )
    free_string_array(sm->sas[i]);
  AM_FREE_ARRAY(sm->sas,string_array_ptr,sm->array_size);
  AM_FREE(sm,string_matrix);
}

int max_string_length_in_column(string_matrix *sm,int col)
{
  int result = 0;
  int row;
  for ( row = 0 ; row < string_matrix_rows(sm) ; row++ )
    result = int_max(result,(int)strlen(string_matrix_ref(sm,row,col)));
  return(result);
}

char *string_before_col(string_array *seps,int col)
{
  char *result;

  if (seps == NULL)
    result = (col == 0) ? "" : " ";
  else
    result = string_array_ref(seps,col);
  return(result);
}

char *rightmost_string(string_array *seps)
{
  char *result = (seps==NULL)?"":
                 string_array_ref(seps,string_array_size(seps)-1);
  return(result);
}

/* MALLOCS and returns a string array in which the i'th entry is
   the i'th row of a plain text representation of the contents
   of the string matrix. The columns of the table are separated by
   the strings in "seps". seps[0] always appears before the left-hand
   column. seps[1] before the second column. seps[cols] always
   appears to the right of the final column. The number of
   entries in "seps" must be one more than the number of columns
   in sm. 

   Alternatively, seps may be NULL in which case one space is
   printed between each column. No spaces to left or two right.
*/
string_array *mk_tabular_string_array(string_matrix *sm,string_array *seps)
{
  int chars_per_line = 0;
  char *buffer;
  int row,col;
  int rows = string_matrix_rows(sm);
  int cols = string_matrix_cols(sm);
  int buffer_size;
  string_array *result = mk_string_array(rows);

  for ( col = 0 ; col < cols ; col++ )
  {
    chars_per_line += max_string_length_in_column(sm,col);
    chars_per_line += strlen(string_before_col(seps,col));
  }
  chars_per_line += strlen(rightmost_string(seps));

  buffer_size = chars_per_line + 1;
  buffer = AM_MALLOC_ARRAY(char,buffer_size);

  for ( row = 0 ; row < rows ; row++ )
  {
    int i = 0;
    int j;
    char *last_sep = rightmost_string(seps);
    for ( col = 0 ; col < cols ; col++ )
    {
      char *sep = string_before_col(seps,col);
      char *con = string_matrix_ref(sm,row,col);
      int chars_after_con = max_string_length_in_column(sm,col) - 
                            strlen(con);
      if ( chars_after_con < 0 ) my_error("uwocbxwloubc");
      for ( j = 0 ; sep[j] != '\0' ; j++ )
      {
        buffer[i] = sep[j];
        i += 1;
      }
      for ( j = 0 ; con[j] != '\0' ; j++ )
      {
        buffer[i] = con[j];
        i += 1;
      }
      for ( j = 0 ; j < chars_after_con ; j++ )
      {
        buffer[i] = ' ';
        i += 1;
      }
    }
    for ( j = 0 ; last_sep[j] != '\0' ; j++ )
    {
      buffer[i] = last_sep[j];
      i += 1;
    }
    if ( i != chars_per_line ) my_error("owbxolwubdcoufcbv");
    buffer[i] = '\0';
    string_array_set(result,row,buffer);
  }

  AM_FREE_ARRAY(buffer,char,buffer_size);

  return(result);
}

void render_string_matrix(FILE *s,char *comment,string_matrix *sm)
{
  string_array *sa = mk_tabular_string_array(sm,NULL);
  int i;
  for ( i = 0 ; i < string_array_size(sa) ; i++ )
    fprintf(s,"%s%s\n",comment,string_array_ref(sa,i));
  free_string_array(sa);
}

void string_matrix_real_set(string_matrix *sm,int row,int col,double value)
{
  char buff[100];
  sprintf(buff,"%g",value);
  string_matrix_set(sm,row,col,buff);
}

/* On entry *ds_ref points to an array of characters of size *power_ref.
   Characters 0 .. (*pos_ref)-1 are filled with characters of a string
   that is being built up.

   IF (*pos_ref < *power_ref)
     then c is placed in the *pos_ref'th character and *pos_ref is incremented
   Else a new array of chars is made of twice the length of previous one,
   the previous array is copied in, the previous array is freed and *ds_ref
   is set to point to the new one. */
/* NULL is a valid string to pass to this function (and probably
   should be).  Use a size of zero for NULL. */
void add_char_to_dstring_s_p(char **ds_ref, char c,
			     int *pos_ref, int *power_ref)
{
  char *new_ds;
  int new_pow = *power_ref;

  if (*pos_ref >= *power_ref)
  {
    while (*pos_ref >= new_pow)
      if (*power_ref == 0)
	new_pow = 32;
      else
	new_pow += new_pow;
    new_ds = AM_MALLOC_ARRAY(char, new_pow);
    if (*ds_ref != NULL)
    {
      memcpy(new_ds, *ds_ref, *power_ref);
      AM_FREE_ARRAY(*ds_ref, char, *power_ref);
    }
    *ds_ref = new_ds;
    *power_ref = new_pow;
  }
  (*ds_ref)[(*pos_ref)++] = c;
}

/* Reads to the end of line, or end of file
   and adding all characters into string. If nothing on
   the line, returns a string of size zero. If nothing
   on line and file ends immediately returns NULL.

   Modified 1-23-96 JS: The reading always ends at an EOL (or EOF), but a
   single or double quote can be used to delimit a string that includes spaces
   or even newlines that should remain part of this string.  Therefore, if
   single quotes, double quotes, and backslashes are actually wanted in the 
   final string they should be specified by \' and \" and \\ respectively.
   Not the whole line is copied as is including all backslashes and quotes.
   All of the checking serves only to find the end of this line as defined
   by these rules.
   Modified 1-31-96 JS: single quote is no longer a valid delimiter
   Rewritten 6-28-96 MD: now resizes string to fit
*/

char *mk_string_from_line(FILE *s)
{
  int size = 0;
  int pos = 0;
  int c = '\0';
  int i = 0;
  char *line = NULL;
  char *retval;
  char *str = NULL;
  bool finished = FALSE;
  bool instring = FALSE;
  bool possible_quoted = FALSE;
#ifdef PC_MVIS_PLATFORM
  bool use_str = s == stdin;
#else
  bool use_str = FALSE;
#endif

  if ( use_str )
  {
    gui_prompt("Enter>",CLI_COMMAND_INPUT,FALSE);
    str = gui_getstring(s);
  }

  while (!finished){
    c = (use_str)? str[i++] : fgetc(s);
    switch (c)
    {
    case EOF:
      finished = TRUE;
      break;
    case '\n':
      if (instring)
      {
	      if (possible_quoted)
	      {
	        add_char_to_dstring_s_p(&line, '\\', &pos, &size);
	        possible_quoted = FALSE;
	      }
	      add_char_to_dstring_s_p(&line, (char) c, &pos, &size);
      }
      else
	      finished = TRUE;
      break;
    case '\\':
      if (possible_quoted)
      {
	      add_char_to_dstring_s_p(&line, (char) c, &pos, &size);
	      possible_quoted = FALSE;
      }
      else
	      possible_quoted = TRUE;
      break;
    case '\"':
      if (!possible_quoted)
	      instring = !instring;
      add_char_to_dstring_s_p(&line, (char) c, &pos, &size);
      break;
    default:
      if (possible_quoted)
      {
	      possible_quoted = FALSE;
	      add_char_to_dstring_s_p(&line, '\\', &pos, &size);
      }
      add_char_to_dstring_s_p(&line, (char) c, &pos, &size);
    }
  }
  if (line != NULL || c == '\n')
    add_char_to_dstring_s_p(&line, '\0', &pos, &size);

  /* Have to do this because the dstring in line may contain
     more characters than the string length. */
  retval = mk_copy_string(line);
  if (line != NULL)
    AM_FREE_ARRAY(line,char,size);
  
  if ( use_str ){
    gui_unprompt();
    free_string(str);
  }

  return (retval);
}

/* Reads to the end of line, or end of file, skipping white space
   and adding each new word into the string array. If nothing on
   the line, returns a string_array of size zero. If nothing
   on line and file ends immediately returns NULL
*/
string_array *mk_string_array_from_line(FILE *s)
{
  char *st = mk_string_from_line(s);
  string_array *sa;
  if ( st == NULL )
    sa = NULL;
  else
    sa = mk_broken_string(st);

  if ( st != NULL ) am_free_string(st);
  return(sa);
}

/* breaks row_string. Error unless the number of space-sparated substrings
   equals number of cols in string matrix. Sets sm(row,i) to i'th substring
   forall i
*/
void string_matrix_row_from_broken_string(
    string_matrix *sm,
    int row,
    char *row_string
  )
{
  string_array *sa = mk_broken_string(row_string);
  int i;
  if ( string_array_size(sa) != string_matrix_cols(sm) )
    my_error("sm from broken row: wrong size");
  for ( i = 0 ; i < string_array_size(sa) ; i++ )
    string_matrix_set(sm,row,i,string_array_ref(sa,i));
  free_string_array(sa);
}

/* PRE: ivecs must be same size.
   Returns true if and only if all components of dx are >= corresponding
   components of dy
*/
bool ivec_weakly_dominates(ivec *dx,ivec *dy)
{
  ivec *diff = mk_ivec_subtract(dx,dy);
  bool result = (ivec_size(dx) == 0) ? TRUE : ivec_min(diff) >= 0.0;
  free_ivec(diff);
  return(result);
}

string_array *mk_string_array_from_ivec(ivec *iv)
{
  string_array *sa;
  
  if ( ivec_size(iv) == 0 )
  {
    sa = mk_string_array(1);
    string_array_set(sa,0,"empty");
  }
  else
  {
    int i;
    sa = mk_string_array(ivec_size(iv));
    for ( i = 0 ; i < ivec_size(iv) ; i++ )
    {
      char buff[100];
      sprintf(buff,"%d",ivec_ref(iv,i));
      string_array_set(sa,i,buff);
    }
  }
  return(sa);
}

/* Makes a string of numbers, each separated by whitespace.
   No quotes or anything. Just numbers. */
char *mk_string_from_ivec(ivec *iv)
{
  string_array *sa = mk_string_array_from_ivec(iv);
  char *s = mk_string_from_string_array(sa);
  free_string_array(sa);
  return(s);
}

int find_index_in_sorted_string_array(string_array *sa,char *s)
{
  int result = -1;
  if ( string_array_size(sa) > 0 )
  {
    int i_lo = 0;
    int i_hi = string_array_size(sa)-1;
    char *s_lo = string_array_ref(sa,i_lo);
    char *s_hi = string_array_ref(sa,i_hi);
    if ( eq_string(s,s_lo) )
      result = i_lo;
    else if ( eq_string(s,s_hi) )
      result = i_hi;
    while ( result < 0 && i_hi > i_lo + 1 )
    {
      int i_mid = (i_lo + i_hi)/2;
      char *s_mid = string_array_ref(sa,i_mid);
      int cmp = strcmp(s,s_mid);
      /* cmp < 0 <=> s is lexicographically before s_mid */
      /* cmp = 0 <=> s equals s_mid */
      /* cmp > 0 <=> s is lexicographically after s_mid */
      if ( cmp == 0 )
        result = i_mid;
      else if ( cmp < 0 )
      {
        i_hi = i_mid;
        s_hi = s_mid;
      }
      else
      {
        i_lo = i_mid;
        s_lo = s_mid;
      }
    }
  }
  return(result);
}

/* This routine is painfully inefficient.  It actually recopies
   (read malloc and free ferociously) the strings all the way down
   the line to make an insertion.  This can be easily cured with 
   some pointer swapping instead.
*/
void insert_in_sorted_string_array(string_array *sa,char *s)
{
  int i = 0;
  int size = string_array_size(sa);
  while ( i < size && strcmp(string_array_ref(sa,i),s) < 0 )
    i += 1;
  if ( i == size )
    add_to_string_array(sa,s);
  else
  {
    int j;
    add_to_string_array(sa,"");
    for ( j = size-1 ; j >= i ; j-- )
      string_array_set(sa,j+1,string_array_ref(sa,j));
    string_array_set(sa,i,s);
  }
}

void maybe_insert_in_sorted_string_array(string_array *sa,char *s)
{
  if ( find_index_in_sorted_string_array(sa,s) < 0 )
    insert_in_sorted_string_array(sa,s);
}

/***** ivec_array time! ****/

#define INITIAL_ivec_ARRAY_SIZE 10

ivec_array *mk_empty_ivec_array()
{
  ivec_array *ivecarr = AM_MALLOC(ivec_array);
  ivecarr -> size = 0;
  ivecarr -> array_size = INITIAL_ivec_ARRAY_SIZE;
  ivecarr -> array = AM_MALLOC_ARRAY(ivec_ptr,ivecarr->array_size);
  return(ivecarr);
}

void add_to_ivec_array(ivec_array *ivecarr,ivec *this_ivec)
/*
   Assume ivec_array is previously of size n. After this it is of size
   n+1, and the n+1'th element is a COPY of this_ivec.
*/
{
  if ( ivecarr -> size == ivecarr -> array_size )
  {
    int new_size = 2 + 2 * ivecarr->array_size;
    ivec **new_array = AM_MALLOC_ARRAY(ivec_ptr,new_size);
    int i;
    for ( i = 0 ; i < ivecarr -> array_size ; i++ )
      new_array[i] = ivecarr->array[i];
    AM_FREE_ARRAY(ivecarr->array,ivec_ptr,ivecarr->array_size);
    ivecarr -> array = new_array;
    ivecarr -> array_size = new_size;
  }
  ivecarr->array[ivecarr->size] = (this_ivec==NULL) ? NULL : mk_copy_ivec(this_ivec);
  ivecarr->size += 1;
}

int ivec_array_size(ivec_array *ivecarr)
{
  return(ivecarr->size);
}

ivec *safe_ivec_array_ref(ivec_array *ivecarr,int index)
/*
     Returns a pointer (not a copy) to the index'th element stored in
   the ivec_array. Error if index < 0 or index >= size
*/
{
  ivec *result;
  if ( index < 0 || index >= ivec_array_size(ivecarr) )
  {
    result = NULL;
    my_error("ivec_array_ref");
  }
  else
    result = ivecarr->array[index];
  return(result);
}

void ivec_array_set(ivec_array *iva, int index, ivec *iv)
{
  if ((index < 0) || (iva == NULL) || (index >= iva->size))
        my_error("ivec_array_set: called with incompatible arguments");
  if (iva->array[index] != NULL)
        free_ivec(iva->array[index]);
  iva->array[index] = (iv == NULL) ? NULL : mk_copy_ivec(iv);
}

void fprintf_ivec_array(FILE *s,char *m1,ivec_array *ivecarr,char *m2)
{
  if ( ivec_array_size(ivecarr) == 0 )
    fprintf(s,"%s = <ivec_array with zero entries>%s",m1,m2);
  else
  {
    int i;
    for ( i = 0 ; i < ivec_array_size(ivecarr) ; i++ )
    {
      char buff[100];
      ivec *iv = ivec_array_ref(ivecarr,i);
      sprintf(buff,"%s[%2d]",m1,i);
      if ( iv == NULL )
        fprintf(s,"%s = NULL%s",buff,m2);
      else
        fprintf_ivec(s,buff,iv,m2);
    }
  }
}

void free_ivec_array(ivec_array *ivecarr)
{
  int i;
  for ( i = 0 ; i < ivec_array_size(ivecarr) ; i++ )
    if ( ivecarr->array[i] != NULL )
      free_ivec(ivecarr->array[i]);
  AM_FREE_ARRAY(ivecarr->array,ivec_ptr,ivecarr->array_size);
  AM_FREE(ivecarr,ivec_array);
}

ivec_array *mk_copy_ivec_array(ivec_array *ivecarr)
{
  ivec_array *new_ar = mk_empty_ivec_array();
  int i;

  for ( i = 0 ; i < ivec_array_size(ivecarr) ; i++ )
    add_to_ivec_array(new_ar,ivec_array_ref(ivecarr,i));

  return(new_ar);
}

ivec_array *mk_array_of_zero_length_ivecs(int size)
{
  ivec_array *iva = mk_empty_ivec_array();
  ivec *iv = mk_ivec(0);
  int i;

  for (i = 0; i < size; i++)
        add_to_ivec_array(iva, iv);
  free_ivec(iv);
  return(iva);
}

void ivec_array_remove(ivec_array *iva,int index)
{
  int i;
  ivec *iv = ivec_array_ref(iva,index);
  if ( iv != NULL ) free_ivec(iv);
  for ( i = index ; i < iva->size-1 ; i++ )
    iva->array[i] = iva->array[i+1];
  iva->array[iva->size-1] = NULL;
  iva->size -= 1;
}

// Added by Artur
//----------------------------------------------
ivec_array *mk_transpose_ivec_array( ivec_array *a )
{
	ivec_array *result  =NULL;
	if( a != NULL )
	{
		int i, j,plen,size=ivec_array_size(a);
		ivec *p;
		result = mk_array_of_zero_length_ivecs(size);
		for(i=0;i<size;i++)
		{
			p = ivec_array_ref(a,i);
			plen = ivec_size(p);
			for(j=0;j<plen;j++)
			{
				add_to_ivec( ivec_array_ref( result, ivec_ref( p, j ) ), i );
			}
		}
	}
	return(result);
}

/* This function is most useful for functions that want
   to accumulate a set of messages, error or otherwise, as they
   go.  This function also works when *err_mess is NULL initially.
 */
void add_to_error_message(char **errmess, char *new_mess)
{
  if (!new_mess) return;
  else
  {
    if (*errmess)
    {
      int len = strlen(*errmess) + strlen(new_mess) + 2;
      char *buf = am_malloc(len * sizeof(char));
      sprintf(buf,"%s%s",*errmess,new_mess);
      am_free_string(*errmess);
      *errmess = mk_copy_string(buf);
      am_free(buf,len*sizeof(char));
    }
    else *errmess = mk_copy_string(new_mess);
  }
}

/* same as above except new_mess is prepended */
void prepend_error_message(char **errmess, char *new_mess)
{
  if (!new_mess) return;
  else
  {
    if (*errmess)
    {
      int len = strlen(*errmess) + strlen(new_mess) + 2;
      char *buf = am_malloc(len * sizeof(char));
      sprintf(buf,"%s%s",new_mess,*errmess);
      am_free_string(*errmess);
      *errmess = mk_copy_string(buf);
      am_free(buf,len*sizeof(char));
    }
    else *errmess = mk_copy_string(new_mess);
  }
}


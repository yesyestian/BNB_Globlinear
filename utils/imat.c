/*
   File:        imat.c
   Author:      Andrew W. Moore
   Created:     Mon Dec 28 10:48:12 EST 1998
   Description: Matrix of integers

   Copyright 1998, Schenley Park Research, Inc
*/

#include "imat.h"

imat *mk_imat(int rows,int cols)
{
  imat *im = AM_MALLOC(imat);
  int i;
  im -> array_size = rows;
  im -> rows = rows;
  im -> cols = cols;
  im -> array = AM_MALLOC_ARRAY(int *,rows);
  for ( i = 0 ; i < rows ; i++ )
    im -> array[i] = AM_MALLOC_ARRAY(int,cols);
  return im;
}

void free_imat(imat *im)
{
  int i;
  for ( i = 0 ; i < im->rows ; i++ )
    AM_FREE_ARRAY(im->array[i],int,im->cols);
  AM_FREE_ARRAY(im->array,int *,im->array_size);
  AM_FREE(im,imat);
}

int slow_imat_ref(imat *im,int row,int col)
{
  if ( row < 0 || row >= im->rows ||
       col < 0 || col >= im->cols )
    my_error("imat_ref");
  return im->array[row][col];
}

void slow_imat_set(imat *im,int row,int col,int value)
{
  if ( row < 0 || row >= im->rows ||
       col < 0 || col >= im->cols )
    my_error("imat_set");
  im->array[row][col] = value;
}

void imat_increment(imat *im,int row,int col,int delta)
{
  if ( row < 0 || row >= im->rows ||
       col < 0 || col >= im->cols )
    my_error("imat_increment");
  im->array[row][col] += delta;
}

int slow_imat_rows(imat *im)
{
  return im->rows;
}

int slow_imat_cols(imat *im)
{
  return im->cols;
}

imat *mk_constant_imat(int rows,int cols,int value)
{
  int i,j;
  imat *im = mk_imat(rows,cols);
  for ( i = 0 ; i < rows ; i++ )
    for ( j = 0 ; j < cols ; j++ )
      imat_set(im,i,j,value);
  return im;
}

imat *mk_zero_imat(int rows,int cols)
{
  return mk_constant_imat(rows,cols,0);
}

void copy_ivec_to_imat_row(ivec *iv,imat *im,int row)
{
  int j;
  for ( j = 0 ; j < imat_cols(im) ; j++ )
    imat_set(im,row,j,ivec_ref(iv,j));
}

void imat_add_multiple_of_second(imat *im1,imat *im2,int scale,imat *result)
{
  int i,j;
  for ( i = 0 ; i < imat_rows(im1) ; i++ )
    for ( j = 0 ; j < imat_cols(im1) ; j++ )
      imat_set(result,i,j,imat_ref(im1,i,j) + scale * imat_ref(im2,i,j));
}

void imat_plus(imat *im1,imat *im2,imat *result)
{
  imat_add_multiple_of_second(im1,im2,1,result);
}

void imat_subtract(imat *im1,imat *im2,imat *result)
{
  imat_add_multiple_of_second(im1,im2,-1,result);
}

bool imat_equal(imat *im1,imat *im2)
{
  int i,j;
  bool result = TRUE;
  for ( i = 0 ; result && i < imat_rows(im1) ; i++ )
    for ( j = 0 ; result && j < imat_cols(im1) ; j++ )
      result = imat_ref(im1,i,j) == imat_ref(im2,i,j);
  return result;
}

int imat_sum(imat *im)
{
  int sum = 0;
  int i,j;
  for ( i = 0 ; i < imat_rows(im) ; i++ )
    for ( j = 0 ; j < imat_cols(im) ; j++ )
      sum += imat_ref(im,i,j);
  return sum;
}

imat *mk_imat_transpose(imat *im)
{
  int tim_rows = imat_cols(im);
  int tim_cols = imat_rows(im);
  imat *tim = mk_imat(tim_rows,tim_cols);
  int i,j;
  for ( i = 0 ; i < tim_rows ; i++ )
    for ( j = 0 ; j < tim_cols ; j++ )
      imat_set(tim,i,j,imat_ref(im,j,i));
  return tim;
}

ivec *mk_sum_of_imat_rows(imat *im)
{
  int cols = imat_cols(im);
  ivec *iv = mk_ivec(cols);
  int i;
  for ( i = 0 ; i < cols ; i++ )
  {
    ivec *col_i = mk_ivec_from_imat_col(im,i);
    ivec_set(iv,i,ivec_sum(col_i));
    free_ivec(col_i);
  }
  return iv;
}

ivec *mk_ivec_from_imat_row(imat *im,int row)
{
  int size = imat_cols(im);
  int i;
  ivec *iv = mk_ivec(size);
  for ( i = 0 ; i < size ; i++ )
    ivec_set(iv,i,imat_ref(im,row,i));
  return iv;
}

ivec *mk_ivec_from_imat_col(imat *im,int col)
{
  int size = imat_rows(im);
  int i;
  ivec *iv = mk_ivec(size);
  for ( i = 0 ; i < size ; i++ )
    ivec_set(iv,i,imat_ref(im,i,col));
  return iv;
}

dym *mk_dym_from_imat(imat *im)
{
  int rows = imat_rows(im);
  int cols = imat_cols(im);
  dym *d = mk_dym(rows,cols);
  int i,j;
  for ( i = 0 ; i < dym_rows(d) ; i++ )
    for ( j = 0 ; j < dym_cols(d) ; j++ )
      dym_set(d,i,j,(double)imat_ref(im,i,j));
  return d;
}

void fprintf_imat(FILE *s,char *m1,imat *im,char *m2)
{
  if ( im == NULL )
    fprintf(s,"%s = (imat *)NULL%s\n",m1,m2);
  else
  {
    dym *d = mk_dym_from_imat(im);
    fprintf_dym(s,m1,d,m2);
    free_dym(d);
  }
}

void pimat(imat *im)
{
  fprintf_imat(stdout,"im",im,"\n");
}

int imat_extreme(imat *im, bool max)
{
  int i,j;
  int res = -7777;
  if ( imat_rows(im)*imat_cols(im) < 1 ) my_error("imat_extreme: empty imat");
  for (i=0;i<im->rows;i++)
    for (j=0;j<im->cols;j++)
      if ((i == 0 && j == 0) || (max && imat_ref(im,i,j) > res)
	  || (!max && imat_ref(im,i,j) < res))
	res = imat_ref(im,i,j);
  return res;
}

int imat_min(imat *im)
{
  return(imat_extreme(im, FALSE));
}

int imat_max(imat *im)
{
  return(imat_extreme(im, TRUE));
}

void imat_add_row(imat *im)
{
  int i;
#ifndef AMFAST
  if ( im->array_size < im->rows ) my_error("imat_add_row");
#endif
  if ( im->array_size == im->rows )
  {
    int newsize = int_max(10,2 * im->rows);
    int **new_array = AM_MALLOC_ARRAY(int *,newsize);
    for ( i = 0 ; i < im->rows ; i++ )
      new_array[i] = im->array[i];
    AM_FREE_ARRAY(im->array,int *,im->rows);
    im->array = new_array;
    im->array_size = newsize;
  }
  im->array[im->rows] = AM_MALLOC_ARRAY(int,imat_cols(im));
  for ( i = 0 ; i < imat_cols(im) ; i++ )
    im->array[im->rows][i] = 0;
  im -> rows += 1;
}

void add_ivec_to_imat_as_new_row(imat *im,ivec *iv)
{
  imat_add_row(im);
  copy_ivec_to_imat_row(iv,im,imat_rows(im)-1);
}
 
void copy_imat(imat *src,imat *dest)
{
  int i,j;
  for ( i = 0 ; i < imat_rows(src) ; i++ )
    for ( j = 0 ; j < imat_cols(src) ; j++ )
      imat_set(dest,i,j,imat_ref(src,i,j));
}

imat *mk_copy_imat(imat *im)
{
  imat *result = mk_imat(imat_rows(im),imat_cols(im));
  copy_imat(im,result);
  return result;
}


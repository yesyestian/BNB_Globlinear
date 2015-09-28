/*
   File:        imat.h
   Author:      Andrew W. Moore
   Created:     Mon Dec 28 10:48:12 EST 1998
   Description: Header for Matrix of integers

   Copyright 1998, Schenley Park Research, Inc
*/


#ifndef IMAT_H
#define IMAT_H
#include "amut.h"

typedef struct imat
{
  int array_size;
  int rows;
  int cols;
  int **array;
} imat;

imat *mk_imat(int rows,int cols);
void free_imat(imat *im);
int slow_imat_ref(imat *im,int row,int col);
void slow_imat_set(imat *im,int row,int col,int value);
int slow_imat_rows(imat *im);
int slow_imat_cols(imat *im);

#ifdef AMFAST
#define imat_ref(im,i,j) ((im)->array[i][j])
#define imat_set(im,i,j,v) ((im)->array[i][j]=(v))
#define imat_rows(im) ((im)->rows)
#define imat_cols(im) ((im)->cols)
#else
#define imat_ref(im,i,j) slow_imat_ref(im,i,j)
#define imat_set(im,i,j,v) slow_imat_set(im,i,j,v)
#define imat_rows(im) slow_imat_rows(im)
#define imat_cols(im) slow_imat_cols(im)
#endif

imat *mk_constant_imat(int rows,int cols,int value);
imat *mk_zero_imat(int rows,int cols);
void copy_ivec_to_imat_row(ivec *iv,imat *im,int row);
void imat_plus(imat *im1,imat *im2,imat *result);
void imat_subtract(imat *im1,imat *im2,imat *result);
bool imat_equal(imat *im1,imat *im2);
int imat_sum(imat *im);
int imat_min(imat *im);
int imat_max(imat *im);
imat *mk_imat_transpose(imat *im);

void imat_increment(imat *im,int row,int col,int delta);
ivec *mk_sum_of_imat_rows(imat *im);
ivec *mk_ivec_from_imat_row(imat *im,int row);
ivec *mk_ivec_from_imat_col(imat *im,int col);

dym *mk_dym_from_imat(imat *im);
void fprintf_imat(FILE *s,char *m1,imat *im,char *m2);
void pimat(imat *im);
void add_ivec_to_imat_as_new_row(imat *im,ivec *iv);
void imat_add_row(imat *im);

void copy_imat(imat *src,imat *dest);
imat *mk_copy_imat(imat *im);

#endif /* #ifndef IMAT_H */

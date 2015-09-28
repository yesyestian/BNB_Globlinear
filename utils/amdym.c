
/*
   File:        amdym.c
   Author:      Andrew W. Moore
   Created:     Thu Sep 15 21:01:13 EDT 1994
   Updated:     amdm was split into amdyv, amdym and svd by Frank Dellaert, Aug 14 1997
   Description: Dynamically allocated and deallocated matrices

   Copyright 1996, Schenley Park Research
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "amdym.h"     /* matrices */
#include "amma.h"      /* Fast, non-fragmenting, memory management */
#include "amar.h"      /* Obvious operations on 1-d arrays */
#include "svd.h"       /* singular value decomposition */


#define DYM_CODE 4508

#define NOTHING_TO_DO

#ifdef AMFAST

#define check_dym_code(d,name) NOTHING_TO_DO

#else /* if AMFAST is not defined... */

void check_dym_code(const dym *d, const char *name)
{
  if ( d == NULL )
  {
    fprintf(stderr,"NULL dym passed in operation %s\n",name);
    my_error("dym data structure");
  }
  if ( d->dym_code != DYM_CODE )
  {
    fprintf(stderr,"Attempt to access a non-allocated DYnamic Matrix\n");
    fprintf(stderr,"This is in the operation %s\n",name);
    my_error("dym data structure error");
  }
}

#endif /* #ifdef AMFAST */

/*
* check_dym_access is only called in safe_dym_xxx
* It used to be non-functional with AMFAST,
* but Frank Dellaert reinstated it June 30 1997
*/

void check_dym_access(const dym *d,int i, int j, const char *name)
{
  check_dym_code(d,name); /* non-functional if AMFAST */

  if ( i < 0 || i >= d->rows || j < 0 || j >= d->cols )
  {
    fprintf(stderr,"In operation \"%s\"\n",name);
    fprintf(stderr,"the dym (dynamic matrix) has rows = %d, cols = %d\n",
            d->rows,d->cols
           );
    fprintf(stderr,"You tried to use indices i=%d j=%d\n",i,j);
    fprintf(stderr,"Here is the dym that was involved:\n");
    fprintf_dym(stderr,"dm",d,"\n");
    my_error("check_dym_access");
  }
}


#ifndef AMFAST

void assert_dym_shape(dym *d,int rows, int cols,char *name)
{
  check_dym_code(d,name);

  if ( rows != d->rows || cols != d->cols )
  {
    fprintf(stderr,"In operation \"%s\"\n",name);
    fprintf(stderr,"the dym (dynamic matrix) has rows = %d, cols = %d\n",
            d->rows,d->cols
           );
    fprintf(stderr,"But should have been predefined with the shape:\n");
    fprintf(stderr,"rows = %d, cols = %d\n",rows,cols);
    my_error("assert_dym_shape");
  }
}

#endif /* #ifdef AMFAST */

int Dyms_mallocked = 0;
int Dyms_freed = 0;

dym *mk_dym(int rows,int cols)
{
  dym *result = AM_MALLOC(dym);
  result -> dym_code = DYM_CODE;
  result -> rows = rows;
  result -> cols = cols;
  result -> rows_allocated = rows;
  result -> tdarr = am_malloc_2d_realnums(rows,cols);
  Dyms_mallocked += 1;
  return(result);
}

double dym_sum(dym *x)
{
  int i,j;
  double sum = 0.0;
  for ( i = 0 ; i < dym_rows(x) ; i++ )
    for ( j = 0 ; j < dym_cols(x) ; j++ )
      sum += dym_ref(x,i,j);
  return sum;
}

bool dym_is_ill_defined(dym *x)
{
  return is_ill_defined(dym_sum(x));
}

void free_dym(dym *d)
{
  int i;
  check_dym_code(d,"free_dym");
#ifndef AMFAST
  if ( dym_is_ill_defined(d) ) my_error("free_dym: dym contained NaN or Inf");
#endif

  for ( i = 0 ; i < d -> rows ; i++ )
  {
    am_free_realnums(d->tdarr[i],d->cols);
    d->tdarr[i] = NULL;
  }

  am_free((char *) (d->tdarr),sizeof(double_ptr) * d->rows_allocated);
  am_free((char *)d,sizeof(dym));

  Dyms_freed += 1;
}

void dym_malloc_report(void)
{
  dyv_malloc_report(); /* for backwards compatibility */
  if ( Dyms_mallocked )
  {
    fprintf(stdout,"# Dynamic Matrices (datatype dym) currently allocated:  %d\n",
           Dyms_mallocked - Dyms_freed
          );
    if ( Dyms_mallocked - Dyms_freed != 0 )
    {
      fprintf(stdout,"#       Number of dym allocations since program start:  %d\n",
             Dyms_mallocked
            );
      fprintf(stdout,"#       Number of dym frees       since program start:  %d\n#\n",
             Dyms_freed
            );
    }
  }
}

void add_row(dym *d)
{
  check_dym_code(d,"add_row");
  if ( d->rows_allocated < d->rows )
    my_error("oujbcowlbucv");
  if ( d->rows_allocated == d->rows )
  {
    int new_size = int_max(100,(int) (2.5 * d->rows_allocated));
    double **new_ar = AM_MALLOC_ARRAY(double_ptr,new_size);
    int i;
    for ( i = 0 ; i < new_size ; i++ )
      new_ar[i] = NULL;
    for ( i = 0 ; i < d->rows ; i++ )
      new_ar[i] = d->tdarr[i];

    am_free((char *)d->tdarr,d->rows_allocated * sizeof(double_ptr));
    d->tdarr = new_ar;
    d->rows_allocated = new_size;
  }
  d->tdarr[d->rows] = am_malloc_realnums(d->cols);
  set_realnums_constant(d->tdarr[d->rows],d->cols,0.0);
  d->rows += 1;
}

double safe_dym_ref(const dym *d, int i,int j)
{
  check_dym_access(d,i,j,"dym_ref");
  return(d->tdarr[i][j]);
}

void safe_dym_set(dym *d,int i,int j,double value)
{
  check_dym_access(d,i,j,"dym_set");
  d->tdarr[i][j] = value;
}

void safe_dym_increment(dym *d,int i,int j,double value)
{
  check_dym_access(d,i,j,"dym_increment");
  d->tdarr[i][j] += value;
}

void copy_dym_to_tdarr(dym *d,double **tdarr)
{
  copy_2d_realnums(d->tdarr,tdarr,d->rows,d->cols);
}
  
double **mk_tdarr_from_dym(dym *d)
{
  double **result;
  check_dym_code(d,"make_copy_tdarr");
  result = am_malloc_2d_realnums(d->rows,d->cols);
  copy_dym_to_tdarr(d,result);
  return(result);
}

void copy_tdarr_to_dym(double **tdarr,int rows,int cols,dym *r_d)
{
  assert_dym_shape(r_d,rows,cols,"copy_tdarr_to_dym");
  copy_2d_realnums(tdarr,r_d->tdarr,rows,cols);
}
  
dym *mk_dym_from_tdarr(double **tdarr,int rows,int cols)
{
  dym *result = mk_dym(rows,cols);
  copy_tdarr_to_dym(tdarr,rows,cols,result);
  return(result);
}

/*
 * Copying dym rows to dym rows
 * Frank Dellaert, Aug 14 1997
 */

void copy_dym_row_to_dym_row(dym *src,int src_row, dym *dst,int dst_row)
{
  int i, ncols;
  check_dym_code(src,"copy_dym_row_to_dym_row");
  check_dym_code(dst,"copy_dym_row_to_dym_row");
  if ( src_row < 0 || src_row >= src->rows )
    my_error("copy_dym_row_to_dym_row: illegal source row");
  if ( dst_row < 0 || dst_row >= dst->rows )
    my_error("copy_dym_row_to_dym_row: illegal destination row");
  if ( src->cols != dst->cols )
    my_error("copy_dym_row_to_dym_row: different number of columns");
  ncols = src->cols;
  for ( i = 0 ; i < ncols ; i++ )
    dst->tdarr[dst_row][i] = src->tdarr[src_row][i];
}

/***** Copying dyvs to and from rows and columns of dyms. And
       Making dyvs from rows and columns of dyms too. *********/

void copy_dyv_to_dym_row(dyv *dv,dym *dm,int row)
{
  int i;
  check_dym_code(dm,"copy_dyv_to_dym_row");
  if ( row < 0 || row >= dm->rows )
    my_error("copy_dyv_to_dym_row: illegal row");
  assert_dyv_shape(dv,dm->cols,"copy_dyv_to_dym_row");
  for ( i = 0 ; i < dm->cols ; i++ )
    dm->tdarr[row][i] = dv->farr[i];
}

void copy_dyv_to_dym_col(dyv *dv,dym *dm,int col)
{
  int i;
  check_dym_code(dm,"copy_dyv_to_dym_col");
  if ( col < 0 || col >= dm->cols )
    my_error("copy_dyv_to_dym_col: illegal col");
  assert_dyv_shape(dv,dm->rows,"copy_dyv_to_dym_col");
  for ( i = 0 ; i < dm->rows ; i++ )
    dm->tdarr[i][col] = dv->farr[i];
}

void copy_dym_row_to_dyv(const dym *dm,dyv *dv,int row)
{
  int i;
  check_dym_code(dm,"copy_dym_row_to_dyv");
  if ( row < 0 || row >= dm->rows )
    my_error("copy_dyv_to_dym_row: illegal row");
  assert_dyv_shape(dv,dm->cols,"copy_dyv_to_dym_row");
  for ( i = 0 ; i < dm->cols ; i++ )
    dv->farr[i] = dm->tdarr[row][i];
}

dyv *mk_dyv_from_dym_row(const dym *dm,int row)
{
  dyv *result = mk_dyv(dm->cols);
  copy_dym_row_to_dyv(dm,result,row);
  return(result);
}

void copy_dym_col_to_dyv(dym *dm,dyv *dv,int col)
{
  int i;
  check_dym_code(dm,"copy_dym_col_to_dyv");
  if ( col < 0 || col >= dm->cols )
    my_error("copy_dym_col_to_dyv: illegal col");
  assert_dyv_shape(dv,dm->rows,"copy_dym_col_to_dyv");
  for ( i = 0 ; i < dm->rows ; i++ )
    dv->farr[i] = dm->tdarr[i][col];
}

dyv *mk_dyv_from_dym_col(dym *dm,int col)
{
  dyv *result = mk_dyv(dm->rows);
  copy_dym_col_to_dyv(dm,result,col);
  return(result);
}

/***** Copying farrs to and from rows and columns of dyms. And
       Making farrs from rows and columns of dyms too. *********/

void copy_farr_to_dym_row(double *farr,dym *dm,int row)
{
  int i;
  check_dym_code(dm,"copy_farr_to_dym_row");
  if ( row < 0 || row >= dm->rows )
    my_error("copy_farr_to_dym_row: illegal row");
  for ( i = 0 ; i < dm->cols ; i++ )
    dm->tdarr[row][i] = farr[i];
}

void copy_farr_to_dym_col(double *farr,dym *dm,int col)
{
  int i;
  check_dym_code(dm,"copy_farr_to_dym_col");
  if ( col < 0 || col >= dm->cols )
    my_error("copy_farr_to_dym_col: illegal col");
  for ( i = 0 ; i < dm->rows ; i++ )
    dm->tdarr[i][col] = farr[i];
}

void copy_dym_row_to_farr(dym *dm,double *farr,int row)
{
  int i;
  check_dym_code(dm,"copy_dym_row_to_farr");
  if ( row < 0 || row >= dm->rows )
    my_error("copy_farr_to_dym_row: illegal row");
  for ( i = 0 ; i < dm->cols ; i++ )
    farr[i] = dm->tdarr[row][i];
}

double *mk_farr_from_dym_row(dym *dm,int row)
{
  double *result = am_malloc_realnums(dm->cols);
  copy_dym_row_to_farr(dm,result,row);
  return(result);
}

void copy_dym_col_to_farr(dym *dm,double *farr,int col)
{
  int i;
  check_dym_code(dm,"copy_dym_col_to_farr");
  if ( col < 0 || col >= dm->cols )
    my_error("copy_dym_col_to_farr: illegal col");
  for ( i = 0 ; i < dm->rows ; i++ )
    farr[i] = dm->tdarr[i][col];
}

double *mk_farr_from_dym_col(dym *dm,int col)
{
  double *result = am_malloc_realnums(dm->rows);
  copy_dym_col_to_farr(dm,result,col);
  return(result);
}

/**** Making whole dyms from dyvs ******/

dym *mk_col_dym_from_dyv(dyv *dv)
{
  dym *result;
  check_dyv_code(dv,"mk_col_dym_from_dyv");
  result = mk_dym(dv->size,1);
  copy_dyv_to_dym_col(dv,result,0);
  return(result);
}

dym *mk_row_dym_from_dyv(dyv *dv)
{
  dym *result;
  check_dyv_code(dv,"mk_row_dym_from_dyv");
  result = mk_dym(1,dv->size);
  copy_dyv_to_dym_row(dv,result,0);
  return(result);
}

dym *mk_diag_dym_from_dyv(dyv *dv)
{
  dym *result;
  int i;

  check_dyv_code(dv,"mk_row_dym_from_dyv");
  result = mk_dym(dv->size,dv->size);
  zero_dym(result);

  for ( i = 0 ; i < dv->size ; i++ )
    result->tdarr[i][i] = dv->farr[i];

  return(result);
}

void add_dyv_to_diag(dym *r_a,dyv *b)
{
  int i;

  check_dyv_code(b,"add_dyv_to_diag");
  assert_dym_shape(r_a,dyv_size(b),dyv_size(b),"add_dyv_to_diag");

  for ( i = 0 ; i < dyv_size(b) ; i++ )
    dym_increment(r_a,i,i,dyv_ref(b,i));
}

dym *mk_add_dyv_to_diag(dym *a,dyv *b)
{
  dym *result = mk_copy_dym(a);
  add_dyv_to_diag(result,b);
  return(result);
}
  
/**** Making whole dyms from farrs ******/

dym *mk_col_dym_from_farr(double *farr,int farr_size)
{
  dym *result;
  result = mk_dym(farr_size,1);
  copy_farr_to_dym_col(farr,result,0);
  return(result);
}

dym *mk_row_dym_from_farr(double *farr, int farr_size)
{
  dym *result;
  result = mk_dym(1,farr_size);
  copy_farr_to_dym_row(farr,result,0);
  return(result);
}

dym *mk_diag_dym_from_farr(double *farr, int farr_size)
{
  dym *result;
  int i;

  result = mk_dym(farr_size,farr_size);
  zero_dym(result);

  for ( i = 0 ; i < farr_size ; i++ )
    result->tdarr[i][i] = farr[i];

  return(result);
}

/***** Simple operations on dyms ******/

int dym_rows(const dym *d)
{
  check_dym_code(d,"dym_rows");
  return(d->rows);
}

int dym_cols(const dym *d)
{
  check_dym_code(d,"dym_cols");
  return(d->cols);
}

int safe_dyv_size(const dyv *d)
{
  check_dyv_code(d,"dyv_size");
  return(d->size);
}

#ifdef NEVER

dym *old_formatted_dym_from_lex(lex *lx,char *fstring,char fcode,char *fname)
{
  int cols = num_char_occurs(fstring,fcode);
  dym *d = mk_dym(0,cols);
  int i;

  for ( i = 0 ; i < lx -> num_lines ; i++ )
  {
    int line_len = lex_line_length(lx,i);
    if ( line_len > 0 && lex_line_ref(lx,i,0)->type == NUMBER )
    {
      int dym_col = 0;
      int lex_col;
      add_row(d);

      for ( lex_col = 0 ; lex_col < line_len ; lex_col++ )
      {
        if ( fstring[lex_col] == fcode )
        {
          lextok *lxt = lex_line_ref(lx,i,lex_col);
          if ( lxt->type != NUMBER )
          {
            fprintf(stderr,"Syntax error attempting to read a number from\n");
            fprintf(stderr,"file \"%s\". Line %d item %d should be number.\n",
                    fname,i+1,lex_col+1);
            my_error("Non-number in datafile");
          }
          dym_set(d,dym_rows(d)-1,dym_col,lxt->number);
          dym_col += 1;
        }
      }

      if ( dym_col != cols )
      {
        int j;
        fprintf(stderr,"Syntax error attempting to read numbers from\n");
        fprintf(stderr,"file \"%s\". Line %d should contain %d number%s.\n",
                fname,i+1,cols,(cols==1)?"":"s");
        fprintf(stderr,"There should be numbers in following line items:\n");
        for ( j = 0 ; fstring[j] != '\0' ; j++ )
          if ( fstring[j] == fcode ) fprintf(stderr," %d ",j+1);
        fprintf(stderr,"\n");
        my_error("Too short a line in datafile");
      }
    }
  }

  return(d);
}
#endif

#ifdef FANCY
char Amdm_error_string[1000];

bool basic_read_io_dyms(
    FILE *s,
    char *fname,
    char *format,
    bool last_col_is_output,
    dym **r_in_dym,
    dym **r_out_dym
  )
/*
   Returns TRUE if it succeeds.
   Returns FALSE if some kind of error, in which case the string 
     Amdm_error_string contains an error message. If returns an error,
     then allocates no memory, and all results are NULL.

   last_col_is_output is IGNORED except when format == NULL

   In this function format may be NULL or else a null-terminated string.

   If format is NULL and last_col_is_output is TRUE then you can treat
   the specification as though the function were called with format = iii..iio
   where the number of i's is one less than the number of numbers on the first
   numerical line of the file.
  
   If format is NULL and last_col_is_output is FALSE then you can treat
   the specification as though the function were called with format = iii..iii
   where the number of i's is the number of numbers on the first
   numerical line of the file.
  
   Let N be the number of characters in format.

   Then we read the file, assuming that every line which starts with a number
   as its first lexical item contains exactly N lexical items all numbers.
   Otherwise we'll signal a syntax error.

   What do we do with these numbers?

   The number on the i'th numerical row of the file, in the j'th numeric
   colum is either ignored, stored in *r_in_dym or stored in *r_out_dym.

   It is stored in *r_in_dym[i,k] if format[j] is the k'th 'i' character in
   format.

   It is stored in *r_out_dym[i,k] if format[j] is the k'th 'o' character in
   format.

   If format contains no i's , no dym is created, and *r_in_dym is set to NULL
   If format contains no o's , no dym is created, and *r_out_dym is set to NULL

    EXAMPLE:
    FILE STARTS HERE:
    .7 6 -9 2.1
    # Line starts with non-numeric so ignored
    Actually, this ignored too
    -1 3 4 3
    
If called with format = ii-o would produce

    *r_in_dym = [  0.7  6.0 ]    *r_out_dym = [ 2.1 ]
                [ -1.0  3.0 ]                 [ 3.0 ]

If called with format = --i- would produce

    *r_in_dym = [  -9.0 ]    *r_out_dym = NULL
                [   4.0 ]                 

If called with format = o-io would produce

    *r_in_dym = [  -9.0 ]    *r_out_dym = [  0.7 2.1 ]
                [   4.0 ]                 [ -1.0 3.0 ]

If called with format = iiio would produce

    *r_in_dym = [  0.7 6.0 -9.0 ]    *r_out_dym = [ 2.1 ]
                [ -1.0 3.0  4.0 ]                 [ 3.0 ]

If called with format = NULL and last_col_is_output == TRUE would also produce

    *r_in_dym = [  0.7 6.0 -9.0 ]    *r_out_dym = [ 2.1 ]
                [ -1.0 3.0  4.0 ]                 [ 3.0 ]

If called with format = iiii would produce

    *r_in_dym = [  0.7 6.0 -9.0 2.1 ]    *r_out_dym = NULL
                [ -1.0 3.0  4.0 3.0 ]                 

If called with format = NULL and last_col_is_output == FALSE would also produce

    *r_in_dym = [  0.7 6.0 -9.0 2.1 ]    *r_out_dym = NULL
                [ -1.0 3.0  4.0 3.0 ]                 


If called with format = o-ioo would produce ERROR because numeric lines
don't contain 5 numbers.

*/
{
  lex lx[1];
  bool file_ended = FALSE;
  dym *din = NULL;
  dym *dout = NULL;
  int line_number = 0;
  int din_cols = (format==NULL) ? 0 : num_char_occurs(format,'i');
  int dout_cols = (format==NULL) ? 0 : num_char_occurs(format,'o');
  int total_cols = (format==NULL) ? 0 : strlen(format);
  bool ok = TRUE;

  sprintf(Amdm_error_string,"Everything's okay");

  while ( ok && !file_ended )
  {
    file_ended = !lex_can_read_line(s,lx);
    line_number += 1;

    if ( !file_ended )
    {
      int line_len = lex_line_length(lx,0);

      if ( (line_number % 1000)==0 ) basic_message(fname,line_number);
 
      if ( line_len > 0 && lex_line_ref(lx,0,0)->type == NUMBER )
      {
        int din_col = 0;
        int dout_col = 0;
        int lex_col;

        if ( format == NULL )
        {
          dout_cols = (last_col_is_output) ? 1 : 0;
          din_cols = line_len - dout_cols;
          total_cols = din_cols + dout_cols;
        }

        if ( din == NULL && din_cols > 0 ) din = mk_dym(0,din_cols);
        if ( dout == NULL && dout_cols > 0 ) dout = mk_dym(0,dout_cols);

        if ( din_cols > 0 ) add_row(din);
        if ( dout_cols > 0 ) add_row(dout);

        for ( lex_col = 0 ; ok && lex_col < int_min(line_len,total_cols) ;
              lex_col++ 
            )
        {
          bool is_input,is_output;

          if ( format != NULL )
            is_input = format[lex_col] == 'i';
          else
            is_input = lex_col < din_cols;

          if ( format != NULL )
            is_output = format[lex_col] == 'o';
          else
            is_output = lex_col == din_cols;

          if ( is_input || is_output )
          {
            lextok *lxt = lex_line_ref(lx,0,lex_col);
            if ( lxt->type != NUMBER )
            {
              ok = FALSE;
              sprintf(Amdm_error_string,
                      "Line %d item %d of file %s unexpected non-number",
                      line_number,lex_col+1,fname
                     );
            }
            else if ( is_input )
            {
              dym_set(din,dym_rows(din)-1,din_col,lxt->number);
              din_col += 1;
            }
            else
            {
              dym_set(dout,dym_rows(dout)-1,dout_col,lxt->number);
              dout_col += 1;
            }
          }
        }

        if ( ok && line_len != total_cols )
        {
          sprintf(Amdm_error_string,
                  "Line %d of file %s should have %d item(s). But it has %d",
                   line_number,fname,total_cols,line_len
                 );
          ok = FALSE;
        }
      }
    }
    lex_free_contents(lx);
  }

  if ( ok && din == NULL && dout == NULL )
  {
    if ( format == NULL )
    {
      sprintf(Amdm_error_string,
              "File %s has no lines with numbers, and there's no format given",
              fname
             );
      ok = FALSE;
    }
    else
    {
      if ( din_cols > 0 ) din = mk_dym(0,din_cols);
      if ( dout_cols > 0 ) dout = mk_dym(0,dout_cols);
    }
  }

  if ( ok && line_number > 1000 ) basic_message(fname,line_number);

  if ( ok )
  {
    *r_in_dym = din;
    *r_out_dym = dout;
  }
  else
  {
    if ( din != NULL ) free_dym(din);
    if ( dout != NULL ) free_dym(dout);
    *r_in_dym = NULL;
    *r_out_dym = NULL;
  }

  return(ok);
}

/* JS 9-7-95 
 * read_io_dyms is changed to allow MATLAB to call our software and pass
 * its matrices in directly.  The new parameter list is the same except
 * for the elimination of the file pointer.  If MATLAB did call this execution
 * and it passed in a matrix directly, fname will now contain 8 characters
 * to indicate it, the next set of 4 characters will be the integer number of
 * matrix rows, the next 4 will be the integer number of colums, and the next
 * 4 will be a double pointer to the values.  The old form of read_io_dyms 
 * exists as read_io_dyms_stream.  The new one should be used unless the caller
 * needs to pre-open the file and re-position the pointer in it before reaching
 * these routines.
 */
/*
   Returns TRUE if it succeeds.
   Returns FALSE if some kind of error, in which case the string 
     Amdm_error_string contains an error message. If returns an error,
     then allocates no memory, and all results are NULL.
*/
int MATLAB;
bool can_read_io_dyms(char *fname,char *format,dym **r_in_dym,dym **r_out_dym)
{
  int i,ocols,icols,mcols,mrows,icount,ocount;
  double *vals;
  char *lformat;
  bool ok = TRUE;

  if (MATLAB&&(!strcmp(fname,"MATLAB"))&&(fname[7]==0x15))
  {
    mrows = *((int *)(fname+8));
    mcols = *((int *)(fname+12));
    vals = *((double **)(fname+16));
    if (format) 
    {
      if ((int) strlen(format) != mcols) 
        my_error("read_io_dyms: format length does not match matrix size");
      lformat = format;
    }
    else
    {
      lformat = am_malloc((mcols+1)*sizeof(char));
      for (i=0;i<mcols-1;i++) lformat[i] = 'i';
      lformat[mcols-1] = 'o';
      lformat[mcols] = '\0';
    }
    for (i=0,icols=0,ocols=0;i< (int) strlen(lformat);i++)
    {
      if (lformat[i] == 'i') icols++;
      if (lformat[i] == 'o') ocols++;
    }
    if (icols) (*r_in_dym) = mk_dym(mrows,icols);
    else       (*r_in_dym) = NULL;
    if (ocols) (*r_out_dym) = mk_dym(mrows,ocols);
    else       (*r_out_dym) = NULL;
    for (i=0,icount=0,ocount=0;i< (int) strlen(lformat);i++)
    {
      if (lformat[i]=='i')
        copy_farr_to_dym_col(vals+(i*mrows),*r_in_dym,icount++);
      if (lformat[i]=='o')
        copy_farr_to_dym_col(vals+(i*mrows),*r_out_dym,ocount++);
    }
    if (!format) am_free(lformat,(mcols+1)*sizeof(char));
  }
  else
  {
    FILE *s = fopen(fname,"r");
    if ( s == NULL )
    {
      sprintf(Amdm_error_string,"Dataset file %s doesn't exist",fname);
      ok = FALSE;
      *r_in_dym = NULL;
      *r_out_dym = NULL;
    }
    else
      ok = read_io_dyms_stream(s,fname,format,r_in_dym,r_out_dym);
  }
  return(ok);
}

/* Just the same as read_io_dyms, except halts program if an error */
void read_io_dyms(char *fname,char *format,dym **r_in_dym,dym **r_out_dym)
{
  bool ok = can_read_io_dyms(fname,format,r_in_dym,r_out_dym);
  if ( !ok )
  {
    fprintf(stderr,"Error reading in dym(s):\n");
    my_error(Amdm_error_string);
  }
}

/*
bool read_io_dyms_stream(
    FILE *s,
    char *fname,
    char *format,
    dym **r_in_dym,
    dym **r_out_dym
  )

   Returns TRUE if it succeeds.
   Returns FALSE if some kind of error, in which case the string 
     Amdm_error_string contains an error message. If returns an error,
     then allocates no memory, and all results are NULL.

   In this function format may be NULL or else a null-terminated string.

   If format is NULL  then you can treat
   the specification as though the function were called with format = iii..iio
   where the number of i's is one less than the number of numbers on the first
   numerical line of the file.
  
   Let N be the number of characters in format.

   Then we read the file, assuming that every line which starts with a number
   as its first lexical item contains exactly N lexical items all numbers.
   Otherwise we'll signal a syntax error.

   What do we do with these numbers?

   The number on the i'th numerical row of the file, in the j'th numeric
   colum is either ignored, stored in *r_in_dym or stored in *r_out_dym.

   It is stored in *r_in_dym[i,k] if format[j] is the k'th 'i' character in
   format.

   It is stored in *r_out_dym[i,k] if format[j] is the k'th 'o' character in
   format.

   *The following is different from basic_read_io_dyms*

   If format contains no o's and no i's, *r_in_dym and *r_out_dym are
   both set to dyms of size 0x0.

   If format contains o's but no i's, *r_in_dym is set to a matrix with
   0 columns, but the same number of rows as *r_out_dym

   If format contains i's but no o's, *r_out_dym is set to a matrix with
   0 columns, but the same number of rows as *r_in_dym

    EXAMPLE:
    FILE STARTS HERE:
    .7 6 -9 2.1
    # Line starts with non-numeric so ignored
    Actually, this ignored too
    -1 3 4 3
    
If called with format = ii-o would produce

    *r_in_dym = [  0.7  6.0 ]    *r_out_dym = [ 2.1 ]
                [ -1.0  3.0 ]                 [ 3.0 ]

If called with format = --i- would produce

    *r_in_dym = [  -9.0 ]    *r_out_dym = [ ]
                [   4.0 ]                 [ ]

If called with format = o-io would produce

    *r_in_dym = [  -9.0 ]    *r_out_dym = [  0.7 2.1 ]
                [   4.0 ]                 [ -1.0 3.0 ]

If called with format = iiio would produce

    *r_in_dym = [  0.7 6.0 -9.0 ]    *r_out_dym = [ 2.1 ]
                [ -1.0 3.0  4.0 ]                 [ 3.0 ]

If called with format = NULL would also produce

    *r_in_dym = [  0.7 6.0 -9.0 ]    *r_out_dym = [ 2.1 ]
                [ -1.0 3.0  4.0 ]                 [ 3.0 ]

If called with format = o-ioo would produce ERROR because numeric lines
don't contain 5 numbers.

*/
bool read_io_dyms_stream(
    FILE *s,
    char *fname,
    char *format,
    dym **r_in_dym,
    dym **r_out_dym
  )
{
  bool last_col_is_output = TRUE;
  bool ok = 
    basic_read_io_dyms(s,fname,format,last_col_is_output,r_in_dym,r_out_dym);

  if ( ok && *r_in_dym == NULL && *r_out_dym == NULL )
  {
    *r_in_dym = mk_dym(0,0);
    *r_out_dym = mk_dym(0,0);
  }
  else if ( ok && *r_out_dym == NULL )
    *r_out_dym = mk_dym(dym_rows(*r_in_dym),0);
  else if ( ok && *r_in_dym == NULL )
    *r_in_dym = mk_dym(dym_rows(*r_out_dym),0);

  return(ok);
}

dym *read_dym(FILE *s,char *fname,char *format)
/*
   In this function format may be NULL or else a null-terminated string.

   If format is NULL then you can treat
   the specification as though the function were called with format = iii..iii
   where the number of i's is the number of numbers on the first
   numerical line of the file.
  
   Let N be the number of characters in format.

   Then we read the file, assuming that every line which starts with a number
   as its first lexical item contains exactly N lexical items all numbers.
   Otherwise we'll signal a syntax error.

   What do we do with these numbers?

   The number on the i'th numerical row of the file, in the j'th numeric
   colum is either ignored, stored in dym *result (the dym we make and return)

   It is stored in result[i,k] if format[j] is the k'th 'i' character in
   format.

   If format contains no i's , no dym is created, and *r_in_dym is set to NULL

    EXAMPLE:
    FILE STARTS HERE:
    .7 6 -9 2.1
    # Line starts with non-numeric so ignored
    Actually, this ignored too
    -1 3 4 3
    
If called with format = ii-- would produce

       result = [  0.7  6.0 ]
                [ -1.0  3.0 ]

If called with format = --i- would produce

       result = [  -9.0 ]    
                [   4.0 ]                 

If called with format = --i- would produce

       result = [  -9.0 ]
                [   4.0 ]

If called with format = iii- would produce

       result = [  0.7 6.0 -9.0 ]
                [ -1.0 3.0  4.0 ]

If called with format = iiii would produce

       result = [  0.7 6.0 -9.0 2.1 ]
                [ -1.0 3.0  4.0 3.0 ]                 

If called with format = NULL  would also produce

       result = [  0.7 6.0 -9.0 2.1 ]
                [ -1.0 3.0  4.0 3.0 ]                 


If called with format = o-ioo would produce ERROR because numeric lines
don't contain 5 numbers.

*/
{
  int i,mrows,mcols,icols,icount;
  double *vals;
  dym *din,*dout;
  bool last_col_is_output = FALSE;
  bool ok;
  char *lformat;

  if (MATLAB&&(!strcmp(fname,"MATLAB"))&&(fname[7]==0x15))
  {
    mrows = *((int *)(fname+8));
    mcols = *((int *)(fname+12));
    vals = *((double **)(fname+16));
    if (format) 
    {
      if ((int) strlen(format) != mcols) 
        my_error("read_io_dyms: format length does not match matrix size");
      lformat = format;
    }
    else
    {
      lformat = am_malloc((mcols+1)*sizeof(char));
      for (i=0;i<mcols;i++) lformat[i] = 'i';
      lformat[mcols] = '\0';
    }
    for (i=0,icols=0;i< (int) strlen(lformat);i++)
      if (lformat[i] == 'i') icols++;
    if (icols) din = mk_dym(mrows,icols);
    else       din = NULL;
    for (i=0,icount=0;i< (int) strlen(lformat);i++)
      if (lformat[i]=='i')
        copy_farr_to_dym_col(vals+(i*mrows),din,icount++);
    if (!format) am_free(lformat,(mcols+1)*sizeof(char));
  }
  else
  {
    ok = basic_read_io_dyms(s,fname,format,last_col_is_output,&din,&dout);  
    if ( !ok )
    {
      fprintf(stderr,"read_dym: Error reading in dym:\n");
      my_error(Amdm_error_string);
    }
    
    if ( dout != NULL ) free_dym(dout);
    /* Need to do this if someone gave us a format with 'o's in it. We
       are meant to ignore outputs in this function.
       */
  }
  return(din);
}
#endif

void save_dym(FILE *s,dym *d)
{
  int i,j;
  for ( i = 0 ; i < dym_rows(d) ; i++ )
    for ( j = 0 ; j < dym_cols(d) ; j++ )
      fprintf(s,"%12g%s",dym_ref(d,i,j),
              (j==dym_cols(d)-1) ? "\n" : " "
             );
}

void save_io_dyms(FILE *s,dym *ins,dym *outs)
{
  int i,j;
  if ( dym_rows(ins) != dym_rows(outs) )
    my_error("owslbcdoacpbifb");

  for ( i = 0 ; i < dym_rows(ins) ; i++ )
  {
    for ( j = 0 ; j < dym_cols(ins) ; j++ )
      fprintf(s,"%12g%s",dym_ref(ins,i,j),
              (j==dym_cols(ins)-1) ? "    " : " "
             );
    for ( j = 0 ; j < dym_cols(outs) ; j++ )
      fprintf(s,"%12g%s",dym_ref(outs,i,j),
              (j==dym_cols(outs)-1) ? "\n" : " "
             );
  }
}

#define DYM_SVD_NULL_THRESH (1e-18)

/* NOTE:  ***SEE DYM NOTES IN AMDYM.H ***** */

void constant_dym(dym *r_d,double v)
{
  int i,j;
  check_dym_code(r_d,"constant_dym");
  for ( i = 0 ; i < r_d -> rows ; i++ )
    for ( j = 0 ; j < r_d -> cols ; j++ )
      r_d->tdarr[i][j] = v;
}

void zero_dym(dym *r_d)
{
  
  check_dym_code(r_d,"zero_dym");
  constant_dym(r_d,0.0);
}

dym *mk_constant_dym(int rows,int cols,double v)
{
  dym *result = mk_dym(rows,cols);
  constant_dym(result,v);
  return(result);
}

dym *mk_zero_dym(int rows,int cols)
{
  dym *result = mk_dym(rows,cols);
  zero_dym(result);
  return(result);
}

/**** Standard operations on dyms ****/

void dym_scalar_mult(dym *d, double alpha, dym *r_d)
{
  int i,j;
  assert_dym_shape(r_d,d->rows,d->cols,"dym_scalar_mult");
  for ( i = 0 ; i < r_d -> rows ; i++ )
    for ( j = 0 ; j < r_d -> cols ; j++ )
      r_d -> tdarr[i][j] = d->tdarr[i][j] * alpha;
}

dym *mk_dym_scalar_mult(dym *d,double alpha)
{
  dym *result;
  check_dym_code(d,"mk_dym_scalar_mult");
  result = mk_dym(d->rows,d->cols);
  dym_scalar_mult(d,alpha,result);
  return(result);
}

void dym_scalar_add(dym *d, double alpha, dym *r_d)
{
  int i,j;
  assert_dym_shape(r_d,d->rows,d->cols,"dym_scalar_add");
  for ( i = 0 ; i < r_d -> rows ; i++ )
    for ( j = 0 ; j < r_d -> cols ; j++ )
      r_d -> tdarr[i][j] = d->tdarr[i][j] + alpha;
}

dym *mk_dym_scalar_add(dym *d,double alpha)
{
  dym *result;
  check_dym_code(d,"mk_dym_scalar_add");
  result = mk_dym(d->rows,d->cols);
  dym_scalar_add(d,alpha,result);
  return(result);
}

void copy_dym(dym *d, dym *r_d)
{
  assert_dym_shape(r_d,d->rows,d->cols,"copy_dym");
  dym_scalar_mult(d,1.0,r_d);
}
    
dym *mk_copy_dym(dym *d)
{
  check_dym_code(d,"mk_copy_dym");
  return(mk_dym_scalar_mult(d,1.0));
}
    
void dym_plus(dym *d_1, dym *d_2, dym *r_d)
{
  int i,j;
  if ( d_1 -> rows != d_2 -> rows ||
       d_1 -> cols != d_2 -> cols 
     )
  {
    fprintf_dym(stderr,"d_1",d_1,"\n");
    fprintf_dym(stderr,"d_2",d_2,"\n");
    my_error("dym_plus: dyms (DYnamic Matrices) different shape");
  }

  assert_dym_shape(r_d,d_1->rows,d_1->cols,"dym_plus");
  for ( i = 0 ; i < r_d -> rows ; i++ )
    for ( j = 0 ; j < r_d -> cols ; j++ )
      r_d -> tdarr[i][j] = d_1->tdarr[i][j] + d_2 -> tdarr[i][j];
}

dym *mk_dym_plus(dym *a,dym *b)
{
  dym *result = mk_dym(a->rows,a->cols);
  dym_plus(a,b,result);
  return(result);
}

void dym_subtract(dym *d_1,dym *d_2,dym *r_d)
{
  dym *a = mk_dym_scalar_mult(d_2,-1.0);
  dym_plus(d_1,a,r_d);
  free_dym(a);
}

dym *mk_dym_subtract(dym *a,dym *b)
{
  dym *result = mk_dym(a->rows,a->cols);
  dym_subtract(a,b,result);
  return(result);
}

void dym_times_dyv(dym *a,dyv *b,dyv *result)
{
  int i;
  dyv *temp = mk_dyv(a->rows); 
             /* We need a copy in case b and result share memory */

  if ( a->cols != b -> size )
    my_error("dym_times_dyv: sizes wrong");
  assert_dyv_shape(result,a->rows,"dym_times_dyv");

  for ( i = 0 ; i < a->rows ; i++ )
  {
    double sum = 0.0;
    int j;
    for ( j = 0 ; j < a->cols ; j++ )
      sum += a->tdarr[i][j] * b->farr[j];
    temp->farr[i] = sum;
  }

  copy_dyv(temp,result);
  free_dyv(temp);
}

dyv *mk_dym_times_dyv(dym *a,dyv *b)
{
  dyv *result = mk_dyv(a->rows);
  dym_times_dyv(a,b,result);
  return(result);
}

void dyv_outer_product(dyv *a, dyv *b, dym *r_d)
{
  int i,j;
  assert_dym_shape(r_d,a->size,b->size,"dyv_outer_product");
  for (i=0;i<a->size;i++)
    for (j=0;j<b->size;j++)
      dym_set(r_d,i,j,dyv_ref(a,i)*dyv_ref(b,j));
}

dym *mk_dyv_outer_product(dyv *a, dyv *b)
{
  dym *result = mk_dym(a->size,b->size);
  dyv_outer_product(a,b,result);
  return result;
}

void dym_mult(dym *d_1, dym *d_2, dym *r_d)
{
  dym *a = mk_dym_mult(d_1,d_2);
             /* Note we have to first do the multiplying to the result
                a, in case the routine was called with d_1's memory
                = r_d's memory or d_2's memory = r_d's memory */
  assert_dym_shape(r_d,d_1 -> rows,d_2 -> cols,"dym_mult");

  copy_dym(a,r_d);
  free_dym(a);
}

dym *mk_dym_mult(dym *a,dym *b)
{
  int nrows = dym_rows(a), ncols = dym_cols(b);
  int acols = dym_cols(a);
  dym *c = mk_dym(nrows,ncols);
  int i,j;

  if ( acols != b->rows )
  {
    fprintf_dym(stderr,"a",a,"\n");
    fprintf_dym(stderr,"b",b,"\n");
    my_error("dym_mult: dyms (DYnamic Matrices) wrong shape\n");
  }

  for ( i = 0 ; i < nrows ; i++ )
    for ( j = 0 ; j < ncols ; j++ )
    {
      double sum = 0.0;
      int k;

      for ( k = 0 ; k < acols ; k++ )
        sum += a->tdarr[i][k] * b->tdarr[k][j];

      c->tdarr[i][j] = sum;
    }

  return c;
}

void dym_transpose(dym *d, dym *r_d)
{
  dym *a = mk_dym(d->cols,d->rows);
             /* Note we have to first do the transpose to the result
                a, in case the routine was called with d's memory
                = r_d's memory */
  int i,j;

  assert_dym_shape(r_d,d->cols,d->rows,"dym_transpose");

  for ( i = 0 ; i < d -> rows ; i++ )
    for ( j = 0 ; j < d -> cols ; j++ )
      a->tdarr[j][i] = d->tdarr[i][j];

  copy_dym(a,r_d);
  free_dym(a);
}

dym *mk_dym_transpose(dym *a)
{
  dym *result = mk_dym(a->cols,a->rows);
  dym_transpose(a,result);
  return(result);
}
void sing_val_w_check(dyv *w_diag,double thresh,dyv *winv_diag)
/*
   Sets winv_diag[i] = 1 / w_diag[i] for each i, EXCEPT, let w_max =
     magnitude of w_diag[i] with largest magnitude.

   Forall w_diag[i] such that w_diag[i] <= thresh * w_max, winv_diag[i] = 0.0
*/
{
  double max_mag = 0.0;
  double zero_below_me;  
  bool be_verbose = Verbosity > 100.0;
  int i;
  
  assert_dyv_shape(winv_diag,w_diag->size,"sing_val_w_check");

  for ( i = 0 ; i < w_diag->size ; i++ )
    max_mag = real_max(max_mag,fabs(w_diag->farr[i]));

  zero_below_me = thresh * max_mag;

  if ( be_verbose )
    printf("dym.c: SVD: thresh=%g, max_mag=%g\n",thresh,max_mag);

  for ( i = 0 ; i < w_diag->size ; i++ )
  {
    if ( be_verbose )
      printf("w_diag->farr[%d] = %g",i,w_diag->farr[i]);

    if ( w_diag->farr[i] == 0.0 || 
         fabs(w_diag->farr[i]) <= zero_below_me
       )
    {
      winv_diag->farr[i] = 0.0;
      if ( be_verbose ) printf("Singular value... zeroing!");
    }
    else
      winv_diag->farr[i] = 1.0 / w_diag->farr[i];

    if ( be_verbose ) printf("\n");
  }
}
  
void dym_scale_rows(dym *d,dyv *w_diag,dym *r_d)
/*
    Diag(w_diag) * d is copied to r_d
*/
{
  int i,j;
  assert_dym_shape(r_d,d->rows,d->cols,"dym_scale_rows::r_d");
  assert_dyv_shape(w_diag,d->rows,"dym_scale_rows::w_diag");
  for ( i = 0 ; i < d->rows ; i++ )
    for ( j = 0 ; j < d->cols ; j++ )
      r_d->tdarr[i][j] = d->tdarr[i][j] * w_diag->farr[i];
}

dym *mk_dym_scale_rows(dym *d,dyv *w_diag)
/*
    Returns Diag(w_diag) * d
*/
{
  dym *result = mk_dym(d->rows,d->cols);
  assert_dyv_shape(w_diag,d->rows,"mk_dym_scale_rows::w_diag");
  dym_scale_rows(d,w_diag,result);
  return(result);
}

void dym_scale_cols(dym *d,dyv *w_diag,dym *r_d)
/*
    d * Diag(w_diag) is copied to r_d
*/
{
  int i,j;
  assert_dym_shape(r_d,d->rows,d->cols,"dym_scale_cols::r_d");
  assert_dyv_shape(w_diag,d->cols,"dym_scale_cols::w_diag");
  for ( i = 0 ; i < d->rows ; i++ )
    for ( j = 0 ; j < d->cols ; j++ )
      r_d->tdarr[i][j] = d->tdarr[i][j] * w_diag->farr[j];
}

dym *mk_dym_scale_cols(dym *d,dyv *w_diag)
/*
    Returns d * Diag(w_diag)
*/
{
  dym *result = mk_dym(d->rows,d->cols);
  assert_dyv_shape(w_diag,d->cols,"mk_dym_scale_cols::w_diag");
  dym_scale_cols(d,w_diag,result);
  return(result);
}

void copy_dym_to_farr(dym *d,double *farr)
/* Copies either a column dym or a 1-row dym to an array of doubles.
   It is an error to call with neither rows nor columns equal to 1
*/
{
  if ( d->rows != 1 && d->cols != 1 )
  {
    fprintf_dym(stderr,"d",d,"\n");
    my_error("dym_to_realnums(): should be 1-columns or 1-row");
  }
  else if ( d->rows == 1 )
    copy_realnums(d->tdarr[0],farr,d->cols);
  else
  {
    int i;
    for ( i = 0 ; i < d->rows ; i++ )
      farr[i] = d->tdarr[i][0];
  }
}

void dym_solve_vector(dym *a, dyv *b, dyv *r_x)
/*
   Sets r_x so that (in the singular value decomp sense)
   it is as close as possible to a solution of
    a * r_x = b
*/
{
  dym *u,*v;
  dyv *w_diag;

  if ( a -> rows != b -> size )
  {
    fprintf(stderr,"amdm.c :: dym_solve, solving a x = b\n");
    fprintf(stderr,"a and b have incompatible numbers of rows\n");
    fprintf_dym(stderr,"a",a,"\n");
    fprintf_dyv(stderr,"b",b,"\n");
    my_error("dym_solve(a,b,r_x)");
  }

  assert_dyv_shape(r_x,a->cols,"dym_solve_vector");

  make_svd_components(a,&u,&w_diag,&v);
  dym_svd_backsub(u,w_diag,v,b,r_x);
  free_dym(u);
  free_dym(v);
  free_dyv(w_diag);
}

dyv *mk_dym_solve_vector(dym *a, dyv *b)
{
  dyv *result = mk_dyv(a->cols);
  dym_solve_vector(a,b,result);
  return(result);
}

void invert_dym(dym *d,dym *r_d)
/* 
   PRECONDITION: d must be a square dym.
     if d is singular or near singular returns the least squares inverse.
*/
{
  check_dym_code(d,"invert_dym");
  assert_dym_shape(r_d,d->rows,d->cols,"invert_dym");

  if ( d->rows != d->cols )
  {
    fprintf_dym(stderr,"d",d,"\n");
    my_error("invert_dym(): the bove dym is not square");
  }

  if ( d->rows == 1 )
  {
    double x = dym_ref(d,0,0);
    dym_set(r_d,0,0,(fabs(x)<DYM_SVD_NULL_THRESH) ? 0.0 : 1.0 / x);
  }
  else
  {
    dym *u,*v;
    dyv *w_diag;
    dyv *winv_diag;

    make_svd_components(d,&u,&w_diag,&v);

    if ( Verbosity > 2000.0 )
    {
      fprintf_dym(stdout,"u",u,"\n");
      wait_for_key();
      fprintf_dyv(stdout,"w_diag",w_diag,"\n");
      wait_for_key();
      fprintf_dym(stdout,"v",v,"\n");
      wait_for_key();
    }

    winv_diag = mk_dyv(d->cols);
    sing_val_w_check(w_diag,DYM_SVD_NULL_THRESH,winv_diag);

    if ( Verbosity > 2000.0 )
    {
      fprintf_dyv(stdout,"winv_diag",winv_diag,"\n");
      wait_for_key();
    }

    dym_transpose(u,r_d);
    if ( Verbosity > 2000.0 )
    {
      fprintf_dym(stdout,"u^T",r_d,"\n");
      wait_for_key();
    }

    dym_scale_rows(r_d,winv_diag,r_d);

    if ( Verbosity > 2000.0 )
    {
      fprintf_dym(stdout,"Diag(winv_diag) * u^T",r_d,"\n");
      wait_for_key();
    }

    dym_mult(v,r_d,r_d);
    if ( Verbosity > 2000.0 )
    {
      fprintf_dym(stdout,"v * winv * u^T",r_d,"\n");
      wait_for_key();
    }

    free_dym(u);
    free_dyv(w_diag);
    free_dym(v);
    free_dyv(winv_diag);
  }
}

dym *mk_invert_dym(dym *d)
{
  dym *result = mk_dym(d->rows,d->cols);
  invert_dym(d,result);
  return(result);
}

bool is_same_realnum(double x,double y)
{
  double eps = 1e-3;

  bool result = TRUE;

  if ( fabs(x) < eps && fabs(y) < eps )
    result = TRUE;
  else if ( fabs(x) < eps && fabs(y) > 2.0 * eps )
    result = FALSE;
  else if ( fabs(y) < eps && fabs(x) > 2.0 * eps  )
    result = FALSE;
  else if ( x < 0.0 && y > 0.0 )
    result = FALSE;
  else if ( y < 0.0 && x > 0.0 )
    result = FALSE;
  else if ( fabs(x) > fabs(y) )
    result = fabs(x) < (1.0 + eps) * fabs(y);
  else
    result = fabs(y) < (1.0 + eps) * fabs(x);

  return(result);
}

bool is_dym_symmetric(dym *d)
{
  int i,j;
  bool result = TRUE;

  if ( d->rows != d->cols )
  {
    fprintf_dym(stderr,"d = ",d,"\n");
    my_error("dym.c is_dym_symmetric(d). Requires a square dym.\n");
  }

  for ( i = 0 ; result && i < d->rows ; i++ )
    for ( j = 0 ; result && j < i ; j++ )
    {
      result = is_same_realnum(d->tdarr[i][j],d->tdarr[j][i]);
      if ( Verbosity > 70.0 )
        printf("i = %d, j = %d, result = %d\n",i,j,result);
    }

  if ( !result ) printf("Not Symmetric\n");

  return(result);
}

/* Solve a lower triangular system of equations Lx=b.
   It's OK if x and b are the same dyv: b will be overwritten
   with the solution.
   The part of L above the main diagonal is ignored. */
void tri_backsub(dym *l, dyv *b, dyv *x)
{
  int i, j;
  double tmp;

  if (dym_rows(l) != dym_cols(l))
    my_error("matrix must be square in tri_backsub");
  if (dyv_size(x) != dym_rows(l))
    my_error("vector must be same size as matrix in tri_backsub");

  for (i = 0; i < dym_rows(l); i++) {
    tmp = 0;
    for (j = 0; j < i; j++)
      tmp += dyv_ref(x, j) * dym_ref(l, i, j);
    dyv_set(x, i, (dyv_ref(b, i) - tmp) / dym_ref(l, i, i));
  }
}

dyv *mk_tri_backsub(dym *l, dyv *b)
{
  dyv *res = mk_dyv(dyv_size(b));
  tri_backsub(l, b, res);
  return res;
}

/* Solve an upper triangular system of equations L'x=b.
   It's OK if x and b are the same dyv: b will be overwritten
   with the solution.
   The part of L above the main diagonal is ignored. */
void tri_trans_backsub(dym *l, dyv *b, dyv *x)
{
  int i, j;
  double tmp;

  if (dym_rows(l) != dym_cols(l))
    my_error("matrix must be square in tri_trans_backsub");
  if (dyv_size(x) != dym_rows(l))
    my_error("vector must be same size as matrix in tri_trans_backsub");

  for (i = dym_cols(l)-1; i >= 0; i--) {
    tmp = 0;
    for (j = i+1; j < dym_rows(l); j++)
      tmp += dyv_ref(x, j) * dym_ref(l, j, i);
    dyv_set(x, i, (dyv_ref(b, i) - tmp) / dym_ref(l, i, i));
  }
}

dyv *mk_tri_trans_backsub(dym *l, dyv *b)
{
  dyv *res = mk_dyv(dyv_size(b));
  tri_trans_backsub(l, b, res);
  return res;
}

/* Solve the symmetric positive definite system of equations
   LL'x = b, where L is lower triangular.
   It's OK if x and b are the same dyv: b will be overwritten
   with the solution.
   The part of L above the main diagonal is ignored. */
void cholesky_backsub(dym *l, dyv *b, dyv *x)
{
  tri_backsub(l, b, x);
  tri_trans_backsub(l, x, x);
}

dyv *mk_cholesky_backsub(dym *l, dyv *b)
{
  dyv *res = mk_dyv(dyv_size(b));
  cholesky_backsub(l, b, res);
  return res;
}

/* Solve the symmetric positive definite system of equations Ax = b.
   If A is not symmetric or positive definite generate an error.
   It's OK if x and b are the same dyv: b will be overwritten
   with the solution. */
void cholesky_solve(dym *a, dyv *b, dyv *x)
{
  dym *l = mk_dym(dym_rows(a), dym_rows(a));
  if (!attempt_cholesky_decomp(a, l))
    my_error("cholesky_solve requires a symmetric +ve definite matrix");
  cholesky_backsub(l, b, x);
  free_dym(l);
}

dyv *mk_cholesky_solve(dym *a, dyv *b)
{
  dyv *res = mk_dyv(dyv_size(b));
  cholesky_solve(a, b, res);
  return res;
}

/*
   PRE: lot is a lower triangular matrix. (i.e. i < j => lot[i,j] = 0)

      [ lot00     0     0  ...  ]
      [ lot10 lot11     0   ... ]
      [ lot20 lot21 lot22  ...  ]
      [  :       :         ...  ]

   This function returns NULL if lot is singular.
   (note: lot is singular iff for some i, lotii == 0)

   Else, this function returns (lot lot^T)^-1 b,
   i.e. it returns x such that

      lot^t lot x = b
*/
dyv *mk_solve_lot_equation(dym *lot,dyv *b)
{
  int size = dym_rows(lot);
  int i;
  dyv *x = NULL;
  bool singular = FALSE;

  if ( size != dyv_size(b) || size != dym_cols(lot) )
    my_error("mk_solve_lot_eqation: precondition violation");

  for ( i = 0 ; i < size && !singular ; i++ )
  {
    if ( dym_ref(lot,i,i) == 0.0 )
      singular = TRUE;
  }

  if ( !singular )
  {
    dyv *y = mk_dyv(size);
    x = mk_dyv(size);

    /* First solve L y = b for y */
    for ( i = 0 ; i < size ; i++ )
    {
      /* Compute y[i]... */
      double sum_wij_yj = 0.0;
      int j;
      for ( j = 0 ; j < i ; j++ )
        sum_wij_yj += dym_ref(lot,i,j) * dyv_ref(y,j);
      dyv_set(y,i, (dyv_ref(b,i) - sum_wij_yj) / dym_ref(lot,i,i));
    }

    /* Now solve L^T x = y for x. Then we'll have L L^T x = b
       as desired */

    for ( i = size-1 ; i >= 0 ; i -= 1 )
    {
      double sum_wij_yj = 0.0;
      int j;
      for ( j = i+1 ; j < size ; j++ )
        sum_wij_yj += dym_ref(lot,j,i) * dyv_ref(x,j);
      dyv_set(x,i, (dyv_ref(y,i) - sum_wij_yj) / dym_ref(lot,i,i));
    }

    free_dyv(y);
  }

  return(x);
}

bool attempt_cholesky_decomp(dym *a,dym *r_l)
/*
   Returns FALSE if a is not symmetric or not positive definite.
   If it is symmetric and positive definite, returns TRUE and also
   returns a lower triangular dym in r_l such that

     r_l r_l^T = a

     [r_l]_ij = 0 for all i,j such that j > i

   As usual, it is fine if a and r_l share the same memory.
*/
{
  bool result = is_dym_symmetric(a);
  dym *l = mk_dym(a->rows,a->cols);
  double maxmat = 1e-6;
  int i,j;

  for ( i = 0 ; result && i < a->rows ; i++ )
    for ( j = i ; result && j < a->rows ; j++ )
      maxmat = real_max(fabs(a->tdarr[i][j]),maxmat);

  zero_dym(l);

  if ( Verbosity > 50.0 )
    printf("is_dym_symmetric = %d\n",result);

  if ( result )
  {
       /* COPIED from numerical recipes in C */
    int k;

    for ( i = 0 ; result && i < a->rows ; i++ )
      for ( j = i ; result && j < a->rows ; j++ )
      {
        double sum = a->tdarr[i][j];
        for ( k = i - 1 ; k >= 0 ; k-- )
          sum -= l->tdarr[i][k] * l->tdarr[j][k];

        if ( Verbosity > 50.0 )
          printf("i = %d , j = %d , sum = %g\n",i,j,sum);

        if ( i == j )
        {
          if ( sum < -1e-6 * maxmat )
            result = FALSE;
          else if ( sum < 0.0 )
            l->tdarr[i][j] = 0.0;
          else
            l->tdarr[i][j] = sqrt(sum);
        }
        else if ( l->tdarr[i][i] == 0.0 )
          result = FALSE;
        else
          l->tdarr[j][i] = sum / l->tdarr[i][i];
      }
  }

  if ( result )
    copy_dym(l,r_l);
  else
    zero_dym(r_l);

  free_dym(l);
  return(result);
}

bool is_dym_symmetric_positive_definite(dym *d)
{
  dym *l = mk_dym(d->rows,d->cols);
  bool result = attempt_cholesky_decomp(d,l);
  free_dym(l);
  return(result);
}

/* Use cholesky decomposition to invert a symmetric,
   positive, definite dym.  Returns TRUE on success.
   Failure likely means the original matrix was not spd.
 */
bool invert_spd_cholesky(dym *d, dym *r_d)
{
  int i,j,k;
  dym *l = mk_dym(dym_cols(r_d),dym_cols(r_d));
  dym *l_inv = mk_zero_dym(dym_cols(r_d),dym_cols(r_d));
  dym *l_inv_t = NULL;
  bool result = attempt_cholesky_decomp(d,l);

  for ( i = 0 ; i < dym_rows(l) ; i++ )
    result = result && dym_ref(l,i,i) != 0.0;

  if (result)
  {
    /* this little chunk taken from sec 2.9 of Num Rec in C */
    for (i=0;i<dym_rows(d);i++)
    {
      dym_set(l_inv,i,i,1.0/dym_ref(l,i,i));
      for (j=i+1;j<dym_rows(d);j++)
      {
        double sum = 0.0;
        for (k=i;k<j;k++) sum -= dym_ref(l,j,k) * dym_ref(l_inv,k,i);
        dym_set(l_inv,j,i,sum/(dym_ref(l,j,j)));
      }
    }
    l_inv_t = mk_dym_transpose(l_inv);
    dym_mult(l_inv_t,l_inv,r_d);
  }
  else zero_dym(r_d);

  free_dym(l); free_dym(l_inv); 
  if (l_inv_t != NULL)
    free_dym(l_inv_t);

  return result;
}

/* returning NULL indicates failure */
dym *mk_invert_spd_cholesky(dym *d)
{
  dym *result = mk_dym(dym_rows(d),dym_cols(d));
  bool ok = invert_spd_cholesky(d,result);
  if (ok) return result;
  else
  {
    free_dym(result);
    return NULL;
  }
}

void enforce_dym_symmetry(dym *d)
/*
   Let d' = matrix after execution
   Let d = matrix before execution

   d' = 0.5 * ( d + d^T )

   so d'[i][j] = 0.5 * (d[i][j] + d[j][i])
*/
{
  int i,j;
  if ( dym_rows(d) != dym_cols(d) )
  {
    fprintf_dym(stderr,"d",d,"\n");
    my_error("enforce_dym_symmetry(): the bove dym is not square");
  }

  for ( i = 0 ; i < d->rows ; i++ )
    for ( j = 0 ; j < i ; j++ )
      dym_set(d,i,j,(dym_ref(d,i,j) + dym_ref(d,j,i))/2.0);

  for ( i = 0 ; i < d->rows ; i++ )
    for ( j = i+1 ; j < d->rows ; j++ )
      dym_set(d,i,j,dym_ref(d,j,i));
}

void invert_symmetric_dym(dym *d,dym *r_d)
/* 
   PRECONDITION: d must be a square symmetric dym.
   if d is singular or near singular returns the least squares inverse.
*/
{
  if ( !is_dym_symmetric(d) )
  {
    fprintf_dym(stderr,"d",d,"\n");
    my_error("invert_symmetric_dym: was given a non-symmetric dym");
  }
  else
  {
    invert_dym(d,r_d);
    enforce_dym_symmetry(r_d);
  }
}

dym *mk_invert_symmetric_dym(dym *d)
{
  dym *result = mk_dym(d->rows,d->cols);
  invert_symmetric_dym(d,result);
  return(result);
}

/******* printing dyms (DYnamic Matrices) *********/

typedef struct buftab_struct
{
  int rows,cols;
  char ***strings;
} buftab;

static char *bufstr(buftab *bt,int i,int j)
{
  char *result;

  if ( i < 0 || i >= bt->rows || j < 0 || j >= bt->cols )
  {
    result = NULL;
    my_error("bufstr()");
  }
  else
    result = bt->strings[i][j];

  if ( result == NULL )
    result = "-";

  return(result);
}

static void fprint_buftab(
    FILE *s,
    buftab *bt
  )
{
  int *widths = AM_MALLOC_ARRAY(int,bt->cols);
  int i,j;
  set_ints_constant(widths,bt->cols,0);

  for ( i = 0 ; i < bt->rows ; i++ )
    for ( j = 0 ; j < bt->cols ; j++ )
      widths[j] = int_max(widths[j],strlen(bufstr(bt,i,j)));

  for ( i = 0 ; i < bt->rows ; i++ )
    for ( j = 0 ; j < bt->cols ; j++ )
    {
      char ford[20];
      sprintf(ford,"%%%ds%s",widths[j],(j==bt->cols-1) ? "\n" : " ");
      fprintf(s,ford,bufstr(bt,i,j));
    }

  am_free((char *)widths,sizeof(int) * bt->cols);
}

static void init_buftab(
    buftab *bt,
    int rows,
    int cols
  )
{
  if ( rows < 0 || cols < 0 )
    my_error("init_buftab()");
  else
  {
    int i,j;
    bt -> rows = rows;
    bt -> cols = cols;
    bt -> strings = AM_MALLOC_ARRAY(char_ptr_ptr,rows);
    for ( i = 0 ; i < rows ; i++ )
      bt->strings[i] = AM_MALLOC_ARRAY(char_ptr,cols);
    for ( i = 0 ; i < rows ; i++ )
      for ( j = 0 ; j < cols ; j++ )
        bt->strings[i][j] = NULL;
  }
}

static void free_buftab_contents(buftab *bt)
{
  int i,j;
  for ( i = 0 ; i < bt->rows ; i++ )
    for ( j = 0 ; j < bt->cols ; j++ )
      if ( bt->strings[i][j] != NULL )
        am_free((char *)(bt->strings[i][j]),
                sizeof(char) * (strlen(bt->strings[i][j]) + 1)
               );

  for ( i = 0 ; i < bt->rows ; i++ )
    am_free((char *)(bt->strings[i]),sizeof(char_ptr) * bt->cols);
    
  am_free((char *)bt->strings,sizeof(char_ptr_ptr) * bt->rows);
}

static void set_buftab(
    buftab *bt,
    int i,
    int j,
    const char *str
  )
{
  if ( i < 0 || i >= bt->rows || j < 0 || j >= bt->cols )
    my_error("set_buftab()");
  else if ( bt->strings[i][j] != NULL )
    my_error("set_buftab: non null string");
  else
    bt->strings[i][j] = make_copy_string(str);
}

void fprintf_dym(FILE *s, const char *m1, const dym *d, const char *m2)
{
  if ( d == NULL )
    fprintf(s,"%s = (dym *)NULL%s",m1,m2);
  else if ( d->dym_code != DYM_CODE )
  {
    fprintf(stderr,"fprintf_dym(s,\"%s\",d,\"\\n\"\n",m1);
    my_error("fprintf_dym called with a non-allocated dym (DYnamic Matrix)");
  }
  else if ( d->rows <= 0 || d->cols <= 0 )
    fprintf(s,"%s = <Dym with %d row%s and %d column%s>%s",
            m1,d->rows,(d->rows==-1)?"":"s",
            d->cols,(d->cols==-1)?"":"s",m2
           );
  else
  {
    int i;
    buftab bt;

    init_buftab(&bt,d->rows,d->cols + 4);

    for ( i = 0 ; i < d->rows ; i++ )
    {
      int j;
      set_buftab(&bt,i,0,(i == (d->rows-1)/2) ? m1 : "");
      set_buftab(&bt,i,1,(i == (d->rows-1)/2) ? "=" : "");
      set_buftab(&bt,i,2,"[");

      for ( j = 0 ; j < d -> cols ; j++ )
      {
        char buff[100];
        sprintf(buff," %g ",d->tdarr[i][j]);
        set_buftab(&bt,i,3+j,buff);
      }

      set_buftab(&bt,i,3+d->cols,"]");
    }

    fprint_buftab(s,&bt);
    free_buftab_contents(&bt);
  }
  fprintf(s,"\n");
}


void pdym(dym *d)
{
  fprintf_dym(stdout,"dym",d,"\n");
}


void fprintf_dym_and_confidence(
    FILE *s,
    char *m1,
    dym *d,
    dym *conf,
    bool huge_uncertainty,
    char *m2
  )
{
  if ( d == NULL )
    fprintf(s,"%s = (dym *)NULL%s",m1,m2);
  else if ( d->rows <= 0 || d->cols <= 0 )
    fprintf(s,"%s = <Dym with %d row%s and %d column%s>%s",
            m1,d->rows,(d->rows==1)?"":"s",
            d->cols,(d->cols==1)?"":"s",m2
           );
  else if ( d->rows != conf->rows || d->cols != conf->cols )
    my_error("fprintf_dym_and_confidence(). d and conf differ in shape");
  else
  {
    int i;
    buftab bt;

    init_buftab(&bt,d->rows,3 * d->cols + 4);

    for ( i = 0 ; i < d->rows ; i++ )
    {
      int j;
      set_buftab(&bt,i,0,(i == (d->rows-1)/2) ? m1 : "");
      set_buftab(&bt,i,1,(i == (d->rows-1)/2) ? "=" : "");
      set_buftab(&bt,i,2,"[");

      for ( j = 0 ; j < d -> cols ; j++ )
      {
        char buff[100];
        sprintf(buff," %g",d->tdarr[i][j]);
        set_buftab(&bt,i,3 + 3 * j,buff);
        set_buftab(&bt,i,4 + 3 * j,"+/-");
        if ( huge_uncertainty )
          sprintf(buff,"huge uncertainty ");
        else
          sprintf(buff,"%g ",conf->tdarr[i][j]);
        set_buftab(&bt,i,5 + 3 * j,buff);
      }

      set_buftab(&bt,i,3 + 3 * d->cols,"]");
    }

    fprint_buftab(s,&bt);
    free_buftab_contents(&bt);
  }
  fprintf(s,"\n");
}

void fprintf_dym_dym(FILE *s,char *m1,dym *d1,char *m2,dym *d2,char *m3)
{
  int maxrows = int_max(dym_rows(d1),dym_rows(d2));
  int crow = int_min(int_min(dym_rows(d1)-1,dym_rows(d2)-1),(maxrows-1)/2);
  buftab bt[1];
  int i;

  init_buftab(bt,maxrows,6 + dym_cols(d1) + dym_cols(d2));

  for ( i = 0 ; i < maxrows ; i++ )
  {
    int j;

    if ( i < dym_rows(d1) )
    {
      set_buftab(bt,i,0,(i == crow) ? m1 : "");
      set_buftab(bt,i,1,"[");

      for ( j = 0 ; j < dym_cols(d1) ; j++ )
      {
        char buff[100];
        sprintf(buff," %g",dym_ref(d1,i,j));
        set_buftab(bt,i,2 + j,buff);
      }

      set_buftab(bt,i,2 + dym_cols(d1),"]");
    }

    if ( i < dym_rows(d2) )
    {
      set_buftab(bt,i,3 + dym_cols(d1),(i == crow) ? m2 : "");
      set_buftab(bt,i,4 + dym_cols(d1),"[");

      for ( j = 0 ; j < dym_cols(d2) ; j++ )
      {
        char buff[100];
        sprintf(buff," %g",dym_ref(d2,i,j));
        set_buftab(bt,i,5 + dym_cols(d1) + j,buff);
      }

      set_buftab(bt,i,5 + dym_cols(d1) + dym_cols(d2),"]");
    }
  }

  fprint_buftab(s,bt);
  free_buftab_contents(bt);
  fprintf(s,"\n");
}

void fprintf_dym_dyv(FILE *s,char *m1,dym *d1,char *m2,dyv *d2,char *m3)
{
  dym *dm = mk_col_dym_from_dyv(d2);
  fprintf_dym_dym(s,m1,d1,m2,dm,m3);
  free_dym(dm);
}

/* defined here because needs buftab stuff */
void fprintf_dyv(FILE *s,char *m1,const dyv *d,char *m2)
{
  if ( d == NULL )
    fprintf(s,"%s = (dyv *)NULL%s",m1,m2);
  else if ( d->dyv_code != DYV_CODE )
  {
    fprintf(stderr,"fprintf_dyv(s,\"%s\",d,\"\\n\"\n",m1);
    my_error("fprintf_dyv called with a non-allocated dyv (DYnamic Vector)");
  }
  else if ( d->size <= 0 )
    fprintf(s,"%s = <Dyv of size %d>%s",m1,d->size,m2);
  else
  {
    int i;
    buftab bt;
    int cols = 1;

    init_buftab(&bt,d->size,cols + 4);

    for ( i = 0 ; i < d->size ; i++ )
    {
      char buff[100];
      set_buftab(&bt,i,0,(i == (d->size-1)/2) ? m1 : "");
      set_buftab(&bt,i,1,(i == (d->size-1)/2) ? "=" : "");
      set_buftab(&bt,i,2,"(");

      sprintf(buff," %g ",d->farr[i]);
      set_buftab(&bt,i,3,buff);

      set_buftab(&bt,i,3+cols,")");
    }

    fprint_buftab(s,&bt);
    free_buftab_contents(&bt);
  }
  fprintf(s,"\n");
}

void pdyv(dyv *d)
{
  fprintf_dyv(stdout,"dyv",d,"\n");
}

double dym_extreme(dym *dv, bool max)
{
  int i,j;
  double res = -77.0;
  if ( dym_rows(dv)*dym_cols(dv) < 1 ) my_error("dym_extreme: empty dym");
  for (i=0;i<dv->rows;i++)
    for (j=0;j<dv->cols;j++)
      if ((i == 0 && j == 0) || (max && dym_ref(dv,i,j) > res)
	  || (!max && dym_ref(dv,i,j) < res))
	res = dym_ref(dv,i,j);
  return res;
}

double dym_min(dym *dv)
{
  return(dym_extreme(dv, FALSE));
}

double dym_max(dym *dv)
{
  return(dym_extreme(dv, TRUE));
}

dym *mk_dym_11(double x00)
{
  dym *result = mk_dym(1,1);
  dym_set(result,0,0,x00);
  return(result);
}

dym *mk_dym_12(double x00, double x01)
{
  dym *result = mk_dym(1,2);
  dym_set(result,0,0,x00);
  dym_set(result,0,1,x01);
  return(result);
}

dym *mk_dym_13(double x00, double x01, double x02)
{
  dym *result = mk_dym(1,3);
  dym_set(result,0,0,x00);
  dym_set(result,0,1,x01);
  dym_set(result,0,2,x02);
  return(result);
}

dym *mk_dym_21(double x00,
               double x10)
{
  dym *result = mk_dym(2,1);
  dym_set(result,0,0,x00);
  dym_set(result,1,0,x10);
  return(result);
}

dym *mk_dym_22(double x00, double x01,
               double x10, double x11)
{
  dym *result = mk_dym(2,2);
  dym_set(result,0,0,x00);
  dym_set(result,0,1,x01);
  dym_set(result,1,0,x10);
  dym_set(result,1,1,x11);
  return(result);
}

dym *mk_dym_23(double x00, double x01, double x02,
               double x10, double x11, double x12)
{
  dym *result = mk_dym(2,3);
  dym_set(result,0,0,x00);
  dym_set(result,0,1,x01);
  dym_set(result,0,2,x02);
  dym_set(result,1,0,x10);
  dym_set(result,1,1,x11);
  dym_set(result,1,2,x12);
  return(result);
}

dym *mk_dym_31(double x00,
               double x10,
               double x20)
{
  dym *result = mk_dym(3,1);
  dym_set(result,0,0,x00);
  dym_set(result,1,0,x10);
  dym_set(result,2,0,x20);
  return(result);
}

dym *mk_dym_32(double x00, double x01,
               double x10, double x11,
               double x20, double x21)
{
  dym *result = mk_dym(3,2);
  dym_set(result,0,0,x00);
  dym_set(result,0,1,x01);
  dym_set(result,1,0,x10);
  dym_set(result,1,1,x11);
  dym_set(result,2,0,x20);
  dym_set(result,2,1,x21);
  return(result);
}

dym *mk_dym_33(double x00, double x01, double x02,
               double x10, double x11, double x12,
               double x20, double x21, double x22)
{
  dym *result = mk_dym(3,3);
  dym_set(result,0,0,x00);
  dym_set(result,0,1,x01);
  dym_set(result,0,2,x02);
  dym_set(result,1,0,x10);
  dym_set(result,1,1,x11);
  dym_set(result,1,2,x12);
  dym_set(result,2,0,x20);
  dym_set(result,2,1,x21);
  dym_set(result,2,2,x22);
  return(result);
}

/* As dym_xt_a_x_value except that a is a diagonal matrix specified as a dyv */
double xt_diag_x_value(dyv *x, dyv *a)
{
  int i;
  double result = 0.0;
  for (i=0;i<x->size;i++) result += dyv_ref(x,i) * dyv_ref(a,i) * dyv_ref(x,i);
  return result;
}

double dym_xt_a_x_value(dyv *x,dym *a)
/* Returns x^T a x */
{
  dyv *ax = mk_dym_times_dyv(a,x);
  double result = dyv_scalar_product(x,ax);
  free_dyv(ax);
  return(result);
}

void dym_ptq(dym *p,dym *q,dym *r_dym)
/*
   Sets r_dym to contain p^T q.

   r_dym must have the following shape: dym_cols(p) rows, and dym_cols(q) cols.
   dym_rows(p) must equal dym_rows(q);

   Note. (RQ)[i][j] = sum-over-k-of R[i][k] * Q[k][j]

   Put P^T = R, then (P^T)[i][k] = P[k][i], so

     (P^T)Q[i][j] = sum-over-k-of P[k][i] * Q[k][j]
*/
{
  dym *a = mk_dym(dym_cols(p),dym_cols(q));
             /* Note we have to first do the multiplying to the result
                a, in case the routine was called with d_1's memory
                = r_d's memory or d_2's memory = r_d's memory */

             /* Efficiency: could make this faster by comparing pointers and
                only foing copy if necessary
                Ditto true for multiply, invert, transpose etc
             */
  int i,j;

  if ( dym_rows(p) != dym_rows(q) )
  {
    fprintf_dym(stderr,"p",p,"\n");
    fprintf_dym(stderr,"q",q,"\n");
    my_error("dym_ptq: dyms (DYnamic Matrices) wrong shape\n");
  }

  assert_dym_shape(r_dym,dym_cols(p),dym_cols(q),"dym_ptq");

  for ( i = 0 ; i < dym_rows(a) ; i++ )
    for ( j = 0 ; j < dym_cols(a) ; j++ )
    {
      double sum = 0.0;
      int k;

      for ( k = 0 ; k < dym_rows(p) ; k++ )
        sum += dym_ref(p,k,i) * dym_ref(q,k,j);

      dym_set(a,i,j,sum);
    }

  copy_dym(a,r_dym);
  free_dym(a);
}

dym *mk_dym_ptq(dym *p,dym *q)
/*
   Makes and returns p^T q
*/
{
  dym *result = mk_dym(dym_cols(p),dym_cols(q));
  dym_ptq(p,q,result);
  return(result);
}

/*
   Sets r_dym to contain p^T q p.

   r_dym must have the following shape: dym_cols(p) rows, and dym_cols(q) cols.
   dym_rows(p) must equal dym_rows(q);

   Note. (RQ)[i][j] = sum-over-k-of R[i][k] * Q[k][j]

   Put P^T = R, then (P^T)[i][k] = P[k][i], so

     (P^T)Q[i][j] = sum-over-k-of P[k][i] * Q[k][j]
*/
void dym_ptqp(dym *p,dym *q,dym *r_dym)
{
  dym *a = mk_dym_ptq(p,q);
  dym_mult(a,p,r_dym);
  free_dym(a);
  return;
}

/* Makes and returns p^T q p */
dym *mk_dym_ptqp(dym *p,dym *q)
{
  dym *result = mk_dym(dym_cols(p),dym_cols(p));
  dym_ptqp(p,q,result);
  return(result);
}

void dym_transpose_times_dyv(dym *p,dyv *v,dyv *r_dyv)
/* Returns p^T v in r_dyv. 
  PRE: dym_rows(p) == dyv_size(v)
       dyv_cols(p) == dyv_size(r_dyv);

  (Rv)[i] = sum-over-k-of R[i][k] v[k]
  if P^T = R then R[i][k] = P[k][i] so

  (P^T v)[i] = sum-over-k-of P[k][i] v[k]
*/
{
  dyv *temp = mk_dyv(dyv_size(r_dyv));
  int i;
  assert_dym_shape(p,dyv_size(v),dyv_size(r_dyv),"dym_transpose_times_dyv");
  for ( i = 0 ; i < dyv_size(r_dyv) ; i++ )
  {
    double sum = 0.0;
    int k;
    for ( k = 0 ; k < dyv_size(v) ; k++ )
      sum += dyv_ref(v,k) * dym_ref(p,k,i);
    dyv_set(temp,i,sum);
  }
  copy_dyv(temp,r_dyv);
  free_dyv(temp);
}

dyv *mk_dym_transpose_times_dyv(dym *p,dyv *v)
/* Returns p^T v */
{
  dyv *result = mk_dyv(dym_cols(p));
  dym_transpose_times_dyv(p,v,result);
  return(result);
}

dym *mk_identity_dym(int dims)
{
  dyv *dv = mk_constant_dyv(dims,1.0);
  dym *result = mk_diag_dym_from_dyv(dv);
  free_dyv(dv);
  return(result);
}

dym *mk_dym_from_string(char *s)
{
  int state = 0;
  int maxcol = 0;
  int icol = 0;
  int irow = 1;
  int index,i = 0;
  double val;
  dym *res;

  /* find the size of the matrix */
  while(s[i])
  {
    if ((s[i] == ';') || (s[i] == '\n'))
    {
      irow++;
      if (state == 1) {icol++; state = 0;}
      if (icol > maxcol) maxcol = icol;
      icol = 0;
    }
    else if (isspace(s[i]))
    {
      if (state == 1) {icol++; state = 0;}
    }
    else state = 1;
    i++;
  }
  if (state == 1) icol++;
  if (icol > maxcol) maxcol = icol;

  res = mk_dym(irow,maxcol);

  /* fill in the array */
  i = index = irow = icol = 0;
  while(s[i])
  {
    if ((s[i] == ';') || (s[i] == '\n'))
    {
      if (state == 1)
      {
	sscanf(s+index,"%lf",&val);
	dym_set(res,irow,icol,val);
      }
      index = i;
      state = 0;
      icol = 0;
      irow++;
    }
    else if (isspace(s[i]))
    {
      if (state == 1) 
      {
	sscanf(s+index,"%lf",&val);
	dym_set(res,irow,icol,val);
	index = i;
	icol++; 
	state = 0;
      }
    }
    else state = 1;
    i++;
  }
  if (state == 1) 
  {
    sscanf(s+index,"%lf",&val);
    dym_set(res,irow,icol,val);
  }
    
  return res;
}

dym *mk_dym_from_args(char *name,int argc,char *argv[],dym *deflt)
/* COPIES in deflt (if so required) */
{
  bool name_there = index_of_arg(name,argc,argv) >= 0;
  dym *result;
  if ( !name_there )
    result = mk_copy_dym(deflt);
  else
  {
    int size = dym_rows(deflt) * dym_cols(deflt);
    dyv *dv = mk_basic_dyv_from_args(name,argc,argv,size);
    int i,j;
    if ( dv == NULL )
    {
      fprintf(stderr,"COMMAND LINE USER ERROR (it's YOUR fault)\n");
      fprintf(stderr,"...when attempting to read a dym identified by\n");
      fprintf(stderr,"the name \"%s\". Perhaps a non-number, or the\n",name);
      fprintf(stderr,"command line finished before all args found?\n");
      fprintf_dym(stderr,"deflt_dym",deflt,"\n");
      my_error("mk_dym_from_args()");
    }

    result = mk_dym(dym_rows(deflt),dym_cols(deflt));
    for ( i = 0 ; i < dym_rows(deflt) ; i++ )
      for ( j = 0 ; j < dym_cols(deflt) ; j++ )
        dym_set(result,i,j,dyv_ref(dv,i * dym_cols(deflt) + j));

    free_dyv(dv);
  }

  return(result);
}

#define TINY 1.0e-20;
bool numrec_ludcmp(double **a, int n, int *indx, double *d)
/* Returns TRUE if and only if it succeeded.
   Returns FALSE if the matrix is singular
*/
{
        int i,j,k;
        int imax = 0;   /* Prevent "used uninitialized" warning */
        double big,dum,sum,temp;
        double *vv;
        bool result = TRUE;

        vv=am_malloc_realnums(n+1);
        *d=1.0;
        for (i=1;i<=n && result;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs(a[i][j])) > big) big=temp;
                if (big == 0.0) result = FALSE;
                vv[i]=1.0/big;
        }
        for (j=1;j<=n && result;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                }
                big=0.0;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                        if ( (dum=vv[i]*fabs(sum)) >= big) {
                                big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum=a[imax][k];
                                a[imax][k]=a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                indx[j]=imax;
                if (a[j][j] == 0.0) a[j][j]=TINY;
                if (j != n) {
                        dum=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] *= dum;
                }
        }
        am_free_realnums(vv,n+1);
  
        return(result);
}
#undef TINY

bool is_determinant_positive(dym *a)
{
  int size = dym_rows(a);
  double **nra = am_malloc_2d_realnums(size+1,size+1);
  int *indx = am_malloc_ints(size+1);
  double d;
  int i,j;
  bool nonsing;
  bool result;

  if ( size != dym_cols(a) ) my_error("dym_determinant: a should be square");

  for ( i = 0 ; i < size ; i++ )
    for ( j = 0 ; j < size ; j++ )
      nra[i+1][j+1] = dym_ref(a,i,j);

  nonsing = numrec_ludcmp(nra,size,indx,&d);

  if ( nonsing )
  {
    for ( j = 1 ; j <= size ; j++ )
      d *= nra[j][j];
    result = d >= 0.0;
  }
  else
    result = TRUE;

  am_free_2d_realnums(nra,size+1,size+1);
  am_free_ints(indx,size+1);
  return(result);
}

double dym_determinant(dym *dm)
/*
   Note this is a pretty expensive way to get a determinant. I use SVD to
   get me the magnitude of the determinant. I then use LU decomp to get me
   the sign. Why go this way? Why not just use LU decomp for the whole thing?
   Because it's
   crrrrap , that's why. det ( 3333 3333 5555 5555 ) = 1.3562 according
   to numrec's ludcmp.

   I should probably be able to deduce the sign from the SVD decomposition.
   If *YOU* are reading this code and *YOU* can figure out how to do it,
   please let me know. (-awm).
*/
{
  dym *u,*v;
  dyv *w_diag;
  double result = 1.0;
  int i;

  if ( dym_rows(dm) != dym_cols(dm) )
    my_error("dym_determinant: you gave me a non-square dym");

  make_svd_components(dm,&u,&w_diag,&v);

  if ( Verbosity > 60.0 )
  {
    fprintf_dym(stdout,"determinant dm",dm,"\n");
    fprintf_dym(stdout,"determinant u",u,"\n");
    fprintf_dyv(stdout,"determinant w_diag",w_diag,"\n");
    fprintf_dym(stdout,"determinant v",v,"\n");
  }

  for ( i = 0 ; i < dyv_size(w_diag) ; i++ )
    result *= dyv_ref(w_diag,i);

  free_dym(u);
  free_dym(v);
  free_dyv(w_diag);

  result = (is_determinant_positive(dm) ? 1.0 : -1.0) * fabs(result);

  return(result);
}

/**** WARNING THE FOLLOWING FUNCTION IS UNIMPLEMENTED AS YET.
      NUMERICAL RECIPES IS VERY UNCLEAR ABOUT PRECISELY WHAT
      ludcmp DOES, so our lu_decomp doesn't know what to do
*****/

bool attempt_lu_decomp(dym *a,dym *l,dym *u)
/* If result is true, the values of dym's l and u are set so that
   l is lower triangular, 
   u is upper triangular,
   a = l u

   Note that a is unchanged.
   It is legal to pass l or u as a (but not both, of course) in which
   case the contents of a are overwritten.

   If result is FALSE then the matrix A was singular.

   PRECONDISTIONS: a is a square matrix, and has same dimensions as l and u
*/
{
  int size = dym_rows(a);
  double **nra = am_malloc_2d_realnums(size+1,size+1);
  int *indx = am_malloc_ints(size+1);
  double d;
  int i,j;
  bool result;

  if ( size != dym_cols(a) ) my_error("lu_decomp: a should be square");

  assert_dym_shape(l,size,size,"attempt_lu_decomp (l)");
  assert_dym_shape(u,size,size,"attempt_lu_decomp (u)");

  for ( i = 0 ; i < size ; i++ )
    for ( j = 0 ; j < size ; j++ )
      nra[i+1][j+1] = dym_ref(a,i,j);

  result = numrec_ludcmp(nra,size,indx,&d);

  zero_dym(l);
  zero_dym(u);

  if ( result )
  {
/*
    for ( i = 0 ; i < size ; i++ )
      for ( j = 0 ; j < size ; j++ )
      {
        int true_i = indx[i+1];
        int true_j = j+1;
        if ( true_i > true_j )
*/
    fprintf_2d_realnums(stdout,"nra",nra,size+1,size+1,"\n");
    fprintf_ints(stdout,"indx",indx,size+1,"\n");
    printf("d = %g\n",d);
  }
  else
    printf("Singular matrix\n");

  my_error("not finished");

  return(result);
}

double **mk_nrecipes_matrix_from_dym(dym *d)
{
  double **tdarr = am_malloc_2d_realnums(d->rows+1,d->cols+1);
  int i,j;
  for ( i = 0 ; i < d->rows ; i++ )
    for ( j = 0 ; j < d->cols ; j++ )
      tdarr[i+1][j+1] = d->tdarr[i][j];
  return(tdarr);
}

void copy_nrecipes_matrix_to_dym(double **nrtdarr,dym *d)
{
  int i,j;
  for ( i = 0 ; i < dym_rows(d) ; i++ )
    for ( j = 0 ; j < dym_cols(d) ; j++ )
      dym_set(d,i,j,nrtdarr[i+1][j+1]);
}

void free_nrecipes_matrix(double **nrtdarr,dym *d)
{
  free_2d_realnums(nrtdarr,d->rows+1,d->cols+1);
}




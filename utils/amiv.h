/* *
   File:        amiv.h
   Author:      Andrew W. Moore
   Created:     Sat Apr  8 18:48:25 EDT 1995
   Updated:     26 Nov 96
   Description: Header for Integer and Boolean Dynamic vectors

   Copyright (C) 1995, Andrew W. Moore
*/

#ifndef AMIV_H
#define AMIV_H

#include "ambs.h"      /* Very basic operations */
#include "amma.h"      /* Fast, non-fragmenting, memory management */
#include "amar.h"      /* Obvious operations on 1-d arrays */

/* An ivec is an array of integers. You can read and write entries in the
   array and you can add or remove array entries. When compiled with AMFAST
   defined, it's all very fast and efficient. Without AMFAST defined it is
   slower but safer (run-time error checking happens).

Return the integer value in the i'th array element
(Precondition: 0 <= i < ivec_size(iv)) 
int ivec_ref(ivec *iv, int i);

Set the integer value in the i'th array element to be "value"
(Precondition: 0 <= i < ivec_size(iv)) 
void ivec_set(ivec *iv,int i,int value);

Increment the integer value in the i'th array element by "value"
(Precondition: 0 <= i < ivec_size(iv)) 
void ivec_increment(ivec *iv,int i,int value);

Return the number of elements in the array
int ivec_size(ivec *iv);

MAKE an ivec with the given number of elements. Eventually, unless 
you are prepared to cause a memory leak, you must free it with
free_ivec(iv). You MUST NOT try freeing with free(iv)
ivec *mk_zero_ivec(int size)

Increase the size of iv by 1, and store value in the new (rightmost) element.
void add_to_ivec(ivec *iv,int value)

Remove the i'th value, decreasing ivec's size by 1.
void ivec_remove(ivec *iv,int index);

Free iv and all its subcomponents. After this you must never access iv again.
void free_ivec(ivec *iv);
*/

typedef struct ivec_struct
{
  int ivec_code;
  int array_size;
  int size;
  int *iarr;
} ivec, *ivec_ptr;

int safe_ivec_ref(const ivec *iv, int i);
void safe_ivec_set(ivec *iv,int i,int value);
void safe_ivec_increment(ivec *iv,int i,int value);
int safe_ivec_size(const ivec *iv);

#ifdef AMFAST

#define ivec_ref(iv,i) ((iv)->iarr[i])
#define ivec_set(iv,i,v) ((iv)->iarr[i] = (v))
#define ivec_increment(iv,i,v) ((iv)->iarr[i] += (v))
#define ivec_size(iv) ((iv)->size)

#else

#define ivec_ref(iv,i) (safe_ivec_ref(iv, i))
#define ivec_set(iv,i,v) (safe_ivec_set(iv, i, v))
#define ivec_increment(iv,i,v) (safe_ivec_increment(iv, i, v))
#define ivec_size(iv) (safe_ivec_size(iv))

#endif

ivec *mk_ivec(int size);
void free_ivec(ivec *iv);
void fprintf_ivec(FILE *s,char *m1,const ivec *iv,char *m2);
void copy_iarr_to_ivec(int *iarr,int size,ivec *r_iv);
ivec *mk_ivec_from_iarr(int *iarr,int size);
void copy_ivec_to_iarr(ivec *iv, int *iarr);
int *mk_iarr_from_ivec(ivec *iv);

/* Makes an ivec of size end - start:
   { start , start+1 , .... end-2 , end-1 } */
ivec *mk_sequence_ivec(int start_value,int end_value);

/*
   Allocates and returns an ivec of size size in which ivec[i] = i
   ivec[0] = 0
   ivec[1] = 1
    .
    .
   ivec[size-1] = size-1
*/
ivec *mk_identity_ivec(int size);

void shuffle_ivec(ivec *iv);
void constant_ivec(ivec *iv,int value);
ivec *mk_constant_ivec(int size,int value);
void zero_ivec(ivec *iv);
ivec *mk_zero_ivec(int size);

void ivec_plus(ivec *iv_1, ivec *iv_2, ivec *r_d);
ivec *mk_ivec_plus(ivec *a,ivec *b);
void ivec_subtract(ivec *iv_1,ivec *iv_2,ivec *r_d);
ivec *mk_ivec_subtract(ivec *a,ivec *b);

void ivec_scalar_mult(ivec *iv,int scale,ivec *riv);
ivec *mk_ivec_scalar_mult(ivec *iv,int scale);
void copy_ivec(ivec *iv,ivec *r_iv);
/* Copies size elements from start */
void copy_ivec_subset(ivec *iv, int start, int size, ivec *r_iv);
ivec *mk_copy_ivec_subset(ivec *iv, int start, int size);
ivec *mk_copy_ivec(ivec *iv);
int num_of_given_value(ivec *iv,int value);
int num_zero_entries(ivec *iv);
int num_nonzero_entries(ivec *iv);
int ivec_min(const ivec *iv);
int ivec_max(const ivec *iv);
bool ivec_equal(const ivec *a,const ivec *b);
ivec *mk_ivec_1(int x0);
ivec *mk_ivec_2(int x0,int x1);
ivec *mk_ivec_3(int x0,int x1,int x2);
ivec *mk_ivec_4(int x0,int x1,int x2,int x3);
ivec *mk_ivec_5(int x0,int x1,int x2,int x3,int x4);
ivec *mk_ivec_6(int x0,int x1,int x2,int x3,int x4,int x5);

void ivec_remove(ivec *iv,int index); /* Reduces size by 1, removes index'th 
                                      element. All elements to right of
                                      delete point copied one to left.
                                      See comments in amdm.c more details */
void ivec_remove_last_element(ivec *iv); /* Reduce size by 1, remove 
                                        last rightmost elt */

/* Remove all occurences of val in iv. Does nothing
   if val doesn't appear in ivec. Keeps remaining
   iv elements in the same order. */
void ivec_remove_value(ivec *iv,int val);

int ivec_sum(ivec *iv);
bool equal_ivecs(ivec *iv1,ivec *iv2);
void ivec_malloc_report();
void add_to_ivec(ivec *iv,int val);
void add_to_sorted_ivec(ivec *siv, int val);

/* Increases the ivec length by one. Inserts val as the index'th element
   in the ivec and moves all items previously with array index greater
   than "index" one to the right.

   Thus iv_after[i] = iv_before[i] if i < index
        iv_after[i+1] = iv_before[i] = i >= index
        iv_after[index] = val */
void ivec_insert(ivec *iv,int index,int val);

/* Find least index i such that value = ivec_ref(iv,i).
  If not found, returns -1
*/
int find_index_in_ivec(const ivec *iv, int value);

/* Finds an index i such that value = ivec_ref(iv,i).
  If not found, returns -1.
  It will start searching at index "hint", and will continue searching
  in indices further and further away from hint, alternating between larger
  and smaller values than hint (modulu ivec_size).
*/
int find_index_in_ivec_hint(const ivec *iv, int value, int hint);

/* Finds leftmost index such that siv[index] > val. If none
   such, returns ivec_size(siv) */
int index_in_sorted_ivec(ivec *siv,int val);


bool is_in_ivec (ivec *iv, int value); /* Checks membership */

/* More sophisticated operations on ivecs, courtesy of jay... */

void ivec_sort(ivec *dv,ivec *r_dv);
ivec *mk_ivec_sort(ivec *iv);

ivec *mk_ivec_union(ivec *v1, ivec *v2);

/*Pre: ivecs are ordered
*/
ivec *mk_ivec_union_ordered(ivec *v1, ivec *v2);
ivec *mk_ivec_diff_ordered(ivec *v1, ivec *v2);

ivec *mk_append_ivecs(ivec *a,ivec *b);

/* And now, the same for arrays of strings. Note that strings internally
   are always copied around and dynamically alolocated and deallocated.
*/

typedef struct string_array_struct
{
  int string_array_code;
  int size;
  int sarr_size;
  char **sarr;
} string_array, *string_array_ptr;

string_array *mk_string_array(int size);
void free_string_array(string_array *sar);
string_array *mk_copy_string_array(string_array *sa);
string_array *mk_string_array_1(char *x1);
string_array *mk_string_array_2(char *x1,char *x2);
string_array *mk_string_array_3(char *x1,char *x2,char *x3);
string_array *mk_string_array_4(char *x1,char *x2,char *x3,char *x4);
string_array *mk_string_array_5(char *x1,char *x2,char *x3,char *x4,char *x5);
void sort_string_array(string_array *sa);
string_array *mk_sort_string_array(string_array *sa);
void fprintf_string_array(FILE *s,char *m1,string_array *sar,char *m2);
void fprintf_string_array_contents(FILE *s,string_array *sar);
void fprintf_string_array_contents_on_one_line(FILE *s,string_array *sar);


char *string_array_ref(string_array *sar, int i);
void string_array_set(string_array *sar,int i,char *value);
int string_array_size(string_array *sar);
string_array *mk_string_array_from_array(char **sarr,int size);
void add_to_string_array(string_array *sa,char *string);
/* Following function adds an element at the pos location.  pos can be
   at the end of the array */
void insert_in_string_array(string_array *sa, int pos, char *string);

/* Following function is synonym for
   add_to_string_array. Not the preferred name, retained
   for compatibility
*/
void string_array_add(string_array *sa,char *string);
string_array *mk_broken_string_using_seppers(char *string,char *seppers);
/* whitespace is removed, and each word is placed as a separate entry
       int the string array.

       If string is NULL or "" or has only whitespace, a string_array of
       length 0 is returned.
    */
string_array *mk_broken_string(char *string);
    
/*  As above except removes all double quote marks before doing the splitting.
*/
string_array *mk_broken_quoteless_string(char *string);

/* Removes all occurences of c from string. If input was a string of
   length m containing n 'c' characters, then after it's a string of
   length m-n. */
char *mk_string_without_character(char *string,char c);
char *mk_quoteless_string(char *s);

/* Puts the last n elements of sa into one long string, in which each element
   is separated by the character sep.  If !sep, then they are not separated
   by anything.*/
char *mk_string_from_last_n_with_separator(string_array *sa,int n,char sep);
char *mk_string_from_last_n(string_array *sa,int n);
char *mk_string_from_string_array_with_separator(string_array *sa,char sep);
char *mk_string_from_string_array(string_array *sa);

char *mk_string_from_line(FILE *s);
string_array *mk_string_array_from_line(FILE *s);
void string_array_malloc_report();

int find_index_in_string_array(string_array *sa,char *string);

int caseless_find_index_in_string_array(string_array *sa,char *string);
bool string_array_member(string_array *sa,char *string);

void string_array_remove(string_array *sa,int index); 
                                   /* Reduces size by 1, removes index'th 
                                      element. All elements to right of
                                      delete point copied one to left.
                                      See comments in amdm.c more details */
void string_array_remove_last_element(string_array *sa); 
                                     /* Reduce size by 1, remove 
                                        last rightmost elt */

/* PRE: 0 <= start_index <= end_index <= size
   Returns a new string array of size end_index - start_index
   consisting of 
     { sa[start_index] , sa[start_index+1] ... sa[end_index-1] }
*/
string_array *mk_string_array_segment(string_array *sa,
                                      int start_index,int end_index);

typedef struct string_matrix
{
  int array_size;
  int rows;
  int cols;
  string_array **sas;
} string_matrix, *string_matrix_ptr;

string_matrix *mk_string_matrix(int rows,int cols);
int string_matrix_rows(string_matrix *sm);
int string_matrix_cols(string_matrix *sm);
char *string_matrix_ref(string_matrix *sm,int row,int col);
void string_matrix_set(string_matrix *sm,int row,int col,char *value);
void fprintf_string_matrix(FILE *s,char *m1,string_matrix *sm,char *m2);
void free_string_matrix(string_matrix *sm);
void string_matrix_add_row(string_matrix *sm);

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
string_array *mk_tabular_string_array(string_matrix *sm,string_array *seps);
void render_string_matrix(FILE *s,char *comment,string_matrix *sm);

void string_matrix_real_set(string_matrix *sm,int row,int col,double value);
void string_matrix_row_from_broken_string(
    string_matrix *sm,
    int row,
    char *row_string
  );

bool ivec_weakly_dominates(ivec *dx,ivec *dy);
string_array *mk_string_array_from_ivec(ivec *iv);
char *mk_string_from_ivec(ivec *iv);

int find_index_in_sorted_string_array(string_array *sa,char *s);
/* place s into string array whether or not a duplicate already exists */
void insert_in_sorted_string_array(string_array *sa,char *s);
/* place s into string array only if it does not already exist there */
void maybe_insert_in_sorted_string_array(string_array *sa,char *s);

/* Ivec Arrays */

typedef struct ivec_array
{
  int size;
  int array_size;
  ivec **array;
} ivec_array;

#ifdef AMFAST
#define ivec_array_ref(iva,i) ((iva)->array[i])
#else
#define ivec_array_ref(iva,i) safe_ivec_array_ref(iva,i)
#endif

// Added by Artur
ivec_array *mk_transpose_ivec_array( ivec_array *iva );

/*Added by Dan: Something I've wanted for a LONG time!*/
#define ivec_array_ref_ref(iva,i,j) ivec_ref(ivec_array_ref(iva,i),j)
#define ivec_array_ref_set(iva,i,j,x) ivec_set(ivec_array_ref(iva,i),j,x)

ivec_array *mk_empty_ivec_array (void);
void add_to_ivec_array (ivec_array *ivecarr, ivec *this_ivec);
int ivec_array_size (ivec_array *ivecarr);
ivec *safe_ivec_array_ref (ivec_array *ivecarr, int index);
void ivec_array_set (ivec_array *iva, int index, ivec *iv);
void fprintf_ivec_array (FILE *s, char *m1, ivec_array *ivecarr, char *m2);
void free_ivec_array (ivec_array *ivecarr);
ivec_array *mk_copy_ivec_array (ivec_array *ivecarr);
ivec_array *mk_array_of_zero_length_ivecs (int size);
void ivec_array_remove(ivec_array *iva,int index);

void add_to_error_message(char **errmess, char *new_mess);
void prepend_error_message(char **errmess, char *new_mess);

#endif /* #ifndef AMIV */

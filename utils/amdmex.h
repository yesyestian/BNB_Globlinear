/* 
   File:        amdmex.h
   Author:      Andrew W. Moore
   Created:     Oct 11, 1996
   Updated:     8 Dec 86
   Description: Extensions and advanced amdm stuff

   Copyright 1996, Schenley Park Research

   This file contains advanced utility functions involving dyvs dyms and
   ivecs. It never accesses the data structures directly, so if the
   underlying representation of dyms and dyvs changes these won't need to.

   The prototypes of these functions used to be declared at the end of amdm.h.
   Now it's amdmex.h
*/

#ifndef AMDMEX_H
#define AMDMEX_H

#include "amiv.h"
#include "amdym.h"
#include "svd.h"

typedef struct dyv_array
{
  int size;
  int array_size;
  dyv **array;
} dyv_array;

typedef struct dym_array
{
  int size;
  int array_size;
  dym **array;
} dym_array;

/* x := x with y appended on the end */
void append_to_dyv(dyv *x,dyv *y);

/* Return x with y appended on the end */
dyv *mk_dyv_append(dyv *x,dyv *y);

/* x := x with y appended on the end */
void append_to_ivec(ivec *x,ivec *y);

/* Return x with y appended on the end */
ivec *mk_ivec_append(ivec *x,ivec *y);

/* Makes a dyv of the same size as lo and hi (which must
   both be the same size such that
   dyv_ref(result,i) is uniformly randomly distributed between
   dyv_ref(lo,i) and dyv_ref(hi,i)
*/
dyv *mk_dyv_range_random(dyv *lo, dyv *hi);

/* Forall (i,j) such that ilo <= i < ihi
                          jlo <= j < jhi (note the strict inequality on right)

   We do a[i][j] += delta
*/
void dym_increment_block(dym *a,int ilo,int jlo,int ihi,int jhi,double delta);

/* Forall i such that ilo <= i < ihi (note the strict inequality on right)

   We do a[i] += delta
*/
void dyv_increment_block(dyv *a,int ilo,int ihi,double delta);


/* After calling, ivec_ref(iv,k) points to the k'th smallest
   element in dv (ties handles arbitrarily) */
void indices_of_sorted_dyv(const dyv *dv,ivec *iv);

/* After calling, ivec_ref(result,k) points to the k'th smallest
   element in dv (ties handles arbitrarily) */
ivec *mk_indices_of_sorted_dyv(const dyv *dv);

ivec *mk_indices_of_sorted_ivec(ivec *v); // Artur

void sorted_eigens_of_spd_dym(dym *d,dym *evectors,dyv *evalues);
dyv *mk_sorted_eigens_of_spd_dym(dym *d);
void random_unit_dyv(dyv *dv);
bool dyv_weakly_dominates(const dyv *dx, const dyv *dy);
bool dyv_equal(dyv *dx,dyv *dy);
bool dym_weakly_dominates(dym *dx,dym *dy);
bool dym_equal(dym *dx,dym *dy);
dyv *mk_normalize_dyv(dyv *d);
void normalize_dyv(dyv *src,dyv *dest);

/* Makes a random dyv with magnitude 1.0 */
dyv *mk_random_unit_dyv(int size);
void dyv_clip_in_unit(dyv *dv);
dyv *mk_dyv_clip_in_unit(dyv *dv);
void random_rotation(dym *d);
/* Declared but not defined - sir 8/6/2000
void test_random_rotation(int argc,char *argv[]); */
void fprintf_oneline_dyv(FILE *s,const char *m1, const dyv *d, const char *m2);

/* Makes a random subset of iv of size "k" */
ivec *mk_random_ivec_subset(ivec *iv,int k);

/* add_to_me := add_to_me appended to to_be_added */
void append_to_ivec(ivec *add_to_me,ivec *to_be_added);

ivec *mk_ivec_from_dyv(dyv *d);
dyv *mk_dyv_from_ivec(ivec *iv);
ivec *mk_ivec_from_args(char *key,int argc,char *argv[],ivec *def);
void copy_dym_row(dym *source,int source_row,dym *dest,int dest_row);
void copy_dym_col(dym *source,int source_col,dym *dest,int dest_col);
void swap_dym_rows(dym *dm,int i,int j);
void swap_dym_cols(dym *dm,int i,int j);
double dyv_sum(const dyv *dv);

/* In the following 3 functions, rows may be NULL denoting USE ALL ROWS */
dyv *mk_sum_of_dym_rows(dym *dm,ivec *rows);
dyv *mk_mean_of_dym_rows(dym *dm,ivec *rows);
dyv *mk_sdev_of_dym_rows(dym *dm,ivec *rows);

ivec *mk_ivec_from_dym_col(dym *a,int col);
ivec *mk_ivec_from_dym_row(dym *a,int row);

void diag_times_dyv(dyv *a,dyv *b,dyv *x);
dyv *mk_diag_times_dyv(dyv *a,dyv *b);

dyv *mk_sum_of_dym_cols(dym *dm);
dyv *mk_mean_of_dym_cols(dym *dm);

dyv_array *mk_empty_dyv_array();
void add_to_dyv_array(dyv_array *da, const dyv *dv);
int dyv_array_size(const dyv_array *da);
dyv *safe_dyv_array_ref(const dyv_array *da,int index);

dyv_array *mk_array_of_zero_length_dyvs(int size);

#ifdef AMFAST
#define dyv_array_ref(dva,i) ((dva)->array[i])
#else
#define dyv_array_ref(dva,i) safe_dyv_array_ref(dva,i)
#endif

void dyv_array_set(dyv_array *dva, int index, const dyv *dv);
void fprintf_dyv_array(FILE *s,char *m1,dyv_array *da,char *m2);
void free_dyv_array(dyv_array *da);
dyv_array *mk_copy_dyv_array(const dyv_array *da);
void dyv_array_remove(dyv_array *dva,int index);
bool is_in_dyv(dyv *dv,double value, double error);

dyv *mk_sum_of_dyv_array(const dyv_array *da);

dym_array *mk_empty_dym_array();
void add_to_dym_array(dym_array *da,dym *dm);
int dym_array_size(dym_array *da);
dym *dym_array_ref(dym_array *da,int index);
void fprintf_dym_array(FILE *s,char *m1,dym_array *da,char *m2);
void free_dym_array(dym_array *da);
dym_array *mk_copy_dym_array(dym_array *da);

dyv *mk_dyv_from_string_array_with_error_message(string_array *sa,
                                                 char *format,
                                                 char **err_mess);
dyv *mk_dyv_from_string(char *string,char *format);
dyv *mk_dyv_from_string_with_error_message(char *string,char *format,char **err_mess);
void mk_io_dyvs_from_string(char *string,char *format,dyv **r_in_dyv,dyv **r_out_dyv);
char *mk_string_from_dyv(dyv *d);
dyv *mk_midpoint_dyv(dyv *a,dyv *b);
int find_index_of_kth_smallest(const dyv *x,int k);
double dyv_kth_smallest(const dyv *d,int k);
double dyv_median(const dyv *d);
void save_dym_to_file(FILE *s, dym *m);

void save_dyv_to_file(FILE *s, dyv *v);
dyv *mk_dyv_from_file(FILE *s, char **r_errmess);
void save_ivec_to_file(FILE *s, ivec *v);
ivec *mk_ivec_from_file(FILE *s, char **r_errmess);
void save_dym_array_to_file(FILE *s, dym_array *ma);
dym_array *mk_dym_array_from_file(FILE *s, char **r_errmess);
dyv *mk_dyv_from_string_array( string_array *sa, char *format );

/* Returns a dyv_array of size "num_dyvs" in which the i'th
   element (forall 0 <= i < num_dyvs) is a vector of
   size "dyv_size" in which each element is 0.0 */
dyv_array *mk_dyv_array_of_zeroed_dyvs(int num_dyvs,int dyv_size);

/* Returns a dym_array of size "num_dyms" in which the i'th
   element (forall 0 <= i < num_dyms) is a matrix of
   size "dym_size x dym_size" in which each element is 0.0 */
dym_array *mk_dym_array_of_zeroed_dyms(int num_dyms,int dym_size);

/* Returns a dym_array of size "num_dyms" in which the i'th
   element (forall 0 <= i < num_dyms) is a matrix of
   size "dym_size_r x dym_size_c" in which each element is 0.0 */
dym_array *mk_dym_array_of_zeroed_nonrect_dyms(int num_dyms,int dym_size_r,
											   int dym_size_c);

string_array *mk_string_array_from_argc_argv(int argc,char *argv[]);
string_array *mk_string_array_from_stream_tokens(FILE *s);
string_array *mk_string_array_from_file_tokens(char *filename);

void make_argc_argv_from_string_array(string_array *sa,int *actual_argc,
                                      char ***actual_argv);

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
                           int *actual_argc,char **r_actual_argv[]);

void free_loaded_argc_argv(int argc,char **argv);

/* rows may be NULL denoting "use all datapoints" */
void make_vector_limits(dym *x,ivec *rows,dyv **xlo,dyv **xhi);
void set_vector_limits_sensible(dyv *xlo,dyv *xhi);

void make_vector_limits_from_dyv_array(const dyv_array *da, dyv **xlo,dyv **xhi);

/* rows may be NULL denoting "use all datapoints" */
void make_sensible_vector_limits(dym *x,ivec *rows,dyv **xlo,dyv **xhi);

/* A line is interesting if its not all white space and
the leftmost non whitespace character isnt # */
bool line_string_is_interesting(char *line_string);

/* Searches the file for the next line that isn't all whitespace and
   that doesn't have # as its first non-whitespace character. 

   If no-such line before file-end, returns NULL */
char *mk_next_interesting_line_string(FILE *s,int *line_number);

/* As above excepts breaks resulting line into a string array of tokens... */
string_array *mk_next_interesting_line(FILE *s,int *line_number);

bool file_exists(char *fname);
void save_dyv_to_file_plain(char *fname,dyv *x,char **r_errmess);
void remove_file(char *fname,char **r_errmess);
void execute_command(char *exec_string,char **r_errmess);
dyv *mk_dyv_from_file_plain(char *fname,int size,char **r_errmess);

/********* DATFILE PARSING UTILITIES ************/

/* int linestyle formats (we recommend AUTO_FORMAT)....

   The following constants determine how lines are read from
   a datafile. COMMA_FORMAT expects commas between each line on
   the datafile. WHITESPACE expects one or more spaces between 
   each item. And AUTO_FORMAT will accept both commas and whitespace
   as separators. */
#define COMMA_FORMAT      0
#define WHITESPACE_FORMAT 1
#define AUTO_FORMAT       2

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

bool contains_a_number(string_array *sa);
bool line_has_unquoted_comma(char *string);
string_array *mk_parse_data_line(char *string,int line_format);
string_array *mk_next_tokens(FILE *s,int *line_number,int line_format);
string_array *mk_default_attribute_names(int num_atts);

void make_attnames_and_dym_from_filename(char *filename,int argc,char *argv[],
					 string_array **r_attnames,
                                         dym **r_x,char **r_errmess);

dym *mk_dym_from_filename(char *filename,char **r_errmess);
dym *mk_dym_from_filename_simple(char *filename);
void save_attnames_and_dym(FILE *s,string_array *attnames,dym *x);
void save_dym_to_filename(char *filename,dym *x);

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

bool is_sivec(const ivec *iv);

/* Makes { 0 , 1 , ... size-1 } */
ivec *mk_identity_sivec(int size);

/* Returns number of elements in sivec */
int sivec_size(const ivec *siv);

/* Returns the minimum value in sivec. 
   Time cost: Constant */
int sivec_min(const ivec *siv);

/* Returns the maximum value in sivec. 
   Time cost: Constant */
int sivec_max(const ivec *siv);

/* Adds the element while maintaining legal siveckiness.
   (If element already there, no change)
   Time cost: O(size) */
void add_to_sivec(ivec *siv,int value);

ivec *mk_add_to_sivec(ivec *siv,int value);

/* Returns -1 if the value does not exist in siv.
   Else returns index such that
      value == ivec_ref(siv,value) 
  Time cost: O(log(size)) */
int index_in_sivec(ivec *siv,int value);

/* Returns true iff siv contains value
   Time cost: O(log(size)) */
bool is_in_sivec(ivec *siv,int value);

/* Does nothing if value is not in siv.
   If value is in siv, the sivec is updated to
   represent siv \ { value } */
void sivec_remove_value(ivec *siv,int value);
  
/* Returns answer to A subset-of B?
   Returns true if and only if the set of integers in a is
   a subset of the set of integers in b */
bool sivec_subset(const ivec *siva,const ivec *sivb);

bool sivec_equal(const ivec *siva,const ivec *sivb);

/* Returns TRUE iff A is a subset of B and A != B */
bool sivec_strict_subset(const ivec *siva,const ivec *sivb);

ivec *mk_sivec_union(const ivec *siva,const ivec *sivb);

/* Returns A \ B.
   This is { x : x in A and x not in B } */
ivec *mk_sivec_difference(const ivec *siva,const ivec *sivb);

ivec *mk_sivec_intersection(const ivec *siva,const ivec *sivb);

/* Returns TRUE iff A intersect B is empty. O(size) time */
bool sivec_disjoint(ivec *siva,ivec *sivb);

ivec *mk_sivec_from_ivec(ivec *iv);

ivec *mk_ivec_from_string(char *s);

/* Turns a space separated string into a sivec.
   Example: "3 1 4 1 5 9" ====> { 1 , 3 , 4 , 5 , 9 } */
ivec *mk_sivec_from_string(char *s);

ivec_array *mk_array_of_empty_ivecs(int size);

int index_of_longest_ivec(ivec_array *iva);

/********** Utilities added by Artur. But is "mk_ranks_fast" just
            the same as "mk_indices_of_sorted_dyv"??? (-AWM) ***/
/********** Artur answers: not neccesarily, because ranks are real valued
				(indices are integers). And ranks are real valued because:
				1) it allows for resolving ties among the original dyv components
				2) the resultant dyv (of ranks) can be straightforwardly used
				    by the correlation stuff in petgui/correlate.c
			  Nota bene mk_ranks_fast calls mk_indices_of_sorted_dyv... */

//----------------------- returns a random permutation of v
void shuffle_dyv(dyv *v);

/*************
This does the same as above, but may be faster in case of large datasets
with very many ties in the processed attributes.
**************/
dyv *mk_ranks_fast( dyv *source );

/***************** sosarray ***************/

bool string_less(char *s1,char *s2);
bool string_greater(char *s1,char *s2);
bool string_leq(char *s1,char *s2);
bool string_geq(char *s1,char *s2);

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

bool is_sosarray(string_array *sa);

/* Returns number of elements in sosarray */
int sosarray_size(string_array *sosarr);

/* If sosarr has 0 elements returns 0
   If value > string_array_max(sosarr) returns size
   If value <= string_array_min(sosarr) returns 0
   Else returns index such that
      value <= string_array_ref(sosarr,index) 
      string_array_ref(sosarr,index-1) < value
      
   It returns the value such that string_array_insert(sa,index,value)
   would represent the set with value added to sa (assuming value
   wasn't already in sa). */
int find_sosarray_insert_index(string_array *sosarr,char *string);

/* Adds the element while maintaining legal sosarraykiness.
   (If element already there, no change)
   Time cost: O(size) */
void add_to_sosarray(string_array *sosarr,char *string);

/* Returns -1 if the string does not exist in sosarr.
   Else returns index such that
      string == string_array_ref(sosarr,string) 
  Time cost: O(log(size)) */
int index_in_sosarray(string_array *sosarr,char *string);

/* Returns true iff sosarr contains string
   Time cost: O(log(size)) */
bool is_in_sosarray(string_array *sosarr,char *string);

void sosarray_remove_at_index(string_array *sosarr,int index);

/* Does nothing if string is not in sosarr.
   If string is in sosarr, the sosarray is updated to
   represent sosarr \ { string } */
void sosarray_remove_string(string_array *sosarr,char *string);
  
/* Returns answer to A subset-of B?
   Returns true if and only if the set of integers in a is
   a subset of the set of integers in b */
bool sosarray_subset(string_array *sosarra,string_array *sosarrb);

bool equal_string_array(string_array *sa1,string_array *sa2);

bool sosarray_equal(string_array *sosarra,string_array *sosarrb);

/* Returns TRUE iff A is a subset of B and A != B */
bool sosarray_strict_subset(string_array *sosarra,string_array *sosarrb);

string_array *mk_sosarray_union(string_array *sosarra,string_array *sosarrb);

/* Returns A \ B.
   This is { x : x in A and x not in B } */
string_array *mk_sosarray_difference(string_array *sosarra,string_array *sosarrb);

string_array *mk_sosarray_intersection(string_array *sosarra,string_array *sosarrb);

/* Returns TRUE iff A intersect B is empty. O(size) time */
bool sosarray_disjoint(string_array *sosarra,string_array *sosarrb);

string_array *mk_sosarray_from_string_array(string_array *sa);

string_array *mk_string_array_from_string(char *s);

/* Turns a space separated string into a sosarray.
   Example: "3 1 4 1 5 9" ====> { 1 , 3 , 4 , 5 , 9 } */
string_array *mk_sosarray_from_string(char *s);

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
string_matrix *mk_ls_style_string_matrix_given_cols(string_array *name,int cols);

/* Returns the max string length in sa */
int string_array_max_length(string_array *sa);

/* Makes an 'LS' style string matrix cleverly designed so that when
   printed it uses less than "max_chars" characters per line. (it auto-chooses
   and sizes the columns) */
string_matrix *mk_ls_style_string_matrix(string_array *names,int max_chars);

/* Prints the contents of string_array cleverly much in the same way
   that "short ls" in unix displays filenames. It's cleverly designed so that when
   printed it uses less than "max_chars" characters per line. (it auto-chooses
   and sizes the columns) */
void display_names(FILE *s,string_array *names);


#endif /* #ifndef AMDMEX_H */



/*
   File:        stats.h
   Author:      Andrew W. Moore
   Created:     Sat Jun 17 10:17:26 EDT 1995
   Description: Header for Cached cumu denisty functions for fast statistics

   Copyright (C) 1995, Andrew W. Moore
*/


#ifndef stats_H
#define stats_H
#include "amut.h"

typedef struct integ
{
  dyv *integral;
  double xlo;
  double xhi;
  double parameter;
  double constant;
} integ, *integ_ptr;

typedef struct minteg
{
  int integ_array_size;
  int integ_size;
  double (*h)(double parameter,double constant,double x);
  double (*compute_constant)(double parameter);
  double (*xlo)(double parameter,double constant);
  double (*xhi)(double parameter,double constant);
  double parameter_lo;
  double parameter_hi;
  integ **integ_array;
}  minteg;

/**** erf and lgamma don't exist on msdos. Here are our re-implementations,
      courtesy of Numerical Recipes in C.  Programmers should call these
      functions rather than erf and lgamma.
      - Mary Soon Lee, 2 Nov 1995
****/

/* This should be called instead of lgamma on Windows NT platforms.
   The code was taken from Numerical Recipes in C (where it is called
   gammln) , amended to use doubles not floats. 

   N.B. I have added a check to make sure the argument is in the correct
   range.
*/
double am_lgamma(double xx);
      
/* This should be called instead of erf on Windows NT platforms.  The
   code was taken from Numerical Recipes in C, amended to use doubles
   not floats. */
double am_erf(double x);



double gamma_pdf(double x,double alpha,double beta);

double gamma_cdf(double x,double alpha,double beta);

double gamma_cdf_inv(double prob,double alpha,double beta);

double t_pdf(double x,double mu,double std_dev,double dof);

double t_cdf(double x,double mu,double std_dev,double dof);

double t_cdf_inv(double prob,double mu,double std_dev,double dof);

void close_statistics();

void stats_malloc_report();
/* Declared but not defined - sir 8/6/2000
void stats_main(int argc,char *argv[]); */

void fprintf_minteg(FILE *s,char *m1,minteg *mit,char *m2);

/* Declared but not defined - sir 8/6/2000
void t_main(int argc,char *argv[]); */

double beta_that_forces_requested_median(double alpha,double median);

/* 
   PRE: size of actual_dist is same as size of hypothesized_dist.

   Given two distributions represented as histograms 
   (actual_dist and hypothesized_dist), how much evidence is there that they are
   from the same distribution? 
   Note that these things must be counts. Each element of actual_dist must
   be an integer. Each element of hypothesized_dist may be non-integer
   because we're talking expected counts there.

   The prob returned by this function answers that question
   using a standard chi-squared test. If it is low (e.g. < 0.05), then it is
   unlikely that they are the same. 

   The "dof" parameter is the mysterious "Degrees Of Freedom" that haunts
   any use of the word "Chi". 
   
       If it is possible for any entry in the dist
       to take any value, then set dof==size.

       If the sum of values is constrained to a certain value
       then set dof==size-1.

       If there are more constraints than that, then subtract
       more from size.
*/
double chi_squared_prob(dyv *actual_dist,dyv *hypothesized_dist,int dof);

/******************** Now some utilities that can help with k-fold-crossvalidation.. *******/

/* Define an ivec as being a natset if and only if it 
   contains a subset of the natural numbers in a
   strictly increasing order. */

/* PRE: a and b are natsets.
   PRE: b is a subset of a

   POST: returns a natset c such that C intersect B = empty
                             C union B     = A

       (i.e. n is in result if and only if n is in
             A but not in B)

   This function is useful if you have a set of rows (record numbers)
   to use as a test set (in b) and you have all the rows for the test
   and train set in a, and you want to get the set of remaining
   rows to use as a training set. */
ivec *mk_ivec_set_difference(ivec *a,ivec *b);

/* If you've got a dataset where you are training and testing
   using the set of rows specified in "train_and_test_rows" and you
   wish to do k-fold cross-validation (k == num_folds) then what 
   subset of train_and_test_rows should you use as the train set
   fold_num'th fold? And which should be the test set? 
   This function tells you.
   It uses a deterministic shuffle and then selects the given rows. 

   PRE: train_and_test_rows is a natset (see above) (it means
        strictly increasing set of natural numbers) (in which case
        num_rows is ignored).

        OR train_and_test_rows may be NULL, indicating "use all
        rows in the set { 0 , 1 , ... num_rows-1 }

   POST: *r_train_set is the fold_num'th training set
         *r_test_set is the fold_num'th test set.
         (Note that the union of these two is the original set of rows).
         The resulting sets are both returned in strictly increasing
         order.

    Note that if you called this function repeatedly with the same number
    of folds, but you varied fold_num from 0 up to num_folds on each
    call, then the union of all the resulting test_sets would contain the
    same set as train_and_test_rows. But all pairs of resulting test_sets 
    would have an empty intersection.
*/
void make_kfold_rows(ivec *train_and_test_rows,int num_rows,int num_folds,int fold_num,
                     ivec **r_train_rows,ivec **r_test_rows);

/* If you've got a dataset where you are training and testing
   using the set of rows specified in "train_and_test_rows" and you
   wish to make a test-set of size "num_test_rows" and a training-set
   of size "num_training_row (= num_tain_and_test - num_test)"
   this returns the training and test rows you should use, obtained by
   shuffling train_and_test.

   PRE: train_and_test_rows is a natset (see above) (it means
        strictly increasing set of natural numbers) (in which case
        num_rows is ignored).

        OR train_and_test_rows may be NULL, indicating "use all
        rows in the set { 0 , 1 , ... num_rows-1 }

   POST: *r_train_set is the training set
         *r_test_set is the test set.
         (Note that the union of these two is the original set of rows).
         The resulting sets are both returned in strictly increasing
         order.
*/
void make_train_and_test_rows(ivec *train_and_test_rows,
			      int num_rows,int num_test_rows,
			      ivec **r_train_rows,ivec **r_test_rows);

/* Makes a dym consisting of a subset of the rows in x. The members of
   of the subset are those rows mentioned in "rows".
   Result will this have "ivec_size(rows)" rows and dym_cols(x) columns */
dym *mk_dym_from_subset_of_rows(dym *x,ivec *rows);

/* Creates two new dyms, both with the same number of columns as the
   input x. *r_test will have "num_test_rows" rows and *r_train will
   have "dym_rows(x) - num_test_rows" rows. *r_test will consist of 
   "num_test_rows" rows of x selected at random without replacement. 
   *r_train will contain the rest. The rows will appear in the same order
   as they did in the original. */
void break_dym_into_train_and_test(dym *x,int num_test_rows,
				   dym **r_train,dym **r_test);

#endif /* #ifndef stats_H */


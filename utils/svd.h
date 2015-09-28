
/*
   File:        svd.h
   Author:      Andrew W. Moore
   Created:     Thu Sep 15 21:01:12 EDT 1994
   Updated:     amdm was split into amdyv, amdym and svd by Frank Dellaert, Aug 14 1997
   Description: Singular Value Decomposition Code

   Copyright (c) 1996,1997, Schenley Park Research
*/

#ifndef SVD_H
#define SVD_H

#include "amdym.h"

bool sing_val_decomp( 
    double **tdarr,
    int rows,
    int cols,
    double **u,
    double *w,
    double **v
  );
/*
   Takes as input a matrix represented by a 2-d array "tdarr".
   The matrix has "rows" rows and "cols" columns, thus the legal entries in
   it are
      tdarr[i][j] where 0 <= i < rows and 0 <= j < cols.

   It computes the singular value decomposition and stores the result in
   u, v and w. (decribed below).

   PRECONDITIONS: tdarr is an array of vectors of doubles. It contains
                  "rows" vectors of doubles. Each vector of doubles contains
                  "cols" elements.
     u is a similar 2-d array with exactly the same dimensions.
     v is a similar structure 2-d array, except its dimensions
        are: "cols" rows AND "cols" columns.

     w is a 1-d array (a vecoir) of doubles with "cols" elements

     IMPORTANT: When this routine is called none of tdarr, u ,v or w
     must share any of the same memory.

  POST-CONDITIONS:

     Returns TRUE if succeeded in performing the computation.

     Returns FALSE, and sets everything to zero if failed.

     If successful:
     Let W = Diag(w) = matrix with W[ii] = w[i], W[ij] = 0 if i != j.

     Then on exit from this routine we will have
        tdarr = u * W * v^T where v^T is the transpose of v and * is matrix *'s
 
     furthermore u will be orthoganal and v will be orthogonal, meaning
     u^T u = Identity, thus, if square, U^-1 = U^T (ditto for V).

     This in turn means that if tdarr is a square matrix then

       tdarr^-1 = v * DW * u^T where

      DW[ii] = 1 / w[i] and DW[ij] = 0 if i != j

     See pages 59-70 of Numerical Recipes in C for further info.
*/

/*

void dym_svd(dym *d,dym *u,dyv *w_diagonal,dym *v);
Places the SVD of d into u, w_diagonal, and v.
If d is a R-row, C-column matrix, then..
   u should be R-row , C-colum.
   w should be sized C
   v should be sized C-row , C-column.

     The SVD of the matrix d has the property that:
     
     Let W = Diag(w_diag) = matrix with W[ii] = w[i], W[ij] = 0 if i != j.

     Then on exit from this routine we will have
        d = u * W * v^T where v^T is the transpose of v and * is matrix *'s
 
     furthermore u will be orthoganal and v will be orthogonal, meaning
     u^T u = Identity, thus, if square, U^-1 = U^T (ditto for V).

     This in turn means that if d is a square matrix then

       d^-1 = v * DW * u^T where

      DW[ii] = 1 / w[i] and DW[ij] = 0 if i != j so DW = Diag(w_diag)^-1

     See pages 59-70 of Numerical Recipes in C for further info.
 /

void make_svd_components(dym *d,dym **u,dyv **w_diagonal,dym **v);
/ 
   The same as dym_svd except u, w_diag and v are dynamically created.
   The caller of this routine only suppiles pointer to pointers to the
   three structures.

    Example:

      dym *d; ....

      dym *u,*v;
      dyv *w_diag;
      make_svd_components(d,&u,&w_diag,&v);
 /

void dym_svd_backsub(dym *u,dyv *w_diag,dym *v,dyv *b,dyv *x);
/ 
   If u, w_diag and v were obtained via SVD from the matrix A,
   then this function solves A x = b, given a known A and B and
   stores the result in r_x. If the equation has too few unknowns
   then it is solved in a least-squares sense. That's the beauty of SVD
 /

*/

void dym_svd(dym *d,dym *u,dyv *w_diagonal,dym *v);

void make_svd_components(dym *d,dym **u,dyv **w_diagonal,dym **v);

void dym_svd_backsub(dym *u,dyv *w_diag,dym *v,dyv *b,dyv *x);

#endif /* #ifndef SVD_H */

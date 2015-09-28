
/*
   File:        svd.c
   Author:      Andrew W. Moore
   Created:     Thu Sep 15 21:01:13 EDT 1994
   Updated:     amdm was split into amdyv, amdym and svd by Frank Dellaert, Aug 14 1997
   Description: Singular Value Decomposition Code

   Copyright (c) 1996,1997, Schenley Park Research
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "svd.h"

#include "ambs.h"      /* Very basic operations */
#include "amma.h"      /* Fast, non-fragmenting, memory management */
#include "amar.h"      /* Obvious operations on 1-d arrays */

   /*************** FIRST SOME STUFF COPIED OUT OF ***************/
     /*************** NUMERICAL RECIPES IN `C'. ***************/

#define NEW_VERSION

#ifdef NEW_VERSION

#define MAX_SVD_ITERS 120

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

static bool dym_copy_of_svdcmp(double **a, int m, int n, double w[], double **v)
{
  double pythag(double a, double b);
  int flag,i,its,j,jj,k;
  int l = 0;  /* To prevent "possible uninitialized variable" warning */
  int nm = 0; /* To prevent "possible uninitialized variable" warning */
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  
  rv1=AM_MALLOC_ARRAY(double,1+n);
  g=scale=anorm=0.0;
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<=m;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=IMIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=MAX_SVD_ITERS;its++) {   /* change */
      flag=1;
      for (l=k;l>=1;l--) {
	nm=l-1;
	if ((double)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((double)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((double)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=1;j<=m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == MAX_SVD_ITERS)    /* change */
      {
        AM_FREE_ARRAY(rv1,double,n+1);
        fprintf(stderr,"*** amdm.c, warning, SVD did not converge\n");
        return(FALSE);
      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  AM_FREE_ARRAY(rv1,double,n+1);
  return(TRUE);
}
#endif

#ifdef OLD_VERSION
static double at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static double maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static int dym_copy_of_svdcmp( double **a, int m, int n, double *w, double **v )
{
        int flag,i,its,j,jj,k,l,nm;
        double c,f,h,s,x,y,z;
        double anorm=0.0,g=0.0,scale=0.0;
        double *rv1 = am_malloc_realnums(n+1);   /* altered */

        if ( Verbosity > 10000.0 )
        {
          printf("Entered dym_copy_of_svdcmp().\n");
          printf("m (number rows) = %d\n",m);
          printf("n (number cols) = %d\n",n);

          fprintf_2d_realnums(stdout,"a",a,m+1,n+1,"\n");
          wait_for_key();
        }

        if (m < n)
          my_error("SVDCMP: You must augment A with extra zero rows");

        for (i=1;i<=n;i++) {
                l=i+1;
                rv1[i]=scale*g;
                g=s=scale=0.0;
                if (i <= m) {
                        for (k=i;k<=m;k++) scale += fabs(a[k][i]);
                        if (scale) {
                                for (k=i;k<=m;k++) {
                                        a[k][i] /= scale;
                                        s += a[k][i]*a[k][i];
                                }
                                f=a[i][i];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                a[i][i]=f-g;
                                if (i != n) {
                                        for (j=l;j<=n;j++) {
                                                for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                                                f=s/h;
                                                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                                        }
                                }
                                for (k=i;k<=m;k++) a[k][i] *= scale;
                        }
                }
                w[i]=scale*g;
                g=s=scale=0.0;
                if (i <= m && i != n) {
                        for (k=l;k<=n;k++) scale += fabs(a[i][k]);
                        if (scale) {
                                for (k=l;k<=n;k++) {
                                        a[i][k] /= scale;
                                        s += a[i][k]*a[i][k];
                                }
                                f=a[i][l];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                a[i][l]=f-g;
                                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                                if (i != m) {
                                        for (j=l;j<=m;j++) {
                                                for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                                                for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                                        }
                                }
                                for (k=l;k<=n;k++) a[i][k] *= scale;
                        }
                }
                anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
        }
        for (i=n;i>=1;i--) {
                if (i < n) {
                        if (g) {
                                for (j=l;j<=n;j++)
                                        v[j][i]=(a[i][j]/a[i][l])/g;
                                for (j=l;j<=n;j++) {
                                        for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                                        for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                                }
                        }
                        for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
                }
                v[i][i]=1.0;
                g=rv1[i];
                l=i;
        }
        for (i=n;i>=1;i--) {
                l=i+1;
                g=w[i];
                if (i < n)
                        for (j=l;j<=n;j++) a[i][j]=0.0;
                if (g) {
                        g=1.0/g;
                        if (i != n) {
                                for (j=l;j<=n;j++) {
                                        for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                                        f=(s/a[i][i])*g;
                                        for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                                }
                        }
                        for (j=i;j<=m;j++) a[j][i] *= g;
                } else {
                        for (j=i;j<=m;j++) a[j][i]=0.0;
                }
                ++a[i][i];
        }
        for (k=n;k>=1;k--) {
                for (its=1;its<=30;its++) {
                        flag=1;
                        for (l=k;l>=1;l--) {
                                nm=l-1;
                                if (fabs(rv1[l])+anorm == anorm) {
                                        flag=0;
                                        break;
                                }
                                if (fabs(w[nm])+anorm == anorm) break;
                        }
                        if (flag) {
                                c=0.0;
                                s=1.0;
                                for (i=l;i<=k;i++) {
                                        f=s*rv1[i];
                                        if (fabs(f)+anorm != anorm) {
                                                g=w[i];
                                                h=PYTHAG(f,g);
                                                w[i]=h;
                                                h=1.0/h;
                                                c=g*h;
                                                s=(-f*h);
                                                for (j=1;j<=m;j++) {
                                                        y=a[j][nm];
                                                        z=a[j][i];
                                                        a[j][nm]=y*c+z*s;
                                                        a[j][i]=z*c-y*s;
                                                }
                                        }
                                }
                        }
                        z=w[k];
                        if (l == k) {
                                if (z < 0.0) {
                                        w[k] = -z;
                                        for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
                                }
                                break;
                        }
                        if (its == 30) {
                          printf("WARNING: No convergence in 30 SVDCMP iterations.  \nUsing mean output.\n");
                          am_free_realnums(rv1,n+1);
                          return 0; /* FAIL */
                        }
                        x=w[l];
                        nm=k-1;
                        y=w[nm];
                        g=rv1[nm];
                        h=rv1[k];
                        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                        g=PYTHAG(f,1.0);
                        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                        c=s=1.0;
                        for (j=l;j<=nm;j++) {
                                i=j+1;
                                g=rv1[i];
                                y=w[i];
                                h=s*g;
                                g=c*g;
                                z=PYTHAG(f,h);
                                rv1[j]=z;
                                c=f/z;
                                s=h/z;
                                f=x*c+g*s;
                                g=g*c-x*s;
                                h=y*s;
                                y=y*c;
                                for (jj=1;jj<=n;jj++) {
                                        x=v[jj][j];
                                        z=v[jj][i];
                                        v[jj][j]=x*c+z*s;
                                        v[jj][i]=z*c-x*s;
                                }
                                z=PYTHAG(f,h);
                                w[j]=z;
                                if (z) {
                                        z=1.0/z;
                                        c=f*z;
                                        s=h*z;
                                }
                                f=(c*g)+(s*y);
                                x=(c*y)-(s*g);
                                for (jj=1;jj<=m;jj++) {
                                        y=a[jj][j];
                                        z=a[jj][i];
                                        a[jj][j]=y*c+z*s;
                                        a[jj][i]=z*c-y*s;
                                }
                        }
                        rv1[l]=0.0;
                        rv1[k]=f;
                        w[k]=x;
                }
        }
        am_free_realnums(rv1,n+1);
        return 1;
}
#endif

#undef SIGN
#undef MAX
#undef PYTHAG

static void dym_copy_of_svbksb( double **u, double w[], double **v, double b[], double x[], int m, int n )
{
        int jj,j,i;
        double s;
        double *tmp = am_malloc_realnums(n+1); /* altered */

        for (j=1;j<=n;j++) {
                s=0.0;
                if (w[j]) {
                        for (i=1;i<=m;i++) s += u[i][j]*b[i];
                        s /= w[j];
                }
                tmp[j]=s;
        }
        for (j=1;j<=n;j++) {
                s=0.0;
                for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
                x[j]=s;
        }
        am_free_realnums(tmp,n+1);
}




     /*************** END OF STUFF COPIED OUT OF ***************/
     /*************** NUMERICAL RECIPES IN `C'. ***************/

/**** Nicer interface to SVD ***/

static bool tdarr_is_legal(double **tdarr,int rows,int cols)
{
  bool result = TRUE;
  int i,j;
  for ( i = 0 ; result && i < rows ; i++ )
    for ( j = 0 ; result && j < cols ; j++ )
      result = (tdarr[i][j] * 0.0) == 0.0;

  return(result);
}

static bool farr_is_legal(double *farr,int size)
{
  bool result = TRUE;
  int i;
  for ( i = 0 ; result && i < size ; i++ )
    result = (farr[i] * 0.0) == 0.0;

  return(result);
}

static bool Svd_self_check = TRUE;
   /* Set to FALSE if you want SVD to NEVER stop and print warnings,
      even if it produceds infs and NaNs and stuff like that.
   */

bool sing_val_decomp( 
    double **tdarr,
    int rows,
    int cols,
    double **u,
    double *w,
    double **v
  )
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
 
     furthermore v with be orthogonal and v will be orthogonal, meaning
     u^T u = Identity, thus, if square, U^-1 = U^T (ditto for V).

     This in turn means that if tdarr is a square matrix then

       tdarr^-1 = v * DW * u^T where

      DW[ii] = 1 / w[i] and DW[ij] = 0 if i != j

     See pages 59-70 of Numerical Recipes in C for further info.




MODIFICATION: Feb 25th 1995, by Andrew Moore.

   When SVD gets called with an a matrix with all very small values
   (e.g. all O(1e-30)), it was returning some NaN values in the u and
   v vectors, as well as a Inf.

   This is due, I suspect, to magnitude errors in intermediate parts of the
   SVD calculation. The final result needn't have over or underflows. To
   fix this, the matrix passed to SVD is normalized so that its maximum
   magnitude component is 1. The resulting w vector is then unnormalized
   by the normalizing factor. This is okay, based on the fact (unproven but
   obviosish) that if the SVD of A is U Diag(W) V^T then the SVD of
   kA is U Diag(kW) V^T where k is a scalar.

*/
{
  int u_rows = rows;
  int u_cols = cols;
  int v_rows = cols;
  int v_cols = cols;
  int w_size = cols;
  bool result;
  double **nra = am_malloc_2d_realnums(rows+1,cols+1);
  double **nrv = am_malloc_2d_realnums(v_rows+1,v_cols+1);
  double *nrw = am_malloc_realnums(w_size+1);
  int i,j;
  double max_mag_tdarr = 0.0; /* See Modification comments Feb 95 above */
  char *problem;

  if ( cols > rows )
    fprintf(stderr,"sing_val_decomp: warning: fewer rows than columns\n");

  for ( i = 0 ; i < rows ; i++ )
    for ( j = 0 ; j < cols ; j++ )
      max_mag_tdarr = real_max(max_mag_tdarr,fabs(tdarr[i][j]));

  if ( max_mag_tdarr <= 0.0 ) max_mag_tdarr = 1.0; 
             /* Hopeless case. The array is all zeros. Let's let numercial
                recipies SVD cope with this */

  for ( i = 0 ; i < rows ; i++ )
    for ( j = 0 ; j < cols ; j++ )
      nra[i+1][j+1] = tdarr[i][j] / max_mag_tdarr;

  if ( Verbosity > 5000.0 )
  {
    printf("SVD: max_mag_tdarr = %g\n",max_mag_tdarr);
    fprintf_2d_realnums(stdout,"svd_tdarr",tdarr,rows,cols,"\n");
    fprintf_2d_realnums(stdout,"svd_nra",nra,rows+1,cols+1,"\n");
    wait_for_key();
  }

  if ( !dym_copy_of_svdcmp(nra,rows,cols,nrw,nrv) )
  {
    /* If SVD fails */
    if ( Svd_self_check && Verbosity > 60.0 )
    {
      fprintf(stderr,"Singular Value Decomposition choked when it was");
      fprintf(stderr,"asked to use the following matrix A:\n");
      fprintf_2d_realnums(stderr,"A",tdarr,rows,cols,"\n");
      wait_for_key();
    }
    set_2d_realnums_constant(u,rows,cols,0.0);
    set_2d_realnums_constant(v,v_rows,v_cols,0.0);
    set_realnums_constant(w,w_size,0.0);
    result = FALSE;
  }
  else
  {
    for ( i = 0 ; i < u_rows ; i++ )
      for ( j = 0 ; j < u_cols ; j++ )
        u[i][j] = nra[i+1][j+1];

    for ( i = 0 ; i < v_rows ; i++ )
      for ( j = 0 ; j < v_cols ; j++ )
        v[i][j] = nrv[i+1][j+1];

    for ( i = 0 ; i < w_size ; i++ )
      w[i] = nrw[i+1] * max_mag_tdarr;
    result = TRUE;
  }

  if ( Verbosity > 5000.0 )
  {
    printf("SVD: result = %s\n",(result) ? "Okay" : "SVD Failed");
    fprintf_2d_realnums(stdout,"svd_nru",nra,rows+1,cols+1,"\n");
    fprintf_2d_realnums(stdout,"svd_u",u,rows,cols,"\n");
    fprintf_2d_realnums(stdout,"svd_nrv",nrv,v_rows+1,v_cols+1,"\n");
    fprintf_2d_realnums(stdout,"svd_v",v,v_rows,v_cols,"\n");
    fprintf_realnums(stdout,"svd_nrw",nrw,cols+1,"\n");
    fprintf_realnums(stdout,"svd_w",w,cols,"\n");
    wait_for_key();
  }

  problem = NULL;

  if ( !tdarr_is_legal(tdarr,rows,cols) )
    problem = "bad-svd-input";
  else if ( !tdarr_is_legal(u,u_rows,u_cols) || 
       !tdarr_is_legal(v,v_rows,v_cols) ||
       !farr_is_legal(w,w_size)
     )
    problem = "bad-svd-output";
  else if ( !result )
    problem = "svd-failed";

  if ( Svd_self_check && problem != NULL )
  {
    int i,j;
    for ( i = 0 ; i < 6 ; i++ )
      for ( j = 0 ; j <= i ; j++ )
        fprintf(stderr,"*%s",(j==i)?"\n":"");
    fprintf(stderr,"******* Singular Value Decomp has a problem.\n");
    for ( i = 0 ; i < 6 ; i++ )
      for ( j = 0 ; j <= (6-i) ; j++ )
        fprintf(stderr,"*%s",(j==(6-i))?"\n":"");

    fprintf(stderr,"\n\n");
    if ( eq_string(problem,"bad-svd-input") )
    {
      fprintf(stderr,"The problem is that SVD was passed a matrix with NaN\n");
      fprintf(stderr,"and/or Inf entries (see tdarr parameter below)\n");
    }
    else if ( eq_string(problem,"bad-svd-output") )
    {
    fprintf(stderr,"The problem is that SVD was produced a matrix with NaN\n");
      fprintf(stderr,"and/or Inf entries (see u,v,w parameter below)\n");
    }
    else if ( eq_string(problem,"svd-failed") )
      fprintf(stderr,"SVD failed to converge on anything.\n");

    fprintf(stderr,"\n\n");
    fprintf_2d_realnums(stdout,"svd_tdarr",tdarr,rows,cols,"\n");
    fprintf_2d_realnums(stderr,"svd_u",u,rows,cols,"\n");
    fprintf_2d_realnums(stderr,"svd_v",v,v_rows,v_cols,"\n");
    fprintf_realnums(stderr,"svd_w",w,cols,"\n");
    wait_for_key();
  }

  free_2d_realnums(nra,rows+1,cols+1);
  free_2d_realnums(nrv,v_rows+1,v_cols+1);
  am_free_realnums(nrw,w_size+1);

  return(result);
}

/****** End of nicer interface etc ***********/

void dym_svd(dym *d,dym *u,dyv *w_diagonal,dym *v)
{
  if ( dym_rows(d) == 1 && dym_cols(d) == 1 )
  {
    double a = dym_ref(d,0,0);
    dym_set(u,0,0,( a < 0.0 ) ? -1.0 : 1.0);
    dyv_set(w_diagonal,0,fabs(a));
    dym_set(v,0,0,1.0);
  }
  else
  {
    bool okay;
    dym *u_temp = mk_dym(d->rows,d->cols); /* In case d and u share memory */

    assert_dym_shape(u,d->rows,d->cols,"dym_svd::u");
    assert_dyv_shape(w_diagonal,d->cols,"dym_svd::w");
    assert_dym_shape(v,d->cols,d->cols,"dym_svd::v");
    
    okay = sing_val_decomp(d->tdarr,d->rows,d->cols,
                           u_temp->tdarr,w_diagonal->farr,v->tdarr);

    copy_dym(u_temp,u);
    free_dym(u_temp);
    if ( !okay )
      fprintf(stderr,"amdm.c :: SVD failed\n");
  }
}

void make_svd_components(dym *d,dym **u,dyv **w_diagonal,dym **v)
{
  *u = mk_dym(d->rows,d->cols);
  *w_diagonal = mk_dyv(d->cols);
  *v = mk_dym(d->cols,d->cols);
  dym_svd(d,*u,*w_diagonal,*v);
}


void dym_svd_backsub(dym *u,dyv *w_diag,dym *v,dyv *b,dyv *x)
/*
   Solves A x = b where A = u * Diag(w_diag) * v^T

   x is the only thing updated.
*/
{
  double **nru = mk_nrecipes_matrix_from_dym(u);
  double *nrw = mk_nrecipes_vector_from_dyv(w_diag);
  double **nrv = mk_nrecipes_matrix_from_dym(v);
  double *nrb = mk_nrecipes_vector_from_dyv(b);
  double *nrx = mk_nrecipes_vector_from_dyv(x);
  assert_dyv_shape(w_diag,u->cols,"dym_svd_backsub::w_diag");
  assert_dym_shape(v,u->cols,u->cols,"dym_svd_backsub::v");
  assert_dyv_shape(b,u->rows,"dym_svd_backsub::b");
  assert_dyv_shape(x,u->cols,"dym_svd_backsub::x");

  if (Verbosity > 2000.0) {
    fprintf_2d_realnums(stdout,"nru",nru,u->rows+1,u->cols+1,"\n");
    fprintf_realnums(stdout,"nrw",nrw,w_diag->size+1,"\n");
    fprintf_2d_realnums(stdout,"nrv",nrv,v->rows+1,v->cols+1,"\n");
    fprintf_realnums(stdout,"nrb",nrb,b->size+1,"\n");
    fprintf_realnums(stdout,"nrx B4",nrx,x->size+1,"\n");
    wait_for_key();
  }
  dym_copy_of_svbksb(nru,nrw,nrv,nrb,nrx,u->rows,u->cols);

  if (Verbosity > 2000.0) {
    fprintf_realnums(stdout,"nrx AFTER",nrx,x->size+1,"\n");
    wait_for_key();
  }
  copy_nrecipes_vector_to_dyv(nrx,x);
  free_nrecipes_vector(nrw,w_diag);
  free_nrecipes_vector(nrb,b);
  free_nrecipes_vector(nrx,x);
  free_nrecipes_matrix(nru,u);
  free_nrecipes_matrix(nrv,v);
}



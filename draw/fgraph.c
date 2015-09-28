/*
   File:        fgraph.c
   Author:      Andrew W. Moore
   Created:     Sep 1st 1997
   Description: New Test Search Algorithm
                

   Copyright (C) 1997, Andrew Moore
*/

#include "fgraph.h"
#include "ongr.h"
#include "drac.h"

fgraph *mk_fgraph(double xlo,double ylo,double xhi,double yhi,
                  double ag_xlo,double ag_ylo,double ag_xhi,double ag_yhi)
{
  fgraph *fg = AM_MALLOC(fgraph);

  fg->x->lo = xlo;
  fg->y->lo = ylo;
  fg->x->hi = xhi;
  fg->y->hi = yhi;
  fg->x->ag_lo = ag_xlo;
  fg->y->ag_lo = ag_ylo;
  fg->x->ag_hi = ag_xhi;
  fg->y->ag_hi = ag_yhi;

  return fg;
}

double f2a(faxis *fx,double val)
{
  double z = (val - fx->lo) / (fx->hi - fx->lo);
  return fx->ag_lo + (fx->ag_hi - fx->ag_lo) * z;
}

double a2f(faxis *fx,double val)
{
  double z = (val - fx->ag_lo) / (fx->ag_hi - fx->ag_lo);
  return fx->lo + (fx->hi - fx->lo) * z;
}

int fg_get_xy(fgraph *fg,double *x,double *y)
{
  double agx,agy;
  int button = ag_get_xy(&agx,&agy);
  *x = a2f(fg->x,agx);
  *y = a2f(fg->y,agy);
  return(button);
}

void fg_dot(fgraph *fg,double x,double y)
{
  ag_dot(f2a(fg->x,x),f2a(fg->y,y));
}

void fg_circle(fgraph *fg,double x,double y,double ag_radius)
{
  ag_circle(f2a(fg->x,x),f2a(fg->y,y),ag_radius);
}

void fg_disc(fgraph *fg,double x,double y,double ag_radius)
{
  ag_disc(f2a(fg->x,x),f2a(fg->y,y),ag_radius);
}

void fg_line(fgraph *fg,double x,double y,double u,double v)
{
  ag_line(f2a(fg->x,x),f2a(fg->y,y),f2a(fg->x,u),f2a(fg->y,v));
}

void fg_box(fgraph *fg,double x,double y,double u,double v)
{
  ag_box(f2a(fg->x,x),f2a(fg->y,y),f2a(fg->x,u),f2a(fg->y,v));
}

void fg_string(fgraph *fg,double x,double y,char *s)
{
  ag_print(f2a(fg->x,x),f2a(fg->y,y),s);
}

void fg_border(fgraph *fg)
{
  int i;
  int size;
  double lo,hi,delta;

  ag_rectangle(fg->x->ag_lo,fg->y->ag_lo,fg->x->ag_hi,fg->y->ag_hi);
  sensible_limits(fg->x->lo,fg->x->hi,&lo,&hi,&delta);

  size = (int) floor(0.5 + (hi - lo)/delta);

  if ( Verbosity > 15.55 )
  {
    printf("lo = %g\n",lo);
    printf("hi = %g\n",hi);
    printf("delta = %g\n",delta);
    printf("size = %d\n",size);
  }

  for ( i = 0 ; i <= size ; i++ )
  {
    double d = lo + delta * i;
    double normal = (d - fg->x->lo) / (fg->x->hi - fg->x->lo);
    if ( normal >= 0.0 && normal <= 1.0 )
    {
      double s =  fg->x->ag_lo * (1 - normal) + fg->x->ag_hi * normal;
      char buff[100];
      ag_line(s,fg->y->ag_lo,s,fg->y->ag_lo - 4.0);
      sprintf(buff,"%g",d);
      ag_print(s,fg->y->ag_lo - 15.0,buff);
    }
  }


  sensible_limits(fg->y->lo,fg->y->hi,&lo,&hi,&delta);
  size = (int) floor(0.5 + (hi - lo)/delta);
  
  if ( Verbosity > 15.55 )
  {
    printf("lo = %g\n",lo);
    printf("hi = %g\n",hi);
    printf("delta = %g\n",delta);
    printf("size = %d\n",size);
  }

  for ( i = 0 ; i <= size ; i++ )
  {
    double d = lo + delta * i;
    double normal = (d - fg->y->lo) / (fg->y->hi - fg->y->lo);
    if ( normal >= 0.0 && normal <= 1.0 )
    {
      double s = fg->y->ag_lo * (1 - normal) + fg->y->ag_hi * normal;
      char buff[100];
      ag_line(fg->x->ag_lo,s,fg->x->ag_lo - 4.0,s);
      sprintf(buff,"%g",d);
      ag_print(fg->x->ag_lo - 40.0,s,buff);
        /* X coordinate of numbers on the y axis is fg->x->ag_lo - 40.0 */
    }
  }
}

void free_fgraph(fgraph *fg)
{
  AM_FREE(fg,fgraph);
}

drac_params *mk_drac_params(fgraph *fg,int grid_size,int num_contours)
{
  drac_params *dp = AM_MALLOC(drac_params);
  dp->num_cols = grid_size;
  dp->num_rows = grid_size;
  dp->num_contours = num_contours;
  sprintf(dp->rfile,"%s","");
  dp->x_low = fg->x->lo;
  dp->y_low = fg->y->lo;
  dp->x_high = fg->x->hi;
  dp->y_high = fg->y->hi;
  sprintf(dp->title,"%s","");
  sprintf(dp->spoints,"%s","");
  return dp;
}

void fg_contours(fgraph *fg,
                     double height_function(char *data,double x,double y),
                     char *data,
                     int grid_size,
                     int num_contours)
{
  drac_params *dp = mk_drac_params(fg,grid_size,num_contours);
  tgrid *tg;
  egrid *eg;
  stage st[1];
  void init_structs_from_drac_params(
    drac_params *par,
    tgrid **tg,
    egrid **eg,
    stage *st,
    double (*h_fn)(char *data,double x,double y),
    char *data
  );
  void draw_contours(stage *st, tgrid *tg, egrid *eg, int n, double *lowh, double *hih);
  void free_tgrid(tgrid *tg);
  void free_egrid(egrid *eg);


  init_structs_from_drac_params(dp,&tg,&eg,st,height_function,data);

  st->bl_window.x = fg->x->ag_lo;
  st->bl_window.y = fg->y->ag_lo;
  st->tr_window.x = fg->x->ag_hi;
  st->tr_window.y = fg->y->ag_hi;

  draw_contours(st, tg, eg, num_contours,NULL,NULL);
  free_tgrid(tg);
  free_egrid(eg);
  AM_FREE(dp,drac_params);
}

/* Is there a point r=(x,y) on the line segment between 
   (xa,ya) (xb,yb) for which (r - pivot).dir == 0?

   If so, return TRUE and put result in *r_x and *r_y
*/
bool fgpivhelp(double xa,double ya,double xb,double yb,
               dyv *pivot,dyv *dir,
               double *r_x,double *r_y)
{
  double xp = dyv_ref(pivot,0);
  double yp = dyv_ref(pivot,1);
  double xd = dyv_ref(dir,0);
  double yd = dyv_ref(dir,1);
  double num = (xp - xa) * xd + (yp - ya) * yd;
  double den = (xb - xa) * xd + (yb - ya) * yd;
  double lambda = (fabs(den) > 1e-6) ? (num / den) :
                  (den < 0.0) ? (num / -1e-6) : (num / 1e-6);
  bool result = lambda >= 0.0 && lambda <= 1.0;

  *r_x = xa + lambda * (xb - xa);
  *r_y = ya + lambda * (yb - ya);

  return(result);
}
 
/* Draws a line indicating the set of points for
   which (x - pivot).dir == 0 */
void fg_pivot_dir(fgraph *fg,dyv *pivot,dyv *dir)
{
  double x1 = -77.7;
  double x2 = -77.7;
  double y1 = -77.7;
  double y2 = -77.7;
  double x,y;
  bool defined1 = FALSE;
  bool defined2 = FALSE;

  if ( fgpivhelp(fg->x->lo,fg->y->lo,fg->x->hi,fg->y->lo,pivot,dir,&x,&y) )
  { if ( defined1 ) { x2 = x; y2 = y; defined2 = TRUE; }
    else { x1 = x; y1 = y; defined1 = TRUE; } 
  }
    
  if ( fgpivhelp(fg->x->hi,fg->y->hi,fg->x->hi,fg->y->lo,pivot,dir,&x,&y) )
  { if ( defined1 ) { x2 = x; y2 = y; defined2 = TRUE; }
    else { x1 = x; y1 = y; defined1 = TRUE; } 
  }

  if ( fgpivhelp(fg->x->hi,fg->y->hi,fg->x->lo,fg->y->hi,pivot,dir,&x,&y) )
  { if ( defined1 ) { x2 = x; y2 = y; defined2 = TRUE; }
    else { x1 = x; y1 = y; defined1 = TRUE; } 
  }
    
  if ( fgpivhelp(fg->x->lo,fg->y->lo,fg->x->lo,fg->y->hi,pivot,dir,&x,&y) )
  { if ( defined1 ) { x2 = x; y2 = y; defined2 = TRUE; }
    else { x1 = x; y1 = y; defined1 = TRUE; } 
  }
    
  if ( defined1 && defined2 )
    fg_line(fg,x1,y1,x2,y2);
}


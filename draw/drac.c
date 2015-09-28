/*
   File:        drac.c
   Authors:     Andrew W. Moore, Mary Soon Lee
   Created:     Fri Jan 21 09:56:04 EST 1994
   Description: DRAwing Contours
   Last update: 23 Jan 97

   Copyright 1997, Schenley Park Research
*/

#include <stdio.h>
#include <math.h>
#include "drac.h"      /* Drawing contours of functions and 2d arrays */

#define RADIUS 10.0
#define NOT_NEIGH 10
#define NUM_CONTOURS 20
#define NUM_COLUMNS 25
#define NUM_ROWS 25
#define DEFAULT_BL_WINDOW_X  0.0
#define DEFAULT_BL_WINDOW_Y  0.0
#define DEFAULT_TR_WINDOW_X  512.0
#define DEFAULT_TR_WINDOW_Y  512.0
#define FRAME_X_BOTTOM 110.0
#define FRAME_Y_BOTTOM 110.0
#define FRAME_X_TOP 30.0
#define FRAME_Y_TOP 30.0
#define DEFAULT_BL_DOMAIN_X  0.0
#define DEFAULT_BL_DOMAIN_Y  0.0
#define DEFAULT_TR_DOMAIN_X  10.0
#define DEFAULT_TR_DOMAIN_Y  10.0


typedef tgpoint *tgpoint_ptr;
typedef tgpoint_ptr *tgpoint_ptr_ptr;

double default_h_fn(char *data,double x, double y)
{
  double result;
  double tem1 = (x - 2.0) * (x - 2.0) + (y - 4.0) * (y - 4.0);
  double tem2 = (2 * x - 8.0) * (x - 8.0) + (y - 6.0) * (y - 6.0);

  if ( data != NULL )
    my_error("default_h_fn expected NULL data");

  if (tem1 < tem2)
    result = tem1;
  else
    result = tem2;
  return(result);
}

tgrid *malloc_tgrid(int cols, int rows)
{
  tgrid *tg = AM_MALLOC(tgrid);  /* see ../amut/amma.h for definition of
                                    AM_MALLOC macro */
  int c,r;

  tg -> num_cols = cols;
  tg -> num_rows = rows;

  tg -> tgpoints = AM_MALLOC_ARRAY(tgpoint_ptr_ptr,rows);

  for ( r = 0 ; r < rows ; r++ )
  {
    tg -> tgpoints[r] = AM_MALLOC_ARRAY(tgpoint_ptr,cols);
    for ( c = 0 ; c < cols ; c++ )
      tg -> tgpoints[r][c] = AM_MALLOC(tgpoint);
  }

  return(tg);
}

void free_tgrid(tgrid *tg)
{
  int cols = tg->num_cols;
  int rows = tg->num_rows;
  int c,r;
  for ( r = 0 ; r < rows ; r++ )
  {
    for ( c = 0 ; c < cols ; c++ )
      am_free((char *)tg -> tgpoints[r][c],sizeof(tgpoint));
    am_free((char *)tg -> tgpoints[r],sizeof(tgpoint_ptr) * cols);
  }

  am_free((char *)tg->tgpoints,sizeof(tgpoint_ptr_ptr) * rows);
  am_free((char *)tg,sizeof(tgrid));
}
  
tgpoint *tg_ref(tgrid *tg, int c, int r)
{
  if ( c < 0 || c >= tg -> num_cols || r < 0 || r >= tg -> num_rows )
/*    my_error("pwicnbap"); */
  {
     printf("ERROR");
     wait_for_key();
     my_error("blob");
   }
  return(tg->tgpoints[r][c]);
}

/* Initializes a stage */
void init_stage(stage *st, drac_params *par)
{
  (st->bl_window).x = DEFAULT_BL_WINDOW_X + FRAME_X_BOTTOM;
  (st->bl_window).y = DEFAULT_BL_WINDOW_Y + FRAME_Y_BOTTOM;
  (st->tr_window).x = DEFAULT_TR_WINDOW_X - FRAME_X_TOP;
  (st->tr_window).y = DEFAULT_TR_WINDOW_Y - FRAME_Y_TOP;
  (st->bl_domain).x = par->x_low;
  (st->bl_domain).y = par->y_low;
  (st->tr_domain).x = par->x_high;
  (st->tr_domain).y = par->y_high;
}

/* Draws in the boundary of the stage */
void frame_stage(stage *st)
{
  double x_min = (st->bl_window).x;
  double y_min = (st->bl_window).y;
  double x_max = (st->tr_window).x;
  double y_max = (st->tr_window).y;
  int i;
  int size;
  double lo,hi,delta;

  ag_rectangle(x_min,y_min,x_max,y_max);

  sensible_limits(st->bl_domain.x,st->tr_domain.x,&lo,&hi,&delta);
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
    double normal = (d - st->bl_domain.x) / (st->tr_domain.x - st->bl_domain.x);
    if ( normal >= 0.0 && normal <= 1.0 )
    {
      double s =  st->bl_window.x * (1 - normal) + st->tr_window.x * normal;
      char buff[100];
      ag_line(s,y_min,s,y_min - 4.0);
      sprintf(buff,"%g",d);
      ag_print(s,y_min - 15.0,buff);
    }
  }


  sensible_limits(st->bl_domain.y,st->tr_domain.y,&lo,&hi,&delta);
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
    double normal = (d - st->bl_domain.y) / (st->tr_domain.y - st->bl_domain.y);
    if ( normal >= 0.0 && normal <= 1.0 )
    {
      double s =  st->bl_window.y * (1 - normal) + st->tr_window.y * normal;
      char buff[100];
      ag_line(x_min,s,x_min - 4.0,s);
      sprintf(buff,"%g",d);
      ag_print(x_min - 40.0,s,buff);
        /* X coordinate of numbers on the y axis is x_min - 40.0 */
    }
  }
}

/* Draws in the boundary of the stage and titles it */
void frame_and_title_stage(stage *st, char *title, char *xlabel, char *ylabel)
{
  frame_stage(st);
  ag_print((st->bl_window).x, (st->tr_window).y + FRAME_Y_TOP - 10.0, title);
  ag_print(st->bl_window.x + 80.0,st->bl_window.y - 40.0,xlabel);
  ag_print_vertical(st->bl_window.x - FRAME_X_BOTTOM + 25.0,
                    st->tr_window.y - 80.0,ylabel);
}

/* Given an x in h_fn coordinates, this returns the corresponding
   x in the graphics window */
double x_tog(stage *st, double x)
{
  double dom_xmin = (st->bl_domain).x;
  double dom_width = (st->tr_domain).x - dom_xmin;
  double win_xmin = (st->bl_window).x;
  double win_width = (st->tr_window).x - win_xmin;
  double result = win_xmin + (((x - dom_xmin) * win_width) / dom_width);
/* This looks unneccessary ...!!!
  if ((x < dom_xmin) || (x > (st->tr_domain).x))
    my_error("x out of bounds");
*/
  return(result);
}

/* Given y in h_fn coordinates, this returns the corresponding
   y in the graphics window */
double y_tog(stage *st, double y)
{
  double dom_ymin = (st->bl_domain).y;
  double dom_height = (st->tr_domain).y - dom_ymin;
  double win_ymin = (st->bl_window).y;
  double win_height = (st->tr_window).y - win_ymin;
  double z= (y - dom_ymin) / dom_height;
  double result = win_ymin + z * win_height;

  if ( z < -0.001 || z > 1.001 )
    my_error("y out of bounds");

  return(result);
}

/* Given an x in graphics coordinates, this returns the corresponding
   x in h_fn coordinates */
double x_fromg(stage *st, double x)
{
  double dom_xmin = (st->bl_domain).x;
  double dom_width = (st->tr_domain).x - dom_xmin;
  double win_xmin = (st->bl_window).x;
  double win_width = (st->tr_window).x - win_xmin;
  double result = dom_xmin + (((x - win_xmin) * dom_width) / win_width);

  if ((x < win_xmin) || (x > (st->tr_window).x))
    my_error("x_g out of bounds");
  return(result);
}

/* Given y in graphics coordinates, this returns the corresponding
   y in h_fn coordinates */
double y_fromg(stage *st, double y)
{
  double dom_ymin = (st->bl_domain).y;
  double dom_height = (st->tr_domain).y - dom_ymin;
  double win_ymin = (st->bl_window).y;
  double win_height = (st->tr_window).y - win_ymin;
  double result = dom_ymin + (((y - win_ymin) * dom_height) / win_height);

  if ((y < win_ymin) || (y > (st->tr_window).y))
    my_error("y_g out of bounds");
  return(result);
}

/* Given a spoint p in h_fn coordinates, this draws it in the right place
   in the graphics window.  */
void drac_dot(stage *st, spoint *p)
{
  ag_dot(x_tog(st,p->x), y_tog(st,p->y));
}

/* Given a spoint p in h_fn coordinates, this highlights the corresponding
   place in the graphics window with a circle.  */
void drac_highlight(stage *st, spoint *p)
{
  ag_circle(x_tog(st,p->x), y_tog(st,p->y), RADIUS);
}

/* Given spoints p1, p2 in h_fn coordinates, this draws a line between the
   corresponding spoints in the graphics window.  */
void drac_line(stage *st, spoint *p1, spoint *p2)
{
  ag_line(x_tog(st,p1->x), y_tog(st,p1->y), 
          x_tog(st,p2->x), y_tog(st,p2->y));
}

/* Draws the spoint represented by tp */
void draw_tgpoint(stage *st, tgpoint *tp)
{
  drac_dot(st, &(tp->spnt));
}

/* Highlights the spoint represented by tp with a circle */
void highlight_tgpoint(stage *st, tgpoint *tp)
{
  drac_highlight(st, &(tp->spnt));
}

/* Initializes the tg_i, tg_j, index, spoint, and height fields in the 
   tg_spoints of a tg_grid.  */
void init_tg_basics(
    drac_params *par, 
    tgrid *tg, 
    double (*h_fn)(char *data,double x,double y),
    char *data
  )
{
  int cols = tg->num_cols;
  int rows = tg->num_rows;
  int i, j;
  double tpx, tpy;
  tgpoint *tp;

  if ((cols <= 2) || (rows <= 2))
    my_error("Too few columns or rows for a tg_grid");
  for (i = 0; i < cols; i++)
    for (j = 0; j < rows; j++)
    {
      tp = tg_ref(tg, i, j);
      tp->tg_i = i;
      tp->tg_j = j;
      tp->index = i + j * cols;
      tpx = par->x_low + 
        (i / ((double) (cols - 1))) * (par->x_high - par->x_low);
      tpy = par->y_low + 
        (j / ((double) (rows - 1))) * (par->y_high - par->y_low);
      (tp->spnt).x = tpx;
      (tp->spnt).y = tpy;
      tp->height = h_fn(data,tpx, tpy);
    }
}

/* Initializes the neighbor fields in the tg_spoints of a tg_grid. 
   This is an ugly function, split into lots of special cases.  */
void init_tg_neighs(tgrid *tg)
{
  int cols = tg->num_cols;
  int rows = tg->num_rows;
  int i, j;
  tgpoint *tp;
  tgpoint **neighs;

  /* first set loop to be FALSE for all edge-spoints, TRUE otherwise */
  for (i = 0; i < cols; i++)
    for (j = 0; j < rows; j++)
    {
      tp = tg_ref(tg, i, j);
      if ((i == 0) || (j == 0) || (i == (cols - 1)) || (j == (rows - 1)))
        tp->loop = FALSE;
      else
        tp->loop = TRUE;
    }

  /* the 4 corner spoints */
  tp = tg_ref(tg, 0, 0);
  tp->num_neighs = 2;
  neighs = tp->neighs;
  neighs[0] = tg_ref(tg, 0, 1);
  neighs[1] = tg_ref(tg, 1, 0);

  tp = tg_ref(tg, 0, rows - 1);
  tp->num_neighs = 3;
  neighs = tp->neighs;
  neighs[0] = tg_ref(tg, 1, rows - 1);
  neighs[1] = tg_ref(tg, 1, rows - 2);
  neighs[2] = tg_ref(tg, 0, rows - 2);

  tp = tg_ref(tg, cols - 1, rows - 1);
  tp->num_neighs = 2;
  neighs = tp->neighs;
  neighs[0] = tg_ref(tg, cols - 1, rows - 2);
  neighs[1] = tg_ref(tg, cols - 2, rows - 1);
  
  tp = tg_ref(tg, cols - 1, 0);
  tp->num_neighs = 3;
  neighs = tp->neighs;
  neighs[0] = tg_ref(tg, cols - 2, 0);
  neighs[1] = tg_ref(tg, cols - 2, 1);
  neighs[2] = tg_ref(tg, cols - 1, 1);

  /* the center spoints */
  for (i = 1; i < cols - 1; i++)
    for (j = 1; j < rows - 1; j++)
    {
      tp = tg_ref(tg, i, j);
      tp->num_neighs = 6;
      neighs = tp->neighs;
      neighs[0] = tg_ref(tg, i, j + 1);
      neighs[1] = tg_ref(tg, i + 1, j);
      neighs[2] = tg_ref(tg, i + 1, j - 1);
      neighs[3] = tg_ref(tg, i, j - 1);
      neighs[4] = tg_ref(tg, i - 1, j);
      neighs[5] = tg_ref(tg, i - 1, j + 1);
    }

  /* the 4 edges */
  for (j = 1; j < rows - 1; j++)    /* left edge */
  {
    tp = tg_ref(tg, 0, j);
    tp->num_neighs = 4;
    neighs = tp->neighs;
    neighs[0] = tg_ref(tg, 0, j + 1);
    neighs[1] = tg_ref(tg, 1, j);
    neighs[2] = tg_ref(tg, 1, j - 1);
    neighs[3] = tg_ref(tg, 0, j - 1);
  }
  for (i = 1; i < cols - 1; i++)    /* top edge */
  {
    tp = tg_ref(tg, i, rows - 1);
    tp->num_neighs = 4;
    neighs = tp->neighs;
    neighs[0] = tg_ref(tg, i + 1, rows - 1);
    neighs[1] = tg_ref(tg, i + 1, rows - 2);
    neighs[2] = tg_ref(tg, i, rows - 2);
    neighs[3] = tg_ref(tg, i - 1, rows - 1);
  }
  for (j = 1; j < rows - 1; j++)    /* right edge */
  {
    tp = tg_ref(tg, cols - 1, j);
    tp->num_neighs = 4;
    neighs = tp->neighs;
    neighs[0] = tg_ref(tg, cols - 1, j - 1);
    neighs[1] = tg_ref(tg, cols - 2, j);
    neighs[2] = tg_ref(tg, cols - 2, j + 1);
    neighs[3] = tg_ref(tg, cols - 1, j + 1);
  }
  for (i = 1; i < cols - 1; i++)    /* bottom edge */
  {
    tp = tg_ref(tg, i, 0);
    tp->num_neighs = 4;
    neighs = tp->neighs;
    neighs[0] = tg_ref(tg, i - 1, 0);
    neighs[1] = tg_ref(tg, i - 1, 1);
    neighs[2] = tg_ref(tg, i, 1);
    neighs[3] = tg_ref(tg, i + 1, 0);
  }
}

/* Initializes the tg_spoints of a tg_grid.  */
void init_tg(
    drac_params *par,
    tgrid *tg,
    double (*h_fn)(char *data,double x,double y),
    char *data
  )
{
  init_tg_basics(par, tg, h_fn, data);
  init_tg_neighs(tg);
}

/* Draws in all the spoints on a tgrid.  */
void draw_tgrid(stage *st, tgrid *tg)
{
  int cols = tg->num_cols;
  int rows = tg->num_rows;
  int i, j;

  for (i = 0; i < cols; i++)
    for (j = 0; j < rows; j++)
    {
      draw_tgpoint(st, tg_ref(tg, i, j));
    }
}

/* Returns the nearest tgpoint to p (p in h_fn coordinates). 
   "nearest" means the tgpoint whose x coordinate is closest to p's
   x coordinate, and similarly for y.  */
tgpoint *nearest_tgpoint(stage *st, tgrid *tg, spoint *p)
{
  int cols = tg->num_cols;  
  int rows = tg->num_rows;
  double dom_xmin = st->bl_domain.x;
  double dom_ymin = st->bl_domain.y;
  double dom_xmax = st->tr_domain.x;
  double dom_ymax = st->tr_domain.y;
  double dom_width = dom_xmax - dom_xmin;
  double dom_height = dom_ymax - dom_ymin;
  int i_l; /* the tg_i of the column of tgpoints just to the left of p */
  int j_b; /* the tg_j of the row of tgpoints just below p */
  tgpoint *p1, *p2;
  int i_needed, j_needed;

  if ((p->x <= dom_xmin) || (p->x >= dom_xmax) ||
      (p->y <= dom_ymin) || (p->y >= dom_ymax))
    my_error("p is outside the region being considered");
  i_l = (int) (floor(((p->x - dom_xmin) * ((double) (cols - 1))) / dom_width));
  j_b = (int) (floor(((p->y - dom_ymin) * ((double)( rows - 1))) / dom_height));

  /* the i we need is either i_l or i_l + 1, similarly for j.... */
  p1 = tg_ref(tg, i_l, j_b);
  p2 = tg_ref(tg, i_l + 1, j_b + 1);
  if ((p->x - p1->spnt.x) < (p2->spnt.x - p->x))
    i_needed = i_l;
  else 
    i_needed = i_l + 1;
  if ((p->y - p1->spnt.y) < (p2->spnt.y - p->y))
    j_needed = j_b;
  else 
    j_needed = j_b + 1;

  return(tg_ref(tg, i_needed, j_needed));
}

/* Waits for the mouse to be clicked within the stage; then returns
   the click's location IN H_FN COORDINATES in p.  */
void drac_use_mouse(stage *st, spoint *p)
{
  double x, y;
  bool not_yet_found = TRUE;
  double win_xmin = st->bl_window.x;
  double win_ymin = st->bl_window.y;
  double win_xmax = st->tr_window.x;
  double win_ymax = st->tr_window.y;

  while (not_yet_found)
  {
    ag_get_xy(&x, &y);
    if ((x <= win_xmin) || (x >= win_xmax) || 
        (y <= win_ymin) || (y >= win_ymax))
      printf("Clicked outside stage....  Try again.\n");
    else
      not_yet_found = FALSE;
  }
  p->x = x_fromg(st, x);
  p->y = y_fromg(st, y);
}

/* Draws lines connecting all the neighbors of tgp.  */
void display_loop(stage *st, tgpoint *tgp)
{
  bool is_loop = tgp->loop;
  int n;
  int max_n = tgp->num_neighs;
  tgpoint **neighs = tgp->neighs;
  tgpoint *tgp1, *tgp2;

  if (is_loop)
  {
    for (n = 0; n < max_n; n++)
    {
      tgp1 = neighs[n];
      if (n == (max_n - 1))
        tgp2 = neighs[0];
      else
        tgp2 = neighs[n + 1];
      drac_line(st, &(tgp1->spnt), &(tgp2->spnt));
    }
  }
  else
  {
    for (n = 0; n < (max_n - 1); n++)
    {
      tgp1 = neighs[n];
      tgp2 = neighs[n + 1];
      drac_line(st, &(tgp1->spnt), &(tgp2->spnt));
    }
  }
}

/* If tgp1 and tgp2 are neighbors, this returns i s.t. 
     tgp1->neighs[i] = tgp2
   else it returns NOT_NEIGH */
int neigh_num(tgpoint *tgp1, tgpoint *tgp2)
{
  int i;
  int result = NOT_NEIGH;
  for (i = 0; i < tgp1->num_neighs; i++)
  {
    if (tgp1->neighs[i] == tgp2)
      result = i;
  }
  return(result);
}

/* This draws an egde */
void draw_tg_edge(stage *st, edge *e)
{
  drac_line(st, &((e->start)->spnt), &((e->end)->spnt));
}

typedef edge *edge_ptr;
typedef edge_ptr *edge_ptr_ptr;

/* This creates the egrid given the number of columns and rows in
   the triangulation grid.  */
egrid *malloc_egrid(int cols, int rows)
{
  egrid *eg = AM_MALLOC(egrid);  /* see ../amut/amma.h for definition of
                                    AM_MALLOC macro */
  int i, tgpoint_index;
  int num_tgpoints = cols * rows;

  eg -> num_cols = cols;
  eg -> num_rows = rows;

  eg -> edges = AM_MALLOC_ARRAY(edge_ptr_ptr, num_tgpoints);

  for (tgpoint_index = 0; tgpoint_index < num_tgpoints; tgpoint_index++)
  {
    eg -> edges[tgpoint_index] = AM_MALLOC_ARRAY(edge_ptr, MAX_NEIGHS);
    for (i = 0; i < MAX_NEIGHS; i++)
      eg -> edges[tgpoint_index][i] = AM_MALLOC(edge);
  }
  return(eg);
}

void free_egrid(egrid *eg)
{
  int rows = eg -> num_rows;
  int cols = eg -> num_cols;
  int i, tgpoint_index;
  int num_tgpoints = cols * rows;

  for (tgpoint_index = 0; tgpoint_index < num_tgpoints; tgpoint_index++)
  {
    for (i = 0; i < MAX_NEIGHS; i++)
      am_free((char *)eg -> edges[tgpoint_index][i],sizeof(edge));
    am_free((char *)eg -> edges[tgpoint_index],MAX_NEIGHS * sizeof(edge_ptr));
  }

  am_free((char *)eg->edges,num_tgpoints * sizeof(edge_ptr_ptr));
  am_free((char *)eg,sizeof(egrid));
}

/* Given two neighboring tgpoints, this returns the edge between them */
edge *e_ref(egrid *eg, tgpoint *tgp1, tgpoint *tgp2)
{
  int i = neigh_num(tgp1, tgp2);

  if (i == NOT_NEIGH)
    my_error("Trying to find edge between non-neighboring spoints");
  return(eg->edges[tgp1->index][i]);
}

/* Given a tgpoint and the neigh_num of its neighbors this returns the edge
   between them with no checks */
edge *fast_e_ref(egrid *eg, tgpoint *tgp1, int i)
{
  return(eg->edges[tgp1->index][i]);
}

/* Given the index for a tgpoint and the neigh_num of its neighbors this
   returns the edge between them with no checks */
edge *vfast_e_ref(egrid *eg, int index, int i)
{
  return(eg->edges[index][i]);
}

/* Initializes edge to run between tgp1 and tgp2, and fills in the
   is_edge and visited fields */
void init_edge_basic(edge *e, tgpoint *tgp1, tgpoint *tgp2)
{
  e->start = tgp1;
  e->end = tgp2;
  e->visited = FALSE;
  e->is_edge = TRUE;
}

/* Initializes the start, end, visited and is_edge fields of the edges in
   an egrid */
void init_eg_basics(tgrid *tg, egrid *eg)
{
  int cols = tg->num_cols;
  int rows = tg->num_rows;
  int i, j, k;
  tgpoint *tgp;
  edge *e;

  for (j = 0; j < rows ; j++)
    for (i = 0; i < cols; i++)
    {
      tgp = tg_ref(tg, i, j);
      for (k = 0; k < tgp->num_neighs; k++)
      {
        e = fast_e_ref(eg, tgp, k);
        init_edge_basic(e, tgp, tgp->neighs[k]);
      }
      for (k = tgp->num_neighs; k < MAX_NEIGHS; k++)
      {
        e = fast_e_ref(eg, tgp, k);
        e->is_edge = FALSE;
      }
    }
}

/* 17 Oct 96: Mary deleted unreferenced formal parameter "stage *st" */
/* Initializes the remaining fields in the edges of an egrid, i.e.
   border_edge, tr1_e1, tr1_e2, tr2_e1, tr2_e2 */
void init_eg_rest(tgrid *tg, egrid *eg)
{
  int num_tgpoints = tg->num_cols * tg->num_rows;
  int i, j, k, l;
  edge *e;
  tgpoint *tgp1, *tgp2;
  bool found_neigh1, found_neigh2;

  for (i = 0; i < num_tgpoints; i++)
    for (j = 0; j < MAX_NEIGHS; j++)
    {
      e = vfast_e_ref(eg, i, j);
      if (e->is_edge == TRUE)
      {
        tgp1 = e->start;
        tgp2 = e->end;
        e->border_edge = TRUE; /* if e is not on the border, this should
                                  later get set to FALSE */
        /* Find the common neighbors of tgp1 and tgp2; there should be
           at least one and at most two.  Once a neighbor is found, fill
           in the appropriate fields in e */
        found_neigh1 = FALSE;
        found_neigh2 = FALSE;
        for (k = 0; k < tgp1->num_neighs; k++)
          for (l = 0; l < tgp2->num_neighs; l++)
          {
            if (tgp1->neighs[k] == tgp2->neighs[l])
            {
              if (found_neigh1 == FALSE)
              {
                found_neigh1 = TRUE;
                e->tr1_e1 = e_ref(eg, tgp1, tgp1->neighs[k]);
                e->tr1_e2 = e_ref(eg, tgp2, tgp1->neighs[k]);
              }
              else if (found_neigh2 == FALSE)
              {  
                found_neigh2 = TRUE;
                e->border_edge = FALSE;
                e->tr2_e1 = e_ref(eg, tgp1, tgp1->neighs[k]);
                e->tr2_e2 = e_ref(eg, tgp2, tgp1->neighs[k]);
              }
              else 
                my_error("Edge spoints should not have three common neighbors");
            }
          }
        if (found_neigh1 == FALSE)
          my_error("Failed to find a common neighbor during init_eg_rest");
      }
    }
}

/* 17 Oct 96: Mary deleted unreferenced formal parameter "stage *st" */
/* Fully initializes an egrid */
void init_eg(tgrid *tg, egrid *eg)
{
  init_eg_basics(tg, eg);
  init_eg_rest(tg, eg);
}

/*  Draws the triangles in the triangulation grid that have e as an
    edge */
void draw_triangles(stage *st, edge *e)
{
  draw_tg_edge(st, e->tr1_e1);
  draw_tg_edge(st, e->tr1_e2);
  if (e->border_edge == FALSE)
  {
    draw_tg_edge(st, e->tr2_e1);
    draw_tg_edge(st, e->tr2_e2);
  }
}

/* Returns TRUE iff edge e crosses height h, i.e. one end-spoint
   is strictly higher than h, the other less than or equal to h.  */
bool crosses_h(edge *e, double h)
{
  double h1 = (e->start)->height;
  double h2 = (e->end)->height;
  if (((h1 <= h) && (h2 > h)) || ((h1 > h) && (h2 <= h)))
    return(TRUE);
  else
    return(FALSE);
}

/* Given an edge, e, whose end tgpoints are on opposite sides of height h,
   this sets cross_pt to be the (x, y) spoint in h_fn coordinates
   that lies on the edge, e, and that linear interpolation gives as 
   having height h */
void set_cross_pt(edge *e, double h, spoint *cross_pt)
{
  tgpoint *tgp1 = e->start;
  tgpoint *tgp2 = e->end;
  double ratio = (tgp1->height - h) / (tgp1->height - tgp2->height);
           /* the ratio d1/d2, where d1 is the distance from tgp1
              to the estimated crossing spoint of the h contour on e,
              and d2 is the distance from tgp2 to tgp1 */
  cross_pt->x = (tgp1->spnt).x + ((tgp2->spnt).x - (tgp1->spnt).x) * ratio;
  cross_pt->y = (tgp1->spnt).y + ((tgp2->spnt).y - (tgp1->spnt).y) * ratio;
}

/* Given an edge, e, that crosses height h, this first performs linear
   interpolation to estimate where the h contour crosses the edge.
   Then for each triangle that e is in (at most two), it finds the
   edge that also crosses the h contour, draws in the estimate of the
   h contour across the triangle, and then iterates to continue tracing
   the contour around.  The algorithm stops trying to extend the
   contour when all the relevant edges have either already been visited
   (as will eventually happen in a loop) or lead off the boundary
   of the grid.  It updates p to be the (x, y) spoint in h_fn coordinates
   where the contour crosses e. */
void draw_contour(stage *st, edge *e, double h, spoint *cross_pt)
{
  edge *e1, *e2;
  spoint p;
  int i;
  int num_triangles_to_try;

  e->visited = TRUE;
  set_cross_pt(e, h, cross_pt);
  if (e->border_edge == TRUE)
     num_triangles_to_try = 1;
  else
      num_triangles_to_try = 2;
  for (i = 0; i < num_triangles_to_try; i++)
  {
    if (i == 0)
    { /* will try first triangle */
      e1 = e->tr1_e1;
      e2 = e->tr1_e2;
    }
    else 
    { /* will try second triangle */
      e1 = e->tr2_e1;
      e2 = e->tr2_e2;
    }
    if (crosses_h(e1, h) == TRUE)
    {
      if (e1->visited == FALSE)
        draw_contour(st, e1, h, &p);
      else  
        /* even if the edge has been visited, we need to find its
           crossing spoint, p, and draw the line segment across this triangle,
           otherwise loops would not be completed */
        set_cross_pt(e1, h, &p);
      drac_line(st, cross_pt, &p);
    }
    else if (crosses_h(e2, h) == TRUE)
    {
      if (e2->visited == FALSE)
        draw_contour(st, e2, h, &p);
      else  
        set_cross_pt(e2, h, &p);
      drac_line(st, cross_pt, &p);
    }
    else
      my_error("h contour should cross one of e1, e2");
  }
}  

/* This marks all the edges in the grid as unvisited */
void clear_edges(egrid *eg)
{
  int num_tgpoints = eg->num_cols * eg->num_rows;
  int i, j;
  edge *e;

  for (i = 0; i < num_tgpoints; i++)
    for (j = 0; j < MAX_NEIGHS; j++)
    {
      e = vfast_e_ref(eg, i, j);
      e->visited = FALSE;
    }
}

/* This draws in the contours at height h.  We progress through the
   edges in the egrid.  Whenever we find an unvisited edge whose
   endspoints are opposite sides of h (so that height h lies in between),
   we look at its neighboring edges and trace the contour round.
   Linear interpolation is used to decide where to draw the contour.  */
void draw_h_contours(stage *st, egrid *eg, double h)
{
  int num_tgpoints = eg->num_cols * eg->num_rows;
  int i, j;
  edge *e;
  spoint p;

  clear_edges(eg);
  for (i = 0; i < num_tgpoints; i++)
    for (j = 0; j < MAX_NEIGHS; j++)
    {
      e = vfast_e_ref(eg, i, j);
      if ((e->is_edge == TRUE) && (e->visited == FALSE))
      {
        if (crosses_h(e, h) == TRUE)
          draw_contour(st, e, h, &p);
      }
    }
}

double tgrid_min(tgrid *tg)
{
  int i,j;
  double result = 0.0;
  for (i = 0; i < tg->num_cols; i++)
    for(j = 0; j < tg->num_rows; j++)
    {
      double cur_h = tg_ref(tg, i, j)->height;
      if ( (i==0 && j==0) || cur_h < result )
        result = cur_h;
    }
  return(result);
}        

double tgrid_max(tgrid *tg)
{
  int i,j;
  double result = 0.0;
  for (i = 0; i < tg->num_cols; i++)
    for(j = 0; j < tg->num_rows; j++)
    {
      double cur_h = tg_ref(tg, i, j)->height;
      if ( (i==0 && j==0) || cur_h > result )
        result = cur_h;
    }
  return(result);
}        

/* Autonomously draws in n contours for the egrid.  It starts by
   finding the lowest and highest heights on the tgrid--low_h and
   high_h, then draws in contours at 
      low_h + d, low_h + 2d, ... , low_h + n * d
   where d = (high_h - low_h)/(n + 1).
       N.B. There is little spoint in wasting effort attempting to
   draw contours at low_h and high_h.
*/
void old_draw_contours(stage *st, tgrid *tg, egrid *eg, int n)
{
  int num_cols = tg->num_cols;
  int num_rows = tg->num_rows;
  double low_height;
  double high_height;
  double d;
  int i;
  char height_info[100];
  int old_pen_color = ag_pen_color();

  /* first some checks */
  if (n < 1) 
    my_error("Asked to draw too few contours");
  if ((num_cols != eg->num_cols) || (num_rows != eg->num_rows))
    my_error("Egrid and tgrid have different dimensions");

  /* now find high and low values */
  low_height = tgrid_min(tg);
  high_height = tgrid_max(tg);

  /* now print the height information and draw the contours */
  d = (high_height - low_height) / (n + 1);
  sprintf(height_info, "Contours from %.3g to %.3g in increments of %.3g",
          low_height + d, low_height + n * d, d);
  ag_print((st->bl_window).x, (st->bl_window).y - FRAME_Y_BOTTOM + 10.0,
           height_info
          );
  for (i = 1; i <= n; i++)
  {
    ag_set_pen_color(ag_spectrum_color(
                     (((double) i) - 1.0) / 
                           (((double) n) - 1.0)));
    draw_h_contours(st, eg, low_height + i * d);
  }
  ag_set_pen_color(old_pen_color);
}

/* Given the minimum and maximum known heights, min_h and max_h, and
   the approximate number n of contours to draw, this updates min_con and 
   max_con to hold the lowest and highest contours.  hdiff is updated to
   hold the difference in height between each contour.
        The goal is to return contours at relatively round heights,
   i.e. 14 rather than 13.82, etc. */
void contour_limits(double min_h, double max_h, int n, 
                    double *min_con, double *max_con, double *hdiff )
{
  double scale,rel_hi,rel_delta, diff;
  int num;

  if ( max_h - min_h < 1e-5 )
  {
    double xmid = (min_h + max_h)/2.0;
    min_h = xmid - 1e-5;
    max_h = xmid + 1e-5;
  }

  scale = pow(10.0,ceil(-log10(max_h - min_h)));
  rel_hi = scale * (max_h - min_h);
  /* rel_hi is now a number in the range [1.0, 10.0) such that 
          (max_h - min_h) = rel_hi / scale 
     and scale is an exact power of 10. */

  rel_delta = ( rel_hi < 2.0 ) ? 0.1 :
              ( rel_hi < 4.0 ) ? 0.2 :
              ( rel_hi < 6.0 ) ? 0.25 :
                0.5;

  diff = rel_delta / scale;
  num = (int) floor((max_h - min_h) / diff);  /* The number of contours
                                            that would result if we drew
                                            them at intervals of hdiff.
                                            Because of how we chose
                                            rel_delta this is in the
                                            range [10, 24]. */
  if (num * 6 > n * 10)
  {
    /* I am prepared to reduce the number of contours to the range
       [2, 5] but no lower. */
    *hdiff = 5.0 * diff;
  }
  else if (num * 8 < n * 5)
  {
    if (rel_delta == 0.25)
      rel_delta = 0.1;
    else 
      rel_delta = rel_delta / 2.0;
    *hdiff = rel_delta / scale;
    /* the number of contours will now be somewhere in the range [20, 60].
       I am not prepared to raise it further.  */
  }
  else
    *hdiff = rel_delta / scale;
  *min_con = *hdiff * ceil(min_h / *hdiff);
  *max_con = *hdiff * floor(max_h / *hdiff);
}

/* Autonomously draws approximately n contours for the egrid.  The
   contours are selected to be at evenly spaced but relatively round
   heights, e.g. 118, 120, ...., 132 rather than old_draw_contours
   which might have produced contours at 117.62, 119.79, etc.

        It starts by finding the lowest and highest heights on the
   tgrid (low_h and high_h) and then uses sensible_limits to decide
   the values of the lowest and highest contours to draw.
        Modified to read its arguments instead if available JGS 11-11-96

        There are some itty bitty details to do with things like
    printing a color key to some of the heights.
*/
void draw_contours(stage *st, tgrid *tg, egrid *eg, int n, double *lowh, double *hih)
{
  int num_cols = tg->num_cols;
  int num_rows = tg->num_rows;
  double low_height;
  double high_height;
  double d, min_con, max_con, cur_con;
  char height_info[100];
  int old_pen_color = ag_pen_color();
  int num_hts_printed = 0;
  int num_to_mention = 6;    /* Num heights to show in color key */

  /* first some checks */
  if (n < 1) 
    my_error("Asked to draw too few contours");
  if ((num_cols != eg->num_cols) || (num_rows != eg->num_rows))
    my_error("Egrid and tgrid have different dimensions");

  /* now find high and low values */
  if (lowh) low_height = *lowh;
  else      low_height = tgrid_min(tg);
  if (hih)  high_height = *hih;
  else      high_height = tgrid_max(tg);

  contour_limits(low_height, high_height, n, &min_con, &max_con, &d);
  sprintf(height_info, "Contours from %.3g to %.3g in increments of %.3g",
          min_con, max_con, d);
  ag_print((st->bl_window).x + 10.0, 
           (st->bl_window).y - FRAME_Y_BOTTOM + 28.0,
           height_info
          );

  ag_set_pen_color(AG_HIGHLIGHT);
  ag_box((st->bl_window).x,
         (st->bl_window).y - FRAME_Y_BOTTOM,
         (st->tr_window).x,
         (st->bl_window).y - FRAME_Y_BOTTOM + 25.0);
  ag_set_pen_color(AG_BKGD);
  ag_box((st->bl_window).x + 3.0,
         (st->bl_window).y - FRAME_Y_BOTTOM + 3.0,
         (st->tr_window).x - 3.0,
         (st->bl_window).y - FRAME_Y_BOTTOM + 22.0);
  ag_set_pen_color(AG_DISPLAY);
  sprintf(height_info, "Key:");
  ag_print((st->bl_window).x + 10.0,
               (st->bl_window).y - FRAME_Y_BOTTOM + 4.0,
               height_info);

  cur_con = min_con;
  while (cur_con <= max_con)
  {
    ag_set_pen_color(ag_spectrum_color((cur_con - min_con) / 
                                       (max_con - min_con)));
    if (((cur_con - min_con) / (max_con - min_con))
        >= ((float) num_hts_printed) / ((float) num_to_mention))
    {
      sprintf(height_info, "%.3g", 
              ((cur_con < 1e-6) && (cur_con > -1e-6)) ? 0.0 : cur_con);
      ag_print((st->bl_window).x + 
               (((st->tr_window).x - (st->bl_window).x) * 
                (((float) num_hts_printed) + 1.0) / 
                ((float) num_to_mention + 2.0)),
               (st->bl_window).y - FRAME_Y_BOTTOM + 4.0,
               height_info);
      num_hts_printed += 1;
    }
    draw_h_contours(st, eg, cur_con);
    cur_con += d;
  }
  ag_set_pen_color(old_pen_color);
}

/* Reads the parameters from the command line arguments and puts them
   into par. */
void drac_params_from_args(drac_params *par, int argc, char **argv)
{
  par->num_cols = int_from_args("-cols",argc,argv,NUM_COLUMNS);
  par->num_rows = int_from_args("-rows",argc,argv,NUM_ROWS);
  par->num_contours = int_from_args("-nc",argc,argv, NUM_CONTOURS);
  sprintf(par->rfile, string_from_args("-res",argc, argv,"contour.ps"));
  par->x_low = double_from_args("-xlow",argc,argv,DEFAULT_BL_DOMAIN_X);
  par->y_low = double_from_args("-ylow",argc,argv,DEFAULT_BL_DOMAIN_Y);
  par->x_high = double_from_args("-xhigh",argc,argv,DEFAULT_TR_DOMAIN_X);
  par->y_high = double_from_args("-yhigh",argc,argv,DEFAULT_TR_DOMAIN_Y);
  sprintf(par->title, string_from_args("-title",argc, argv,"CONTOURS"));
  sprintf(par->spoints, string_from_args("-spoints",argc, argv,""));
}

/* Given the name of a file, f, this scans through what should be sets
   of (x, y, z) triples to determine the number of rows and columns in
   the grid.  The spoints must occur in the following order: scanning
   along rows from the minimum x to the maximum x, starting with the
   bottom row (minimum y) and rising to the top row.  The spoints must
   lie on a rectangular grid, but there is no check for this. */
void size_of_grid_in_file(char *f, int *cols, int *rows)
{
  FILE *s = fopen(f, "r");
  double cur_x, cur_y, cur_z;
  double min_x = 0.0; /* Just to prevent unitialized warning */
  int i = 0;
  int j = 0;
  int n = 0; /* number of triples scanned */
  int q;
  bool NOTHING_READ = TRUE;
  bool STOP = FALSE;
  if (s == NULL)
    my_error("Failed when trying to open a file to read in a grid");
  while (STOP == FALSE)
  {
    q = fscanf(s, "%lf %lf %lf", &cur_x, &cur_y, &cur_z);
    if (q == EOF) /* test for end of file */
      STOP = TRUE;
    else if (q == 3) /* found a new triple */
    {    
      n++;
      if (NOTHING_READ == TRUE) /* this was the 1st triple; set min_x */
      {
        min_x = cur_x;
        NOTHING_READ = FALSE;
        i = 1;
      }
      else if (min_x == cur_x) /* at the start of a new row */
      {
        i = 1;
        j++;
      }
      else
        i++;
    }
    else
      my_error("Error in format of grid file");
  }      
  if (n != i * (j + 1))
    my_error("Wrong number of spoints in grid file");
  *cols = i;
  *rows = j + 1;
  fclose(s);
}

/* Initializes a tgrid from a file.  The dimensions of the tgrid
   (i.e. num_cols, num_rows) should already be correctly set. */
void init_tg_from_file(tgrid *tg, char *f, drac_params *par)
{
  int cols = tg->num_cols;
  int rows = tg->num_rows;
  int i, j;
  tgpoint *tp;
  FILE *s = fopen(f, "r");
  double cur_x, cur_y, cur_z;

  if (s == NULL)
    my_error("Failed when trying to re-open a file to read in a grid");
  if ((cols <= 2) || (rows <= 2))
    my_error("Too few columns or rows for a tg_grid");

  for (j = 0; j < rows; j++)
    for (i = 0; i < cols; i++)
    {
      fscanf(s, "%lf %lf %lf", &cur_x, &cur_y, &cur_z);
      tp = tg_ref(tg, i, j);
      tp->tg_i = i;
      tp->tg_j = j;
      tp->index = i + j * cols;
      (tp->spnt).x = cur_x;
      (tp->spnt).y = cur_y;
      tp->height = cur_z;
    }
  init_tg_neighs(tg);
  /* Now set the bl_domain and tr_domain spoints correctly in par */
  par->x_low = tg_ref(tg, 0, 0)->spnt.x;
  par->y_low = tg_ref(tg, 0, 0)->spnt.y;
  par->x_high = tg_ref(tg, cols - 1, 0)->spnt.x;
  par->y_high = tg_ref(tg, 0, rows - 1)->spnt.y;
}

/* Initializes the structures from the parameters.  Note that if
   the spoints are being read in from a grid in a file, then parameters
   such as num_cols are deduced from that file, not taken from the
   command line.  */
void init_structs_from_drac_params(
    drac_params *par,
    tgrid **tg,
    egrid **eg,
    stage *st,
    double (*h_fn)(char *data,double x,double y),
    char *data
  )
{
  bool read_grid_from_file;

  if (strcmp(par->spoints, "") != 0)
     read_grid_from_file = TRUE;
  else
     read_grid_from_file = FALSE;
  if (read_grid_from_file == TRUE)
    size_of_grid_in_file(par->spoints, &par->num_cols, &par->num_rows);
  *tg = malloc_tgrid(par->num_cols, par->num_rows);
  *eg = malloc_egrid(par->num_cols, par->num_rows);
  if (read_grid_from_file == TRUE)
    init_tg_from_file(*tg, par->spoints, par);
  else
    init_tg(par, *tg, h_fn,data);
  init_stage(st, par);
  init_eg(*tg, *eg);
}

/* Gives the user help on using the program */
void drac_help(void)
{
  printf("\n     This is a contour drawing program.  There are two ways to\n");
  printf("use it.  Either you can define the height function by changing\n");
  printf("h_fn(x, y) in the file ~awm/w/drac/drac.c, or else you can specify\n");
  printf("a file that contains the (x,y,z) coordinates of the spoints on a\n");
  printf("rectangular grid.\n");
  printf("     The program constructs a triangulation grid, finds the\n");
  printf("contours, displays them on the screen, and saves the diagram\n");
  printf("to a postscript file.\n");
  printf("     To run the program, type drac followed by any arguments:-\n");
  printf("-cols   the number of columns in the grid (default = %d)\n",
          NUM_COLUMNS);
  printf("-rows   the number of rows in the grid    (default = %d)\n",
          NUM_ROWS);
  printf("-nc     the number of contours to draw    (default = %d)\n",
          NUM_CONTOURS);
  printf("-res    the file to which to save the postscript (default = contour.ps)\n");
  printf("-title  the title for the diagram         (default = CONTOURS)\n");
  printf("-xlow   the minimum x value to be mapped  (default = %g)\n",
         DEFAULT_BL_DOMAIN_X);
  printf("-ylow   the minimum y value to be mapped  (default = %g)\n",
         DEFAULT_BL_DOMAIN_Y);
  printf("-xhigh  the maximum x value to be mapped  (default = %g)\n",
         DEFAULT_TR_DOMAIN_X);
  printf("-yhigh  the maximum y value to be mapped  (default = %g)\n",
         DEFAULT_TR_DOMAIN_Y);
  printf("E.g. 'drac -nc 4' would produce a diagram with 4 contours\n");
  printf("     If you use a file to give the (x, y, z) coordinates, the\n");
  printf("spoints must be given in the following order: scan along the\n");
  printf("rows from the bottom row (minimum y) to the top row, going\n");
  printf("from the minimum x to the maximum x along each row.  For each\n");
  printf("spoint, give the x, y, z values in turn.  The spoints must lie\n");
  printf("on a rectangular grid.  Values for the number of columns, the\n");
  printf("minimum x value to be mapped, etc, are then deduced from the\n");
  printf("file and override any values that were specified in the\n");
  printf("command line.\n");
}


void fprintf_dproj(FILE *s, char *m1, dproj *dp, char *m2)
{
  char buff[100];
  fprintf(s,"%s -> screen_x_min = %g%s",m1,dp->screen_x_min,m2);
  fprintf(s,"%s -> screen_x_max = %g%s",m1,dp->screen_x_max,m2);
  fprintf(s,"%s -> screen_y_min = %g%s",m1,dp->screen_y_min,m2);
  fprintf(s,"%s -> screen_y_max = %g%s",m1,dp->screen_y_max,m2);
  fprintf(s,"%s -> indim = %d%s",m1,dp->indim,m2);
  sprintf(buff,"%s -> in_min",m1);
  fprintf_realnums(s,buff,dp->in_min,dp->indim,m2);
  sprintf(buff,"%s -> in_max",m1);
  fprintf_realnums(s,buff,dp->in_max,dp->indim,m2);
  fprintf(s,"%s -> height_min = %g%s",m1,dp->height_min,m2);
  fprintf(s,"%s -> height_max = %g%s",m1,dp->height_max,m2);
}

double ag_value(
    double x,
    double x_min,
    double x_max,
    double screen_min,
    double screen_max
  )
{
  return(screen_min + 
         (screen_max - screen_min) * (x - x_min) / (x_max - x_min)
        );
}

/* 17 Oct 96: Mary deleted unreference formal parameter "double height" */
double dproj_to_ag_x(dproj *dp,double *farr)
{
  if ( dp->indim < 1 || dp -> indim > 2 )
    my_error("drac.c :: dprooj_to_ag_x: only do this for 1 or 2 dimensions");
  return(ag_value(farr[0],dp->in_min[0],dp->in_max[0],
                  dp->screen_x_min,dp->screen_x_max
                 )
        );
}

double dproj_to_ag_y(dproj *dp,double *farr,double height)
{
  double y,y_min,y_max;
  double result;
  if ( dp->indim < 1 || dp -> indim > 2 )
    my_error("drac.c :: dprooj_to_ag_y: only do this for 1 or 2 dimensions");
  
  if ( dp->indim == 1 )
  {
    y = height;
    y_min = dp->height_min;
    y_max = dp->height_max;
  }
  else
  {
    y = farr[1];
    y_min = dp->in_min[1];
    y_max = dp->in_max[1];
  }

  result = ag_value(y,y_min,y_max,dp->screen_y_min,dp->screen_y_max);
  return(result);
}

void dproj_dot(dproj *dp,double *farr,double hgt)
{
/*  fprintf_realnums(stdout,"dproj_dot, farr = ",farr,dp->indim,"\n"); */
/* Correct version :
  ag_dot(dproj_to_ag_x(dp,farr),dproj_to_ag_y(dp,farr,hgt)); */
/* Whizzy version : */
  int old_color = ag_pen_color();
  ag_set_pen_color(AG_RED);
  ag_disc(dproj_to_ag_x(dp,farr),dproj_to_ag_y(dp,farr,hgt),3.0);
  ag_set_pen_color(AG_DISPLAY);
  ag_circle(dproj_to_ag_x(dp,farr),dproj_to_ag_y(dp,farr,hgt),3.0);
  ag_set_pen_color(old_color);
}

void dproj_simple_dot(dproj *dp,double *farr,double hgt)
{
  ag_dot(dproj_to_ag_x(dp,farr),dproj_to_ag_y(dp,farr,hgt));
}

void dproj_line(dproj *dp,double *farr1,double hgt1,double *farr2,double hgt2)
{
  ag_line(dproj_to_ag_x(dp,farr1),dproj_to_ag_y(dp,farr1,hgt1),
          dproj_to_ag_x(dp,farr2),dproj_to_ag_y(dp,farr2,hgt2)
         );
}

void dproj_rectangle(
    dproj *dp,
    double *farr1,
    double hgt1,
    double *farr2,
    double hgt2
  )
{
  ag_rectangle(dproj_to_ag_x(dp,farr1),dproj_to_ag_y(dp,farr1,hgt1),
               dproj_to_ag_x(dp,farr2),dproj_to_ag_y(dp,farr2,hgt2)
              );
}

void dproj_box(
    dproj *dp,
    double *farr1,
    double hgt1,
    double *farr2,
    double hgt2
  )
{
  ag_box(dproj_to_ag_x(dp,farr1),dproj_to_ag_y(dp,farr1,hgt1),
         dproj_to_ag_x(dp,farr2),dproj_to_ag_y(dp,farr2,hgt2)
        );
}

void dproj_string(dproj *dp,double *farr,double hgt,char *string)
{
  ag_print(dproj_to_ag_x(dp,farr),dproj_to_ag_y(dp,farr,hgt),
           string
          );
}

void dproj_small_cross(dproj *dp,double *farr,double hgt)
{
  double x = dproj_to_ag_x(dp,farr);
  double y = dproj_to_ag_y(dp,farr,hgt);
  double r = 5.0;

  ag_line(x-r,y-r,x+r,y+r);
  ag_line(x-r,y+r,x+r,y-r);
}

void dproj_small_circle(dproj *dp,double *farr,double hgt)
{
  double x = dproj_to_ag_x(dp,farr);
  double y = dproj_to_ag_y(dp,farr,hgt);
  double r = 5.0;

  ag_circle(x,y,r);
}

void dproj_circle(dproj *dp,double *farr,double hgt,double radius)
{
  double x = dproj_to_ag_x(dp,farr);
  double y = dproj_to_ag_y(dp,farr,hgt);
  double farr2[2];
  double r;

  farr2[0] = farr[0] + radius;
  farr2[1] = farr[1];

  r = dproj_to_ag_x(dp,farr2) - x;

  ag_circle(x,y,r);
}


void dproj_dyv_dot(dproj *dp,dyv *d,double hgt)
{
  double *farr = mk_farr_from_dyv(d);
  dproj_dot(dp,farr,hgt);
  am_free_realnums(farr,dyv_size(d));
}

void dproj_dyv_line(dproj *dp,dyv *d1,double hgt1,dyv *d2,double hgt2)
{
  double *farr1 = mk_farr_from_dyv(d1);
  double *farr2 = mk_farr_from_dyv(d2);
  dproj_line(dp,farr1,hgt1,farr2,hgt2);
  am_free_realnums(farr1,dyv_size(d1));
  am_free_realnums(farr2,dyv_size(d2));
}

void dproj_dyv_rectangle(dproj *dp,dyv *d1,double hgt1,dyv *d2,double hgt2)
{
  double *farr1 = mk_farr_from_dyv(d1);
  double *farr2 = mk_farr_from_dyv(d2);
  dproj_rectangle(dp,farr1,hgt1,farr2,hgt2);
  am_free_realnums(farr1,dyv_size(d1));
  am_free_realnums(farr2,dyv_size(d2));
}

void dproj_dyv_box(dproj *dp,dyv *d1,double hgt1,dyv *d2,double hgt2)
{
  double *farr1 = mk_farr_from_dyv(d1);
  double *farr2 = mk_farr_from_dyv(d2);
  dproj_box(dp,farr1,hgt1,farr2,hgt2);
  am_free_realnums(farr1,dyv_size(d1));
  am_free_realnums(farr2,dyv_size(d2));
}

void dproj_dyv_string(dproj *dp,dyv *d,double hgt,char *string)
{
  double *farr = mk_farr_from_dyv(d);
  dproj_string(dp,farr,hgt,string);
  am_free_realnums(farr,dyv_size(d));
}

void dproj_dyv_small_cross(dproj *dp,dyv *dv,double hgt)
{
  double *farr = mk_farr_from_dyv(dv);
  dproj_small_cross(dp,farr,hgt);
  am_free_realnums(farr,dyv_size(dv));
}

void dproj_dyv_small_circle(dproj *dp,dyv *dv,double hgt)
{
  double *farr = mk_farr_from_dyv(dv);
  dproj_small_circle(dp,farr,hgt);
  am_free_realnums(farr,dyv_size(dv));
}

void dproj_dyv_circle(dproj *dp,dyv *dv,double hgt,double radius)
{
  double *farr = mk_farr_from_dyv(dv);
  dproj_circle(dp,farr,hgt,radius);
  am_free_realnums(farr,dyv_size(dv));
}

void draw_a_height(
    dproj *dp,
    double (*height_function)(char *data,  double x, double y),
    char *data,
    double *farr
  )
{
  char buff[100];
  double height = height_function(data,farr[0],farr[1]);
  if ( fabs(height) < 1e-5 ) height = 0.0;

  sprintf(buff,"%.3g",height);
/*  dproj_dot(dp,farr,height); */
  dproj_string(dp,farr,height,buff);
}

void draw_nine_heights(
    dproj *dp,
    double (*height_function)(char *data,  double x, double y),
    char *data
  )
{
  int i,j;
  for ( i = -1 ; i <= 1 ; i++ )
    for ( j = -1 ; j <= 1 ; j++ )
    {
      double farr[2];
      farr[0] = (dp->in_min[0] + dp->in_max[0]) / 2.0 +
                0.45 * i * (dp->in_max[0] - dp->in_min[0]);
      farr[1] = (dp->in_min[1] + dp->in_max[1]) / 2.0 +
                0.45 * j * (dp->in_max[1] - dp->in_min[1]);
      draw_a_height(dp,height_function,data,farr);
    }
}

void basic_draw_2d_function(
    bool draw_borders_only,
    double (*height_function)(char *data, double x, double y),
    char *data,
    int num_contours,
    int grid_size,
    char *title,
    char *xlabel,
    char *ylabel,
    double x_low,
    double y_low,
    double x_high,
    double y_high,
    dproj *result_dproj
  )
/*
   *result_dproj MUST NOT be NULL
*/
{
  int fake_argc = 0;
  char **fake_argv = NULL;
  drac_params pa;
  tgrid *tg;
  egrid *eg;
  stage st;

  drac_params_from_args(&pa,fake_argc,fake_argv);

  pa.num_cols = grid_size;              /* Number of columns in the grids */
  pa.num_rows = grid_size;              /* Number of rows in the grids */
  pa.num_contours = num_contours;       /* The number of contours to draw */
  pa.x_low = x_low;
  pa.x_high = x_high;
  pa.y_low = y_low;
  pa.y_high = y_high;                   /* Domains for x and y */

  init_structs_from_drac_params(&pa,&tg,&eg,&st,height_function,data);

  frame_and_title_stage(&st,
                        (title == NULL) ? "Contours" : title,
                        (xlabel == NULL) ? "" : xlabel,
                        (ylabel == NULL) ? "" : ylabel);

  result_dproj -> screen_x_min    = st.bl_window.x;
  result_dproj -> screen_y_min    = st.bl_window.y;
  result_dproj -> screen_x_max    = st.tr_window.x;
  result_dproj -> screen_y_max    = st.tr_window.y;
  result_dproj -> indim           = 2;
  result_dproj -> in_min[0]       = st.bl_domain.x;
  result_dproj -> in_max[0]       = st.tr_domain.x;
  result_dproj -> in_min[1]       = st.bl_domain.y;
  result_dproj -> in_max[1]       = st.tr_domain.y;
  result_dproj -> height_min      = tgrid_min(tg);
  result_dproj -> height_max      = tgrid_max(tg);

  if ( Verbosity > 15.55 )
    fprintf_dproj(stdout,"dproj",result_dproj,"\n");

  if ( !draw_borders_only )
  {
    draw_contours(&st, tg, eg, pa.num_contours,NULL,NULL);
    draw_nine_heights(result_dproj,height_function,data);
  }

  free_egrid(eg);
  free_tgrid(tg);
}

void draw_2d_function(
    double (*height_function)(char *data, double x, double y),
    char *data,
    int num_contours,
    int grid_size,
    char *title,
    double x_low,
    double y_low,
    double x_high,
    double y_high,
    dproj *result_dproj
  )
/*
   *result_dproj may be NULL, denoting that we don't want dproj info.
*/
{
  dproj spare_dproj[1];
  dproj *dp;

  if ( result_dproj == NULL )
    dp = spare_dproj;
  else
    dp = result_dproj;

  basic_draw_2d_function(FALSE,
                         height_function,
                         data,
                         num_contours,
                         grid_size,
                         title,
                         NULL,
                         NULL,
                         x_low,
                         y_low,
                         x_high,
                         y_high,
                         dp
                        );
}


double dummy_1d_height_fn(char *data,double x)
{
  return(0.0);
}

void simple_draw_2d_function(
    double (*height_function)(char *data, double x, double y),
    char *data,
    double x_low,
    double y_low,
    double x_high,
    double y_high
  )
{
  draw_2d_function(height_function,data,
                   10,21,(char *)NULL,
                   x_low,y_low,x_high,y_high,NULL);
}

void draw_1d_function(
    double (*height_function)(char *data, double x),
    char *data,
    int grid_size,
    char *xlabel,
    char *ylabel,
    double x_low,
    double x_high,
    dproj *result_dproj
  )
/*
  result_dproj may be null in which case it is ignored
*/
{
  int num_spoints = grid_size+1;
  int i;
  ongr on[1];
  frame fr[1];

  double *x_arr = am_malloc_realnums(num_spoints);
  double *y_arr = am_malloc_realnums(num_spoints);

  for ( i = 0 ; i < num_spoints ; i++ )
  {
    double z = i / (double) (num_spoints-1);
    double x = x_low + (x_high - x_low) * z;
    double y = height_function(data,x);
    x_arr[i] = x;
    y_arr[i] = y;
  }

  full_screen_frame(fr);
  clear_ongr(on);
  compute_axis_limits(x_arr,num_spoints,&on->x_axis);
  compute_axis_limits(y_arr,num_spoints,&on->y_axis);
  compute_axes_details(on,fr);
  if (xlabel) sprintf(on->x_axis.label,"%s",xlabel);
  if (ylabel) sprintf(on->y_axis.label,"%s",ylabel);
  draw_axes(fr,on);
  plot_in_frame(fr,on,x_arr,y_arr,num_spoints,"LN");

  if ( result_dproj != NULL )
  {
    result_dproj -> screen_x_min    = on->x_axis.start.x;
    result_dproj -> screen_y_min    = on->y_axis.start.y;
    result_dproj -> screen_x_max    = on->x_axis.end.x;
    result_dproj -> screen_y_max    = on->y_axis.end.y;
    result_dproj -> indim           = 1;
    result_dproj -> in_min[0]       = on->x_axis.lo;
    result_dproj -> in_max[0]       = on->x_axis.hi;
    result_dproj -> height_min      = on->y_axis.lo;
    result_dproj -> height_max      = on->y_axis.hi;
  }
      
  am_free_realnums(x_arr,num_spoints);
  am_free_realnums(y_arr,num_spoints);
}

double dummy_height_fn(char *data,double x,double y)
{
  return(0.0);
}

void draw_2d_border(dproj *dp,double xlo,double ylo,double xhi,double yhi)
{
  basic_draw_2d_function(TRUE,
                         dummy_height_fn,
                         (char *)NULL,
                         3,
                         3,
                         "",
                         NULL,
                         NULL,
                         xlo,
                         ylo,
                         xhi,
                         yhi,
                         dp
                        );
}


void simple_draw_1d_function(
    double (*height_function)(char *data, double x),
    char *data,
    double x_low,
    double x_high
  )
{
  dproj dummy[1];
  draw_1d_function(height_function,data,
                   200,"x","f(x)",x_low,x_high,dummy
                  );
}

double test_height_fn(char *data,double x,double y)
{
  double result = (x / 100) * (x / 100) * sin(y);
  return(result);
}

/* =================================================================== */
/* Surgraphs --- 3D surface graphs defined on a regular grid           */
/* =================================================================== */

/* Make a surgraph with a regular (x,z) grid whose heights are given
   by points.  */
surgraph *mk_surgraph(dym *points)
{
  surgraph *s = AM_MALLOC(surgraph);
  s->points = mk_copy_dym(points);
  s->title = NULL;
  s->x = mk_empty_laxis();
  s->y = mk_empty_laxis();
  s->z = mk_empty_laxis();
  s->dots = mk_empty_dym_array();
  s->num_contours = 14;
  deduce_surgraph_height_limits(s);
  return(s);
}

/* Makes a surgraph with a grid of size 0x0.  This surgraph can
   still be used to store and display irregularly positioned
   datapoints (by calling surgraph_include_dot). */
surgraph *mk_surgraph_no_points()
{
  surgraph *s = AM_MALLOC(surgraph);
  s->points = mk_dym(0, 0);
  s->title = NULL;
  s->x = mk_empty_laxis();
  s->y = mk_empty_laxis();
  s->z = mk_empty_laxis();
  s->dots = mk_empty_dym_array();
  s->num_contours = 14;
  deduce_surgraph_height_limits(s);
  return(s);
}

void free_surgraph(surgraph *s)
{
  if (s->points != NULL ) free_dym(s->points);
  if (s->dots != NULL ) free_dym_array(s->dots);
  if (s->title != NULL ) free_string(s->title);
  if (s->x != NULL ) free_laxis(s->x);
  if (s->y != NULL ) free_laxis(s->y);
  if (s->z != NULL ) free_laxis(s->z);
  AM_FREE(s, surgraph);
}

surgraph *mk_copy_surgraph(surgraph *s)
{
  surgraph *ns = AM_MALLOC(surgraph);

  ns->points = (s->points==NULL) ? NULL : mk_copy_dym(s->points);
  ns->dots = (s->dots==NULL) ? NULL : mk_copy_dym_array(s->dots);
  ns->title = (s->title==NULL) ? NULL : mk_copy_string(s->title);
  ns->x = (s->x==NULL) ? NULL : mk_copy_laxis(s->x);
  ns->y = (s->y==NULL) ? NULL : mk_copy_laxis(s->y);
  ns->z = (s->z==NULL) ? NULL : mk_copy_laxis(s->z);
  ns->num_contours = s->num_contours;
  return(ns);
}

void fprintf_surgraph(FILE *s,char *m1,surgraph *x,char *m2)
{
  char buff[100];

  sprintf(buff,"%s -> title",m1);
  fprintf_string(s,buff,x->title,m2);
  sprintf(buff,"%s -> points",m1);
  fprintf_dym(s,buff,x->points,m2);
  sprintf(buff,"%s -> x",m1);
  fprintf_laxis(s,buff,x->x,m2);
  sprintf(buff,"%s -> y",m1);
  fprintf_laxis(s,buff,x->y,m2);
  sprintf(buff,"%s -> z",m1);
  fprintf_laxis(s,buff,x->z,m2);
  sprintf(buff,"%s -> dots",m1);
  fprintf_dym_array(s,buff,x->dots,m2);
  fprintf(s,"%s -> num_contours = %d%s",m1, x->num_contours, m2);
}

void set_surgraph_title(surgraph *sg,char *label)
{
  if (sg->title != NULL) free_string(sg->title);
  sg->title = mk_copy_string(label);
}

void set_surgraph_x_min(surgraph *sg,double xmin)
{
  sg->x->min_defined = TRUE;
  sg->x->min_value = xmin;
}

void set_surgraph_x_max(surgraph *sg,double xmax)
{
  sg->x->max_defined = TRUE;
  sg->x->max_value = xmax;
}

void set_surgraph_x_axis_label(surgraph *sg,char *label)
{
  maybe_replace(&(sg->x->name),label);
}

void set_surgraph_y_min(surgraph *sg,double ymin)
{
  sg->y->min_defined = TRUE;
  sg->y->min_value = ymin;
}

void set_surgraph_y_max(surgraph *sg,double ymax)
{
  sg->y->max_defined = TRUE;
  sg->y->max_value = ymax;
}

void set_surgraph_y_axis_label(surgraph *sg,char *label)
{
  maybe_replace(&(sg->y->name),label);
}

void set_surgraph_z_min(surgraph *sg,double zmin)
{
  sg->z->min_defined = TRUE;
  sg->z->min_value = zmin;
}

void set_surgraph_z_max(surgraph *sg,double zmax)
{
  sg->z->max_defined = TRUE;
  sg->z->max_value = zmax;
}

void set_surgraph_z_axis_label(surgraph *sg,char *label)
{
  maybe_replace(&(sg->z->name),label);
}

void set_surgraph_num_contours(surgraph *sg, int num_contours)
{
  if (num_contours < 0) 
    my_error("set_surgraph_num_contours: num_contours cannot be negative"); 
  sg->num_contours = num_contours;
}

int surgraph_num_contours(surgraph *sg)
{
  return(sg->num_contours);
}

int surgraph_num_dot_styles(surgraph *sg)
{
  return(dym_array_size(sg->dots));
}

void increase_surgraph_num_styles(surgraph *sg)
{
  dym *zero_by_three = mk_dym(0,3);

  if ( dym_array_size(sg->dots) > SURGRAPH_MAX_DOT_STYLES )
    my_error("Must be fewer than SURGRAPH_MAX_DOT_STYLES surgraph dot styles");

  add_to_dym_array(sg->dots,zero_by_three);
  free_dym(zero_by_three);
}

dym *surgraph_dots_dym_ref(surgraph *sg,int dot_style)
{
  if ( dot_style < 0 || dot_style >= SURGRAPH_MAX_DOT_STYLES )
    my_error("dot_style < 0 or dot_style >= SURGRAPH_MAX_DOT_STYLES");
  while ( surgraph_num_dot_styles(sg) <= dot_style )
    increase_surgraph_num_styles(sg);
  return(dym_array_ref(sg->dots,dot_style));
}

/* Add the datapoint with height ht at (x, z) to be displayed as a 
   separate point (i.e. not part of the sg->points grid from which the 
   interpolation is performed). 

   PRE: 0 <= style < SURGRAPH_MAX_DOT_STYLES. See the comment in drac.h for the meaning of dot styles 
*/
void surgraph_include_dot_with_style(surgraph *sg, double x, double z, double ht,int style)
{
  dym *dots = surgraph_dots_dym_ref(sg,style);

  add_row(dots);
  dym_set(dots, dym_rows(dots)-1, 0, x);
  dym_set(dots, dym_rows(dots)-1, 1, ht);
  dym_set(dots, dym_rows(dots)-1, 2, z);
}

/* Add the datapoint with height ht at (x, z) to be displayed as a 
   separate point (i.e. not part of the sg->points grid from which the 
   interpolation is performed). */
void surgraph_include_dot(surgraph *sg, double x, double z, double ht)
{
  surgraph_include_dot_with_style(sg,x,z,ht,0);
}

/* Returns the height (i.e. y value) of the (x_index, z_index)th point
   in the surgraph points grid.  It is an error to call it with
   (x_index, z_index) out of range. */
double surgraph_point_height(surgraph *sg, int x_index, int z_index)
{
  dym *points = sg->points;

  if (x_index < 0 || x_index >= dym_rows(points) ||
      z_index < 0 || z_index >= dym_cols(points))
    my_error("surgraph_point_height: indices out of range");
  return(dym_ref(points, x_index, z_index));
}

/* Sets the height (i.e. y value) of the (x_index, z_index)th point
   in the surgraph points grid.  It is an error to call it with
   (x_index, z_index) out of range. */
void set_surgraph_point_height(surgraph *sg, int x_index, 
                               int z_index, double height)
{
  dym *points = sg->points;

  if (x_index < 0 || x_index >= dym_rows(points) ||
      z_index < 0 || z_index >= dym_cols(points))
    my_error("set_surgraph_point_height: indices out of range");
  dym_set(points, x_index, z_index, height);
}

/* Sets the y laxis limits to reflect the minimum/maximum heights in sg.  */
void deduce_surgraph_height_limits(surgraph *sg)
{
  double min_ht = 100.0;
  double max_ht = -100.0;
  int i, j;
  dym *points = sg->points;
  bool seen_first_pt = FALSE;

  for (i = 0; i < dym_rows(points); i++)
    for (j = 0; j < dym_cols(points); j++)
    {
      double ht = dym_ref(points, i, j); 
      if (i + j == 0)
      {
        /* Looking at the very first point. */
        min_ht = ht;
        max_ht = min_ht;
        seen_first_pt = TRUE;
      }
      if (ht < min_ht) min_ht = ht;
      if (ht > max_ht) max_ht = ht;
    }

#ifdef NEVER
  for (i = 0; i < dym_rows(dots); i++)
  {
    double ht = dym_ref(dots, i, 1); 
    if (i == 0 && !seen_first_pt)
    {
      /* Looking at the very first point. */
      min_ht = ht;
      max_ht = min_ht;
      seen_first_pt = TRUE;
    }
    if (ht < min_ht) min_ht = ht;
    if (ht > max_ht) max_ht = ht;
  }
#endif

  set_surgraph_y_min(sg, min_ht);
  set_surgraph_y_max(sg, max_ht);
}

/* --------------------------------------------------------- */
/* Drawing contours from a surgraph                          */
/* --------------------------------------------------------- */

/* This initializes the values in the tg grid from the dym points.  In
   the special case where points is NULL or has too few elements we set
   all heights in tg to be zero.  This is useful when we are just
   going to produce a frame without drawing any contours. 

   Note that plenty of confusion is likely to result from the
   fact that tgrid calls the 2nd horizontal direction the y-axis,
   but surgraphs call this the z axis.  

   tg->num_rows is the number of rows in the tg grid, i.e. the number
   of distinct z values in the grid.  

   tg->num_cols is the number of columns in the tg grid, i.e. the
   number of distinct x values in the grid.
*/
void init_tg_from_surgraph(tgrid *tg, surgraph *sg, drac_params *par)
{
  int num_x_indices = tg->num_cols;
  int num_z_indices = tg->num_rows;
  int x_index, z_index;
  double tpx, tpz;
  tgpoint *tp;
  bool use_points = TRUE;

  if ((sg->points == NULL) || 
      (dym_rows(sg->points) < num_x_indices) || 
      (dym_cols(sg->points) < num_z_indices))
  {
    use_points = FALSE;
  }

  if ((num_x_indices <= 2) || (num_z_indices <= 2))
    my_error("Too few columns or rows for a tg_grid");
  for (x_index = 0; x_index < num_x_indices; x_index++)
    for (z_index = 0; z_index < num_z_indices; z_index++)
    {
      tp = tg_ref(tg, x_index, z_index);
      tp->tg_i = x_index;
      tp->tg_j = z_index;
      tp->index = x_index + z_index * num_x_indices;
      tpx = par->x_low + 
            (x_index / 
              ((double) (num_x_indices - 1))) * (par->x_high - par->x_low);
      tpz = par->y_low + 
            (z_index / 
              ((double) (num_z_indices - 1))) * (par->y_high - par->y_low);
      (tp->spnt).x = tpx;
      (tp->spnt).y = tpz;
      tp->height = (use_points) ? 
                     surgraph_point_height(sg, x_index, z_index) : 0.0;
    }
}

/* This makes an apict from a surgraph.

   sg->points defines a regular grid in the horizontal x, z
   dimensions, with the heights (y coordinates) of the (xi, zi)th grid
   point being dym_ref(sg->points, xi, zi).  Provided the (x, z) grid is at
   least 3 x 3, mk_apict_from_surgraph plots contours through the (x, y, z)
   coordinates defined by this grid.

   Whether or not sg->points is used, mk_apict_from_surgraph will plot
   a dot at any (x, y, z) datapoint specified in sg->dots.
   
   It is a horrible piece of ugliness that y represents the height
   dimension, and x and z are the horizontal directions. 

   It is another horrible piece of ugliness that the number of rows
   in sg->points is the number of distinct x values in the grid,
   but in the drac_params and tgrid structures this is called the
   number of columns.  ARGH!

   With all this confusion, I think I'll add some more (JGS 11-11-96)
   The data points were not actually showing up when requested.  I noted
   that switching the data point drawing section back to the x,y->z
   form rather than its previous (and consistent with the above comments)
   x,z->y form made the dots again appear correctly.  I think the code
   that built the dots structure already swapped y and z for us, but I
   couldn't be sure.

   Used to use the min and max height out of all the points given it as
   the range on the contours.  Change to look at the axis limit values
   in the surgraph (this enable user specified limits to affect the graph)

   AWM - The dots now come up in one of four different colors depending on their style
*/
apict *mk_apict_from_surgraph(surgraph *sg)
{
  apict *ap;
  int fake_argc = 0;
  char **fake_argv = NULL;
  drac_params pa;
  tgrid *tg;
  egrid *eg;
  stage st;
  dproj dp[1];
  int num_contours = sg->num_contours;
  int num_x_indices = dym_rows(sg->points); 
  int num_z_indices = dym_cols(sg->points);
  char *title = sg->title;
  char *xlabel = sg->x->name;
  char *zlabel = sg->z->name;
  dyv *dv = mk_dyv_2(0.0, 0.0);
  int i;
  int saved_ag_color;

  apict_on();
  drac_params_from_args(&pa,fake_argc,fake_argv);

  pa.num_cols = (num_x_indices >= 3) ? num_x_indices : 3; /* #columns in the grids */
  pa.num_rows = (num_z_indices >= 3) ? num_z_indices : 3; /* #rows in the grids */
  pa.num_contours = num_contours;       /* The number of contours to draw */
  pa.x_low = (sg->x->min_defined) ? sg->x->min_value : 0.0;
  pa.x_high = (sg->x->max_defined) ? sg->x->max_value : 1.0;
  pa.y_low = (sg->z->min_defined) ? sg->z->min_value : 0.0;
  pa.y_high = (sg->z->max_defined) ? sg->z->max_value : 1.0;

  tg = malloc_tgrid(pa.num_cols, pa.num_rows);
  eg = malloc_egrid(pa.num_cols, pa.num_rows);
  init_tg_from_surgraph(tg, sg, &pa);
  init_tg_neighs(tg);
  init_stage(&st, &pa);
  init_eg(tg, eg);
  frame_and_title_stage(&st,
                        (title == NULL) ? "3DGraph" : title,
                        (xlabel == NULL) ? "" : xlabel,
                        (zlabel == NULL) ? "" : zlabel);

  dp -> screen_x_min    = st.bl_window.x;
  dp -> screen_y_min    = st.bl_window.y;
  dp -> screen_x_max    = st.tr_window.x;
  dp -> screen_y_max    = st.tr_window.y;
  dp -> indim           = 2;
  dp -> in_min[0]       = st.bl_domain.x;
  dp -> in_max[0]       = st.tr_domain.x;
  dp -> in_min[1]       = st.bl_domain.y;
  dp -> in_max[1]       = st.tr_domain.y;
  dp -> height_min      = tgrid_min(tg);
  dp -> height_max      = tgrid_max(tg);

  if (num_x_indices >= 3 && num_z_indices >=3)
  {
    draw_contours(&st, tg, eg, pa.num_contours,
                  (sg->y->min_defined) ? &sg->y->min_value : NULL,
                  (sg->y->max_defined) ? &sg->y->max_value : NULL);
  }

  saved_ag_color = ag_pen_color();

  /* Now draw the points in dots... */
  for (i = 0; i < surgraph_num_dot_styles(sg); i++)
  {
    int j;
    dym *dots = surgraph_dots_dym_ref(sg,i);
    int k = i % 5;
    int render_color = (k==0)?AG_RED:(k==1)?AG_BLUE:(k==2)?AG_GREEN:(k==3)?AG_BLACK:AG_PURPLE;
    ag_set_pen_color(render_color);
    for ( j = 0 ; j < dym_rows(dots) ; j++ )
    {
      double farr[2];
      farr[0] = dym_ref(dots,j,0);
      farr[1] = dym_ref(dots,j,1);
      dproj_simple_dot(dp, farr, dym_ref(dots, j, 2)); /* 23 Jan 97: corrected from
                                                          dym_ref(dots, i, 2) by Mary */
    }
  }

  ag_set_pen_color(saved_ag_color);

  free_dyv(dv);
  free_egrid(eg);
  free_tgrid(tg);
  ap = apict_off();
  return(ap);
}

surgraph *mk_test_surgraph()
{
  dym *points = mk_constant_dym(5, 10, 2.0);
  surgraph *sg;
  int i, j;

  sg = mk_surgraph(points);
  for (i = 0; i < 5; i++)
    for (j = 0; j < 10; j++)
    {
      set_surgraph_point_height(sg, i, j, 8 + 10 - j);
    }
  set_surgraph_title(sg, "A Tilted Plane");
  set_surgraph_y_axis_label(sg, "Height");
  set_surgraph_x_min(sg, 6);
  set_surgraph_x_max(sg, 9);
  set_surgraph_z_min(sg, -1);
  set_surgraph_z_max(sg, 10);
  set_surgraph_x_axis_label(sg, "Walnuts");
  set_surgraph_z_axis_label(sg, "Peanuts");
  set_surgraph_num_contours(sg, 3);
  deduce_surgraph_height_limits(sg);
  surgraph_include_dot(sg, 6, 1, 12);
  surgraph_include_dot(sg, 6, 5, 16);
  surgraph_include_dot(sg, 9, 9, 8);
  surgraph_include_dot(sg, 8, 0.5, 17);
  free_dym(points);

  return(sg);
}

void render_surgraph(surgraph *sg)
{
  apict *ap = mk_apict_from_surgraph(sg);
  render_apict(ap);
  free_apict(ap);
}

/* --------------------------------------------------------- */
/* Making a surgraph from a 2d function                      */
/* --------------------------------------------------------- */

/* 
   mk_surgraph_from_2d_function creates a surgraph on a grid
   of the given size (which must be at least 3), where the
   x,z coordinates have the given range, and the heights are
   given by the height_function.

   N.B.  It is a horrible piece of ugliness that in surgraphs y
   represents the height dimension, and x and z are the horizontal
   directions.  
*/
surgraph *mk_surgraph_from_2d_function(
    double (*height_function)(char *data, double x, double z),
    char *data,
    int num_contours,
    int grid_size,
    char *title,
    char *xlabel,
    char *zlabel,
    double x_low,
    double z_low,
    double x_high,
    double z_high
  )
{
  dym *points;
  surgraph *sg;
  char sg_title[200];

  if (grid_size >= 3)
  {
    int xi, zi;

    points = mk_constant_dym(grid_size, grid_size, -77);
    sg = mk_surgraph(points);
    for (xi = 0; xi < grid_size; xi++)
      for (zi = 0; zi < grid_size; zi++)
      {
        double tpx, tpz;

        tpx = x_low + 
              (xi / ((double) (grid_size - 1))) * (x_high - x_low);
        tpz = z_low + 
              (zi/ ((double) (grid_size - 1))) * (z_high - z_low);

        set_surgraph_point_height(sg, xi, zi,
                                  height_function((char *) data, tpx, tpz));
    }
  }
  else
  {
    /* Do not set up a grid from which to interpolate a surface */
    points = mk_dym(0, 0);
    sg = mk_surgraph(points);
  }
  free_dym(points);

  sprintf(sg_title, "%s%s%s%s%s%s",
          (title == NULL) ? "3DGraph" : title,
          (xlabel && zlabel) ? " of (" : "",
          (xlabel && zlabel) ? xlabel : "",
          (xlabel && zlabel) ? "," : "",
          (xlabel && zlabel) ? zlabel : "",
          (xlabel && zlabel) ? ")" : "");
  set_surgraph_title(sg, sg_title);
  set_surgraph_x_min(sg, x_low);
  set_surgraph_x_max(sg, x_high);
  set_surgraph_z_min(sg, z_low);  /* Remember horrible fact...*/
  set_surgraph_z_max(sg, z_high);
  deduce_surgraph_height_limits(sg);
  set_surgraph_num_contours(sg, num_contours);

  return(sg);
}


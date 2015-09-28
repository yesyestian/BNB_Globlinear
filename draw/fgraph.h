/*
   File:        fgraph.h
   Author:      Andrew W. Moore
   Created:     Sep 1st 1997
   Description: New Test Search Algorithm


   Copyright (C) 1997, Andrew Moore
*/

#ifndef FGRAPH_H
#define FGRAPH_H
#include "amgr.h"
#include "./utils/amut.h"

typedef struct faxis
{
  double lo,hi;
  double ag_lo,ag_hi;
} faxis;

typedef struct fgraph
{
  faxis x[1];
  faxis y[1];
} fgraph;

fgraph *mk_fgraph(double xlo,double ylo,double xhi,double yhi,
                  double ag_xlo,double ag_ylo,double ag_xhi,double ag_yhi);

void fg_dot(fgraph *fg,double x,double y);

void fg_circle(fgraph *fg,double x,double y,double ag_radius);

void fg_disc(fgraph *fg,double x,double y,double ag_radius);

void fg_box(fgraph *fg,double x,double y,double u,double v);

void fg_line(fgraph *fg,double x,double y,double u,double v);

void fg_string(fgraph *fg,double x,double y,char *s);

void fg_border(fgraph *fg);

int fg_get_xy(fgraph *fg,double *x,double *y);

void free_fgraph(fgraph *fg);

void fg_contours(fgraph *fg,
                     double height_function(char *data,double x,double y),
                     char *data,
                     int grid_size,
                     int num_contours);

/* Draws a line indicating the set of points for
   which (x - pivot).dir == 0 */
void fg_pivot_dir(fgraph *fg,dyv *pivot,dyv *dir);

#endif





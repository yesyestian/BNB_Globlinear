/*
   File:        drac.h
   Authors:     Andrew W. Moore, Mary Soon Lee
   Created:     Fri Jan 21 09:56:03 EST 1994
   Description: Header for DRAwing Contours
   Last edited: 7 October 96

   19 May 96: Surgraphs added.  Surgraphs are 3D surface graphs with
   points on a regular grid.

   Copyright (C) 1996, Mary Soon Lee
*/

#ifndef DRAC_H
#define DRAC_H

#include "./utils/ambs.h"      /* Very basic operations */
#include "amgr.h"      /* Basic (0,512)x(0,512) Graphics window */
#include "./utils/amma.h"      /* Fast, non-fragmenting, memory management */
#include "./utils/amar.h"      /* Obvious operations on 1-d arrays */
#include "ongr.h"      /* Simple 1-d graph drawing */
#include "./utils/amdym.h"
#include "apict.h"
#include "lingraph.h"  /* Needed for laxis structure */

#ifndef DOUBLE_FUNCTION_DEFINED
typedef double (*double_function)();
#define DOUBLE_FUNCTION_DEFINED
#endif

#define MAX_NEIGHS 6

typedef struct spoint
{
  double x;
  double y;
} spoint;

/* A spoint in the triangulation grid, complete with its neighbors.  */
typedef struct tgpoint_struct
{
  int tg_i;               /* i index into array */
  int tg_j;               /* j index into array */
  int index;              /* tg_i + tg_j * num_cols; a unique single
                             number index for the tgpoint, used in
                             finding the edges from it */
  spoint spnt;            /* the spoint in h_fn coordinates */
  double height;          /* the height of the spoint (from h_fn) */
  int num_neighs;         /* number of neighbours of spoint in the grid */
  struct tgpoint_struct *neighs[MAX_NEIGHS];
                          /* its neighbors in clockwise order */
  bool loop;              /* TRUE iff the neighbors form a complete loop.
                             i.e. spoint not on an edge. */
} tgpoint;

typedef struct tgrid_struct
{
  int num_cols,num_rows;
  tgpoint ***tgpoints;
} tgrid;

/* This defines the sub-window within the AG window in which contours
   will be drawn, plus the region over which to map the height function. */
typedef struct stage_struct
{
  spoint bl_window;         /* bottom left spoint of sub-window in AG coords */
  spoint tr_window;         /* top right spoint of sub-window in AG coords */
  spoint bl_domain;         /* bottom left spoint of region of h_fn to map */
  spoint tr_domain;         /* top right spoint of region of h_fn to map */
} stage;

/* This represents an edge, including spointers to the edges that are
   in the same triangles in the triangulation of the grid.  If the edge
   lies on the border of the grid, it will only be in one triangle,
   otherwise it will be in two.  */
typedef struct edge_struct
{
  bool is_edge;            /* TRUE iff represents a real edge */
  tgpoint *start, *end;    /* start and end spoints of the edge (order
                              arbitrary) */
  bool border_edge;        /* TRUE iff the edge is on the border of grid */
  struct edge_struct *tr1_e1, *tr1_e2;
                           /* the other two edges in the first adjoining
                              triangle */
  struct edge_struct *tr2_e1, *tr2_e2;
                           /* the other two edges in the second adjoining
                              triangle.  Invalid if the main edge is on
                              the border of the grid.  */
  bool visited;            /* for use when creating contours.  True iff
                              the edge has been visited during the process */
} edge;

/* This is an array of all the edges in a triangulation grid.  To find
   the edges from a particular tgpoint, tg, first find that tgpoint's
   index number, n = tg_j * num_rows + tg_i.  Then edges[n][i] holds the
   edge (if any) going from tg to tg->neighs[i].  */
typedef struct egrid_struct
{
  int num_cols, num_rows;   /* the number of cols + rows in the grid */
  edge ***edges;
} egrid;

typedef struct drac_params_struct
{
  int num_cols;     /* Number of columns in the grids */
  int num_rows;     /* Number of rows in the grids */
  int num_contours; /* The number of contours to draw */
  char rfile[100];  /* The file to which to send the postscript diagram */
  char spoints[100]; /* File from which to read the grid */
  char title[200];  /* The title for the diagram */
  double x_low, x_high, y_low, y_high; /* Domains for x and y */
} drac_params;

void drac_main(int argc, char **argv);   /* (argc, argv) */

typedef struct dproj_struct
{
  double screen_x_min;
  double screen_y_min;
  double screen_x_max;
  double screen_y_max;
  int indim;
  double in_min[2];
  double in_max[2];
  double height_min;
  double height_max;
} dproj;

void fprintf_dproj(FILE *s, char *m1, dproj *dp, char *m2);

void dproj_dot(dproj *dp,double *farr,double hgt);
void dproj_line(dproj *dp,double *farr1,double hgt1,double *farr2,double hgt2);
void dproj_rectangle(
    dproj *dp,
    double *farr1,
    double hgt1,
    double *farr2,
    double hgt2
  );
void dproj_box(
    dproj *dp,
    double *farr1,
    double hgt1,
    double *farr2,
    double hgt2
  );
void dproj_string(dproj *dp,double *farr,double hgt,char *string);
void dproj_small_cross(dproj *d,double *farr,double hgt);
void dproj_small_circle(dproj *d,double *farr,double hgt);
void dproj_circle(dproj *d,double *farr,double hgt,double radius);

void dproj_dyv_dot(dproj *dp,dyv *d,double hgt);
void dproj_dyv_line(dproj *dp,dyv *d1,double hgt1,dyv *d2,double hgt2);
void dproj_dyv_rectangle(dproj *dp,dyv *d1,double hgt1,dyv *d2,double hgt2);
void dproj_dyv_box(dproj *dp,dyv *d1,double hgt1,dyv *d2,double hgt2);
void dproj_dyv_string(dproj *dp,dyv *d,double hgt,char *string);
void dproj_dyv_small_cross(dproj *d,dyv *dv,double hgt);
void dproj_dyv_small_circle(dproj *d,dyv *dv,double hgt);
void dproj_dyv_circle(dproj *d,dyv *dv,double hgt,double radius);

/*
   Note: After you've drawn your 2-d graph or 1-d graph,
   you might desire to draw some
   graphics objects on it of your own choosing. That can be done.
   First you'll need a data structure which remembers the scaling that
   was used for drawing the graph. That's a dproj. Declare a dproj, and pass
   it in your code when you call a graphing function. Then, for example,
   to draw a dot at location x where x is a 1-d or 2-d vector and at
   height h, call dproj_dot(dp,x,h) where dp spoints to the the dproj.
*/

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
  );

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
  );

void draw_2d_border(dproj *dp,double xlo,double ylo,double xhi,double yhi);

void simple_draw_2d_function(
    double (*height_function)(char *data, double x, double y),
    char *data,
    double x_low,
    double y_low,
    double x_high,
    double y_high
  );

void draw_1d_function(
    double (*height_function)(char *data, double x),
    char *data,
    int grid_size,
    char *xlabel,
    char *ylabel,
    double x_low,
    double x_high,
    dproj *result_dproj
  );

void simple_draw_1d_function(
    double (*height_function)(char *data, double x),
    char *data,
    double x_low,
    double x_high
  );


/* =================================================================== */
/* Surgraphs --- 3D surface graphs defined on a regular grid           */
/* =================================================================== */
/*

   Surgraphs are the structures used to define 3D graphs.  When a surgraph
   is included in a report, a 2D contour-plot of the surgraph will be displayed
   in the main scrolling window.

   Typically a surgraph contains height values defined on a regular
   horizontal grid.  The 2D rendering of the surgraph draws contour
   lines to show the shape of the function.

   The programmer can also specify individual datapoints, not on the
   regular grid, which will be shown with a dot.

   An ugly fact
   ------------
   In surgraphs y represents the height dimension, and x and z are the
   horizontal directions.

   Creating a surgraph
   -------------------
   The following functions are provided to help programmers create
   surgraphs:

      surgraph *mk_surgraph(dym *points)
            the dym points specifies the values of the surgraph on a
            regular grid.  The height (y coordinate) of the point at
            the (i,j) horizontal grid position is given by
            dym_ref(points, i, j).

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
	  );
           creates a surgraph on a grid of the given size (which must
           be at least 3), where the x,z coordinates have the given
           range, and the heights are given by the height_function.
           num_contours is the number of contours to use in the
           2D rendering of the surgraph (note that this number will
           be treated as an approximate guide).

     surgraph *mk_surgraph_from_2d_agspec(sidat *si, facode *fc,
                                          char *datafile, agspec *ag)
           creates a surgraph based on the agspec [[Andrew: document]]

     surgraph *mk_surgraph_no_points()
           makes a surgraph with a 0x0 grid (this can still be used to
           store and display irregularly positioned datapoints by
           calling surgraph_include_dot)

     void surgraph_include_dot(surgraph *sg, double x, double z, double ht)
           adds the datapoint with height ht at (x, z) to be displayed
           as a separate point (i.e. not part of the regular grid
           from which the interpolation is performed).

           The dot_style of a dot included this way is zero.

     void surgraph_include_dot_with_style(surgraph *sg,double x,double z,double ht,int style);
           Includes a dot that MIGHT be rendered with different appearence.

           style may be any integer between 0 and SURGRAPH_MAX_DOT_STYLES. When a surgraph is rendered,
           different styles may come out looking different (e.g. different colors).
           If there are fewer styles renderable than actual styles the styles
           will come out in the form rendered style = style % (number renderable styles)


   Specifying extra information for a surgraph
   -------------------------------------------
   To specify a label for the x axis, use:
      void set_surgraph_x_axis_label(surgraph *sg,char *label);
   and similarly for the y and z axes.

   To title the surgraph use:
      void set_surgraph_title(surgraph *sg,char *label);

   To add an isolated point (not on the surface) to the surgraph use:
      void surgraph_include_dot(surgraph *sg, double x, double z, double ht);
   N.B.  It is legal to have no regular grid points at all, just datapoints
   added using surgraph_include_dot.

   To specify the number of contours to use when rendering the surgraph
   in apict form use:
      void set_surgraph_num_contours(surgraph *sg, int num_contours);


   Other interface functions for surgraphs
   ---------------------------------------
   To free a surgraph use:
      void free_surgraph(surgraph *x);

   To copy a surgraph use:
      surgraph *mk_copy_surgraph(surgraph *x);

   To make an apict corresponding to a surgraph use:
      apict *mk_apict_from_surgraph(surgraph *sg);

   To display a surgraph use:
      void render_surgraph(surgraph *sg);
   Note that the display mechanism may vary depending on the platform.

   To find the height of a point in the main grid:
      double surgraph_point_height(surgraph *sg, int x_index, int z_index);

   To set the height of a point in the main grid:
      void set_surgraph_point_height(surgraph *sg, int x_index,
                                     int z_index, double height);

   To find the number of contours for the surgraph:
      int surgraph_num_contours(surgraph *sg);

   To set the height range correctly based on the points in the surgraph:
      void deduce_surgraph_height_limits(surgraph *sg);
*/

#define SURGRAPH_MAX_DOT_STYLES 1000

typedef struct surgraph
{
  char *title; /* NULL denotes undefined */
  dym *points; /* The regular grid */
  laxis *x;    /* The x-axis limits for the grid-points. */
  laxis *y;    /* The y-axis limits for the grid-points. */
  laxis *z;    /* The z-axis limits for the grid-points. */
  dym_array *dots;  /* Individual datapoints. ith entry for ith style */
  int num_contours; /* # of contours to use when rendering as an apict */
} surgraph;

surgraph *mk_surgraph(dym *points);
surgraph *mk_surgraph_no_points();
void free_surgraph(surgraph *x);
surgraph *mk_copy_surgraph(surgraph *x);
void fprintf_surgraph(FILE *s,char *m1,surgraph *x,char *m2);
void set_surgraph_title(surgraph *sg,char *label);
void set_surgraph_x_axis_label(surgraph *sg,char *label);
void set_surgraph_x_min(surgraph *sg,double xmin);
void set_surgraph_x_max(surgraph *sg,double xmax);
void set_surgraph_y_axis_label(surgraph *sg,char *label);
void set_surgraph_y_min(surgraph *sg,double ymin);
void set_surgraph_y_max(surgraph *sg,double ymax);
void set_surgraph_z_axis_label(surgraph *sg,char *label);
void set_surgraph_z_min(surgraph *sg,double ymin);
void set_surgraph_z_max(surgraph *sg,double ymax);
void surgraph_include_dot(surgraph *sg, double x, double z, double ht);
void surgraph_include_dot_with_style(surgraph *sg,double x,double z,double ht,int style);

double surgraph_point_height(surgraph *sg, int x_index, int z_index);
void set_surgraph_point_height(surgraph *sg, int x_index,
                               int z_index, double height);
apict *mk_apict_from_surgraph(surgraph *sg);
void render_surgraph(surgraph *sg);
void set_surgraph_num_contours(surgraph *sg, int num_contours);
int surgraph_num_contours(surgraph *sg);
void deduce_surgraph_height_limits(surgraph *sg);

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
  );

#endif /* #ifndef DRAC_H */

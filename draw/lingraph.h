/*
   File:        lingraph.h
   Author:      Andrew W. Moore
   Created:     Wed May  8 19:16:35 EDT 1996
   Description: Header for Graphs with lines on them

   Modified:    sir - Fri Jun  9 13:32:50 EDT 2000
                APICT_FROM_LINGRAPH_IS_ACTUALLY_NEEDED_BY_SOMEONE
		idef was missing from .h.

   Copyright (C) 1996, A. W. Moore
*/


#ifndef LINGRAPH_H
#define LINGRAPH_H
#include "./utils/amut.h"
#include "apict.h"
#include "amgr.h"

#define JOINED_STYLE 0
#define DOTTED_STYLE 1
#define PIXEL_STYLE 2

typedef struct laxis
{
  char *name;   /* NULL denotes undefined */
  bool min_defined;
  double min_value;
  bool max_defined;
  double max_value;
} laxis;

typedef struct lingraph
{
  dym_array *lines;
    /* Number of lines = dym_array_size(lines).
       Number of points on jth line = dym_rows(dym_array_ref(lines,j))
       x coord of ith point on jth line = dym_ref(dym_array_ref(lines,j),i,0)
       y coord of ith point on jth line = dym_ref(dym_array_ref(lines,j),i,1)
    */
  ivec *color_codes;
  ivec *styles;
  string_array *labels;  /* string_array_ref(labels,i) = label for i'th line.
                             NULL denotes undefined */
  char *title; /* NULL denotes undefined */

  dym *string_locations;
  string_array *strings;
      /* The above two fields are used for strings which
         can appear within the actual diagram. */

  laxis *x;
  laxis *y;
} lingraph;

lingraph *mk_empty_lingraph();
int lingraph_num_lines(lingraph *lg);
void free_lingraph(lingraph *x);
lingraph *mk_copy_lingraph(lingraph *x);
void fprintf_lingraph(FILE *s,char *m1,lingraph *x,char *m2);
void add_to_lingraph(lingraph *lg,int line_num,double x,double y);
void add_string_to_lingraph(lingraph *lg,double x,double y,char *string);
void set_lines_same_color(lingraph *lg,int line_num1,int line_num2);
void set_line_style_dotted(lingraph *lg,int line_num);
void set_line_style_pixel(lingraph *lg,int line_num);
void set_x_axis_label(lingraph *lg,char *label);
void set_y_axis_label(lingraph *lg,char *label);
void set_lingraph_title(lingraph *lg,char *label);
void set_line_label(lingraph *lg,int line_num,char *label);
void set_x_min(lingraph *lg,double xmin);
void set_x_max(lingraph *lg,double xmax);
void set_y_min(lingraph *lg,double ymin);
void set_y_max(lingraph *lg,double ymax);

#ifdef APICT_FROM_LINGRAPH_IS_ACTUALLY_NEEDED_BY_SOMEONE
apict *mk_apict_from_lingraph(lingraph *lg);
#endif

void render_lingraph_in_agbox(agbox *agb,lingraph *lg);
void render_lingraph(lingraph *lg);

/*------------------------------------------
  This adds an ellipse 'linenum' to lingraph 'lg'
  ellipse is defined by its center coordinates 'x0' 'y0',
  half-axes lengths 'a' and 'b' and rotation angle 'phi'
  (phi in degrees!!!)
  'resolution' determines how many line segments will
  compose the ellipse. Put resolution=72 to get one line
  segment per a 5 degree angle slice (360/72 = 5).
  Coded by Artur Dubrawski, 03 Aug. 1996
------------------------------------------*/
extern void add_ellipse_to_lingraph( lingraph* /*lg*/, int /*linenum*/,
				     double /*x0*/, double /*y0*/,
				     double /*a*/, double /*b*/,
				     double /*phi*/,
				     int /*resolution*/ );

/* The following functions are NOT for general use. */
laxis *mk_empty_laxis();
void free_laxis(laxis *x);
laxis *mk_copy_laxis(laxis *x);
void fprintf_laxis(FILE *s,char *m1,laxis *x,char *m2);
void maybe_replace(char **x,char *label);

void add_dyv_to_lingraph(lingraph *lg,int line,dyv *dv);
lingraph *mk_lingraph_from_dyv(char *title,dyv *dv);
void draw_dyv(char *title,dyv *dv);

#endif /* #ifndef LINGRAPH_H */



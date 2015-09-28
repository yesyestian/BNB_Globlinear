
/*
   File:          ongr.h
   Author:        Andrew W. Moore
   Created:       August 19th 1992

   Description:   Some simple graph-drawing
*/

#ifndef ONGR_H
#define ONGR_H

#include "./utils/ambs.h" /* Header file for a set of basic convenience routines */
#include "amgr.h" /* Header file for Simple X and PostScript graphcs */
#include "./utils/amar.h" /* Simple array handling */
#include "./utils/amut.h"

typedef struct onpoint_struct
{
  double x,y;
} onpoint;

typedef struct axis_struct
{
/* The following are computed by clear_ongr() */
  bool vertical;
  char label[100];
/* label may be subsequently altered by user */

/* The following are computed by compute_limits() or may be set directly
   by user
*/
  double lowest_val,highest_val;

/* The following are computed by compute_axis_markings, and may then
   be altered by user
*/
  double lo,hi,delta;
  onpoint start,end;    /* In andrews-graphics coords */
  int marks;
  onpoint rel_mark_number; /* Where, relative to the axis mark, does the
                                 printed number go? In ag coords */
} axis;

typedef struct ongr_struct
{
  axis x_axis;
  axis y_axis;
} ongr;

typedef struct frame_struct
/* Draw graph relative to these coords */
{
  onpoint bottom_left;
  onpoint top_right;
} frame;

extern void clear_ongr( ongr *on );
/* Do this before you use ongr */

extern void compute_axis_limits( double *farr, int size, axis *ax );
/* Sets ax->low_value and ax->high_value to be the extremes of the
   array. If you prefer, just set ax->low_value and ax->high_value directly.
*/

extern void maybe_expand_axis_limits(double *farr, int size, axis *ax );
/*
   Like the previous function, only makes sure the axis range doesn't
   shrink below its current value.
*/

extern void compute_axes_details( ongr *on, frame *fr );
/*
   If the `frame' is defined, and both axes have low_value and high_value
   defined, then this works out all the other trifling details of noth
   axis, like what are good round number limits and markings, and where to
   put the axes.
*/

extern void copy_frame();
extern void full_screen_frame( frame *fr );
/* Sets the frame to be the whole of the graphics window */

void sub_frame(
    frame *super_frame,
    int i,
    int j,
    int rows,
    int cols,
    frame *fr
  );
/*
   Obvious and boring, but useful
   *** Beware!!! Stupid. i denotes column number
                         j denotes row number
       (Why studpid? Cos in a different order from rows,cols)
*/

extern void print_and_shrink( frame *fr, char *mess );
/*
    Does 2 things. First it draws the message at the top left of
    frame. Second, it shrinks frame, so that it's top is below
    the message.
*/

extern void draw_axes( frame *fr, ongr *on );
/*
  Draws both axes. Doesn't clear the screen before. We hope you've done that
*/

void draw_hist_axes(frame *fr,ongr *on,string_array *barnames);

extern void plot_in_frame( frame *fr, ongr *on, double *x_arr, double *y_arr, int size, char *mark_code );
/*
    Markcode is a 2-character string.
      First character is line style N => no lines joining points
                                    L => lines joining points

      Second character is dot style N => nothing
                                    D => dots
                                    C => little circles
                                    L => lines down to axis
*/

void simple_plot_linestyle(
    double *x_arr,
    double *y_arr,
    int size,
    char *xlabel,
    char *ylabel,
    char *linestyle
  );
void simple_plot(
    double *x_arr,
    double *y_arr,
    int size,
    char *xlabel,
    char *ylabel
  );

extern void ongr_basic_histogram( frame *fr, double *x, int size );

/* Returns an ag_color associated with code. code can be as large as
   you like but after about 12, it wraps around. Use this if you want, say,
   different bars in a chart to have different colors. You simply use
   color_code_to_ag_color(i) as the color of the i'th bar. */
int color_code_to_ag_color(int color_code);

/* Draws a histogram with N bars where N = dym_rows(freqs)
   Each bar is compartmented into K subbars, on top of each other, where
   K = dym_cols(freqs) and the j'th subcomponent of the i'th bar has
   height = dym_ref(freqs,i,j) and is on top of the j-1'th subcomponent.
   The i'th bar is centered on x=xlo + (i+0.5) * (xhi - xlo) / num_bars */
void ongr_plot_hist(frame *fr,ongr *on,dym *freqs,bool ratio);


#endif /* #ifndef ONGR_H */

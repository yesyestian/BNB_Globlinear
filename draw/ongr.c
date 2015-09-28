
/*
   File:          ongr.c
   Author:        Andrew W. Moore
   Created:       August 19th 1992

   Description:   Some simple graph-drawing
*/

#include <stdio.h>
#include <math.h>
#include "ongr.h" /* Header for simple graph drawing */

#include <string.h>

#define DOT_RADIUS 4.0

void  r2a( frame *f, onpoint *rel_onpoint, onpoint *abs_onpoint )
/*
   Relative coords to abs coords within the frame.
*/
{
  abs_onpoint->x = f->bottom_left.x + rel_onpoint -> x;
  abs_onpoint->y = f->bottom_left.y + rel_onpoint -> y;
}

double frame_width( frame *fr )
{
  return(fr->top_right.x - fr->bottom_left.x);
}

double frame_height( frame *fr )
{
  return(fr->top_right.y - fr->bottom_left.y);
}

void ongr_dot( frame *f, onpoint *p )
{
  onpoint abs;
  double bw = real_min(frame_width(f),frame_height(f))/150.0;
  bw = real_max(DOT_RADIUS * 0.75,real_min(1.5 * DOT_RADIUS,bw));
  r2a(f,p,&abs);
  ag_disc(abs.x,abs.y,bw);
}

void ongr_pixel( frame *f, onpoint *p )
{
  onpoint abs;
  r2a(f,p,&abs);
  ag_pixel(abs.x,abs.y);
}

void ongr_circle( frame *f, onpoint *p, double r )
{
  onpoint abs;
  r2a(f,p,&abs);
  ag_circle(abs.x,abs.y,r);
}

void ongr_disc( frame *f, onpoint *p, double r )
{
  onpoint abs;
  r2a(f,p,&abs);
  ag_disc(abs.x,abs.y,r);
}

void ongr_line( frame *f, onpoint *p1, onpoint *p2 )
{
  onpoint abs1,abs2;
  r2a(f,p1,&abs1);
  r2a(f,p2,&abs2);
  ag_line(abs1.x,abs1.y,abs2.x,abs2.y);
}

void ongr_bar( frame *f, onpoint *p1, onpoint *p2 )
{
  onpoint abs1,abs2;
  r2a(f,p1,&abs1);
  r2a(f,p2,&abs2);
  ag_rectangle(abs1.x,abs1.y,abs2.x,abs2.y);
}

void ongr_cross( frame *f, onpoint *p )
{
  onpoint abs;
  r2a(f,p,&abs);
  ag_line(abs.x-2.0,abs.y,abs.x+2.0,abs.y);  
  ag_line(abs.x,abs.y-2.0,abs.x,abs.y+2.0);  
}

void ongr_print( frame *f, onpoint *p, char *string )
{
  onpoint abs;
  r2a(f,p,&abs);
  ag_print(abs.x,abs.y,string);
}

void clear_axis( axis *ax, bool is_vert )
{
  ax->vertical = is_vert;
  sprintf(ax->label,"%s",(is_vert) ? "Output" : "Input");
  ax -> lowest_val = 0.0;
  ax -> highest_val = 1.0;
}

void clear_ongr( ongr *on )
{
  clear_axis(&on->x_axis,FALSE);
  clear_axis(&on->y_axis,TRUE);
}

void compute_axis_limits( double *farr, int size, axis *ax )
{
  if ( size < 1 )
    printf("compute_axis_limits: no way dude. Axis is empty!\n");
  else
  {
    ax->lowest_val = doubles_min(farr,size);
    ax->highest_val = doubles_max(farr,size);
  }

  if ( Verbosity > 10.0 )
    printf("%s axis limits: (%g , %g)\n",
           (ax->vertical) ? "Vertical" : "Horizontal",
           ax->lowest_val,
           ax->highest_val
          );
}

void maybe_expand_axis_limits( double *farr, int size, axis *ax )
{
  if ( size >= 1 )
  {
    ax->lowest_val = real_min(ax->lowest_val,doubles_min(farr,size));
    ax->highest_val = real_max(ax->highest_val,doubles_max(farr,size));
  }
}

#define CW 8.0
#define CH 18.0
#define CH_IN_NUM 4

/* 17 Oct 96: Mary deleted unreferenced formal parameters "double w, double h" */
void compute_base_start(onpoint *base)
{
  base -> x = (2.0 + CH_IN_NUM) * CW;
  base -> y = 2.5 * CH;
}

/* 17 Oct 96: Mary deleted unreferenced formal parameters "double h" */
/* w and h are the width and height of the frame we'll be using. */
/* Takes into account the space which the labels will take up. */
void compute_hor_axis_markings(axis *ax, double w)
{
  compute_base_start(&ax->start);
  ax->end.x = w - CW;
  ax->end.y = ax->start.y;
  sensible_limits(ax->lowest_val,ax->highest_val,
                  &ax->lo,&ax->hi,&ax->delta
                 );
  if ( ax->end.x - ax->start.x < 150.0 )
    ax->delta = ax->hi - ax -> lo;

  ax->marks = my_irint( (ax->hi - ax->lo) / ax->delta );

  ax->rel_mark_number.x = -CH_IN_NUM/2.0 * CW;
  ax->rel_mark_number.y = -CH;
}

/* 17 Oct 96: Mary removed unreferenced formal parameters "double w" */
void compute_vert_axis_markings(axis *ax, double h)
/* w and h are the width and height of the frame we'll be using. */
/* Takes into account the space which the labels will take up. */
{
  compute_base_start(&ax->start);
  ax->end.x = ax->start.x;
  ax->end.y = h - 2.0 * CH;

  sensible_limits(ax->lowest_val,ax->highest_val,
                  &ax->lo,&ax->hi,&ax->delta
                 );
  if ( ax->end.y - ax->start.y < CH * 3 )
    ax->delta = ax->hi - ax -> lo;

  ax->marks = my_irint( (ax->hi - ax->lo) / ax->delta );

  ax->rel_mark_number.x = -(1 + CH_IN_NUM) * CW;
  ax->rel_mark_number.y = -0.5 * CH;
}

void compute_axes_details(ongr *on, frame *fr)
/*
    Alters horizontal and vertical axis position and labelling
    according to frame details.
*/
{
  double w = frame_width(fr);
  double h = frame_height(fr);
  compute_hor_axis_markings(&on->x_axis,w);
  compute_vert_axis_markings(&on->y_axis,h);
}

void copy_frame( frame *fr1, frame *fr2 )
{
  fr2->bottom_left.x = fr1->bottom_left.x;
  fr2->bottom_left.y = fr1->bottom_left.y;
  fr2->top_right.x = fr1->top_right.x;
  fr2->top_right.y = fr1->top_right.y;
}

void full_screen_frame( frame *fr )
{
  fr->bottom_left.x = 0.0;
  fr->bottom_left.y = 0.0;
  fr->top_right.x = 512.0;
  fr->top_right.y = 512.0;
}

void sub_frame(
    frame *super_frame,
    int i,
    int j,
    int rows,
    int cols,
    frame *fr
  )
{
  double w = frame_width(super_frame) / cols;
  double h = frame_height(super_frame) / rows;
  fr->bottom_left.x = super_frame->bottom_left.x + w * i;
  fr->bottom_left.y = super_frame->bottom_left.y + h * j;
  fr->top_right.x = fr->bottom_left.x + w;
  fr->top_right.y = fr->bottom_left.y + h;
}

void print_and_shrink( frame *fr, char *mess )
{
  onpoint p;
  p.x = CW/2.0;
  p.y = frame_height(fr) - CH * 0.8;
  ongr_print(fr,&p,mess);
  fr -> top_right.y -= CH;
}

void draw_one_axis( frame *fr, axis *ax )
{
  int i;
  onpoint p;
  int color = ag_pen_color();
  ag_set_pen_color(AG_RED);
  ongr_line(fr,&ax->start,&ax->end);

  for ( i = 0 ; i <= ax->marks ; i++ )
  {
    char mark[100];
    double tick;
    p.x = ax->start.x + ((ax->end.x - ax->start.x) * i) / ax->marks;
    p.y = ax->start.y + ((ax->end.y - ax->start.y) * i) / ax->marks;
    ag_set_pen_color(AG_RED);
    ongr_cross(fr,&p);
    p.x += ax->rel_mark_number.x;
    p.y += ax->rel_mark_number.y;
    tick = ax->lo + ax->delta * i;
    sprintf(mark,"%g", fabs(tick) < 1.0e-15 ? 0.0 : tick); /*fix bad rounding*/
    ag_set_pen_color(color);
    ongr_print(fr,&p,mark);
  }
  ag_set_pen_color(color);

/* Finally, we do the label for the axis */

  if ( ax -> vertical )
  {
    p.x = CW/2.0;
    p.y = frame_height(fr) - CH;
  }
  else
  {
    p.x = real_max(0.0,real_min(ax->start.x,frame_width(fr)-CW*strlen(ax->label)));
    p.y = CH/2.0;
  }

  ongr_print(fr,&p,ax->label);
}

void draw_hist_axis( frame *fr, axis *ax , string_array *bar_labels )
{
  int i;
  onpoint p;
  int color = ag_pen_color();
  int num_bars = string_array_size(bar_labels);
  ag_set_pen_color(AG_RED);
  ongr_line(fr,&ax->start,&ax->end);

  for ( i = 0 ; i < num_bars ; i++ )
  {
    char* mark = string_array_ref(bar_labels,i);
    p.x = ax->start.x + ((ax->end.x - ax->start.x) * (i+0.5)) / num_bars;
    p.y = ax->start.y + ((ax->end.y - ax->start.y) * (i+0.5)) / num_bars;
    ag_set_pen_color(AG_RED);
    ongr_cross(fr,&p);
    p.x += ax->rel_mark_number.x;
    p.y += ax->rel_mark_number.y;
    ag_set_pen_color(color);
    if ( mark != NULL ) ongr_print(fr,&p,mark);
  }
  ag_set_pen_color(color);

/* Finally, we do the label for the axis */

  p.x = (frame_width(fr)) / 40.0;
  p.y = CH/2.0;
  ongr_print(fr,&p,ax->label);
}

void draw_axes( frame *fr, ongr *on )
{
  int color = ag_pen_color();
//  ag_set_pen_color(AG_GREEN);
//  ag_rectangle(fr->bottom_left.x,fr->bottom_left.y,
  //             fr->top_right.x-1.0,fr->top_right.y-1.0
    //          );
  ag_set_pen_color(color);
  draw_one_axis(fr,&on->x_axis);
  draw_one_axis(fr,&on->y_axis);
}

void draw_hist_axes(frame *fr,ongr *on,string_array *barnames)
{
  int color = ag_pen_color();
  ag_set_pen_color(AG_GREEN);
  ag_rectangle(fr->bottom_left.x,fr->bottom_left.y,
               fr->top_right.x-1.0,fr->top_right.y-1.0
              );
  ag_set_pen_color(color);
  draw_hist_axis(fr,&on->x_axis,barnames);
  draw_one_axis(fr,&on->y_axis);
}

void plot_graphic_in_frame(
    frame *fr, 
    ongr *on, 
    double *x_arr,
    double *y_arr,
    int size,
    char *mark_code,
    double bar_width,
    int amut_col
  )
/*
    Markcode is a 2-character string.
      First character is line style N => no lines joining onpoints
                                    L => lines joining onpoints

      Second character is dot style N => no dots
                                    D => dots

    If bar_width is negative its ignored. But if positive draws a histogram
    bar at each point, of width "bar_width" (in x axis units).
*/
{
  int i;
  onpoint last;
  axis *xax = &on->x_axis;
  axis *yax = &on->y_axis;
  int color = ag_pen_color();
  ag_set_pen_color(amut_col);

  for ( i = 0 ; i < size ; i++ )
  {
    onpoint this_onpoint;
    this_onpoint.x = xax->start.x + 
             (xax->end.x - xax->start.x) *
             (x_arr[i] - xax->lo) / (xax->hi - xax->lo);
    this_onpoint.y = yax->start.y + 
             (yax->end.y - yax->start.y) *
             (y_arr[i] - yax->lo) / (yax->hi - yax->lo);
    if ( mark_code[1] == 'D' )
    {
      ongr_dot(fr,&this_onpoint);
    }
    else if ( mark_code[1] == 'P' )
    {
      ongr_pixel(fr,&this_onpoint);
    }
    else if ( mark_code[1] == 'C' )
    {
      ongr_disc(fr,&this_onpoint,DOT_RADIUS);
      ongr_circle(fr,&this_onpoint,DOT_RADIUS);
    }
    else if ( mark_code[1] == 'L' )
    {
      onpoint bot;
      bot.x = this_onpoint.x;
      bot.y = on->x_axis.start.y;
      ongr_line(fr,&this_onpoint,&bot);
    }
    
    if ( bar_width > 0.0 )
    {
      onpoint op1,op2;
      double bw = bar_width * (xax->end.x - xax->start.x) / (xax->hi - xax->lo);

      bw *= 0.97;

      op1.x = this_onpoint.x - bw/2.0;
      op1.y = on->x_axis.start.y;
      op2.x = this_onpoint.x + bw/2.0;
      op2.y = this_onpoint.y;
      ongr_bar(fr,&op1,&op2);
    }

    if ( mark_code[0] == 'L' && i > 0 )
    {
      ongr_line(fr,&last,&this_onpoint);
    }

    last.x = this_onpoint.x;
    last.y = this_onpoint.y;
  }
  ag_set_pen_color(color);
}

void plot_in_frame( 
    frame *fr,
    ongr *on,
    double *x_arr,
    double *y_arr,
    int size,
    char *mark_code 
  )
/*
    Markcode is a 2-character string.
      First character is line style N => no lines joining onpoints
                                    L => lines joining onpoints

      Second character is dot style N => no dots
                                    D => dots
                                    C => little circles
                                    L => lines down to axis
				    P => pixels
*/
{
  int col = (eq_string(mark_code,"LN")) ? AG_BLUE :
            (eq_string(mark_code,"LD")) ? AG_PURPLE :
            (eq_string(mark_code,"ND")) ? AG_RED :
            (eq_string(mark_code,"NL")) ? AG_GREEN : AG_BLACK;
            
  plot_graphic_in_frame(fr,on,x_arr,y_arr,size,mark_code,-1.0,col);
}

void simple_plot_linestyle(
    double *x_arr, 
    double *y_arr, 
    int size, 
    char *xlabel, 
    char *ylabel,
    char *linestyle
  )
/*
    linestyle is a 2-character string.
      First character is line style N => no lines joining points
                                    L => lines joining points

      Second character is dot style N => nothing
                                    D => dots
                                    C => little circles
                                    L => lines down to axis
*/
{
  ongr on;
  frame fr;
  clear_ongr(&on);
  compute_axis_limits(x_arr,size,&on.x_axis);
  compute_axis_limits(y_arr,size,&on.y_axis);
  sprintf(on.x_axis.label,xlabel);
  sprintf(on.y_axis.label,ylabel);
  full_screen_frame(&fr);
  compute_axes_details(&on,&fr);
  draw_axes(&fr,&on);
  plot_in_frame(&fr,&on,x_arr,y_arr,size,linestyle);
}

void simple_plot(
    double *x_arr, 
    double *y_arr, 
    int size, 
    char *xlabel, 
    char *ylabel 
  )
{
  simple_plot_linestyle(x_arr,y_arr,size,xlabel,ylabel,"ND");
}

#define MAX_X_ARR_SIZE 300

void very_simple_plot( double *y_arr, int size, char *xlabel, char *ylabel )
{
  double x_arr[MAX_X_ARR_SIZE];
  int i;

  if ( size > MAX_X_ARR_SIZE )
    my_error("ongr.c: MAX_X_ARR_SIZE too small.\n");

  for ( i = 0 ; i < size ; i++ )
    x_arr[i] = (double) i;
  simple_plot(x_arr,y_arr,size,xlabel,ylabel);
}

#define MAX_BARS 500

void ongr_basic_histogram( frame *fr, double *x, int size )
{
  double rough_num_cols = int_min(MAX_BARS/6-1,(int) ceil(sqrt((double) size)));
  double xlo,xhi,delta;
  int bars;
  double freq[MAX_BARS],bases[MAX_BARS];
  int i;
  ongr on;

  sensible_limits(doubles_min(x,size),doubles_max(x,size),&xlo,&xhi,&delta);
/* &&& */
/*
  xlo = 0.0;
  xhi = 10.0;
  delta = 0.2;
*/

  bars = (int) ceil((xhi - xlo) / delta);

  while ( bars < rough_num_cols )
  {
    bars *= 2;
    if ( bars < rough_num_cols )
    {
      bars /= 2;
      bars *= 5;
    }
    if ( bars < rough_num_cols )
      bars *= 2;
  }

  if ( bars > MAX_BARS ) my_error("MAX_BARS too small");

  delta = (xhi - xlo) / bars;
  bars += 3;
  xhi += delta;
  xlo -= delta;

  set_realnums_constant(freq,bars,0.0);
  printf("delta = %g, bars = %d\n",delta,bars);
  
  for ( i = 0 ; i < size ; i++ )
  {
    int bar_num = (int) floor(0.5 + (x[i] - xlo)/delta);
    if ( bar_num < 0 || bar_num >= bars ) my_error("oisdncoisna");
    freq[bar_num] += 1.0;
  }

  for ( i = 0 ; i < bars ; i++ )
    bases[i] = xlo + i * (xhi - xlo) / (bars-1);
  clear_ongr(&on);
  compute_axis_limits(bases,bars,&on.x_axis);
  compute_axis_limits(freq,bars,&on.y_axis);
  sprintf(on.y_axis.label,"Histogram Frequency");
  on.y_axis.lowest_val = 0.0;
  
  compute_axes_details(&on,fr);
  draw_axes(fr,&on);
  plot_graphic_in_frame(fr,&on,bases,freq,bars,"NN",delta,AG_PURPLE);

  if ( size < MAX_X_ARR_SIZE )
  {
    double fake_y[MAX_X_ARR_SIZE];
    double v = 0.5;
    set_realnums_constant(fake_y,size,v);
    plot_in_frame(fr,&on,x,fake_y,size,"ND");
  }  
}

int Amut_cols[] = { AG_BLUE , AG_RED , AG_BLACK , AG_GREEN , AG_PURPLE ,
		  AG_CYAN , AG_GRAY , AG_DARKGREEN , AG_MAGENTA , AG_YELLOW ,
		  AG_DARKRED , AG_OLIVE , AG_DARKBLUE };

/* Returns an ag_color associated with code. code can be as large as
   you like but after about 12, it wraps around. Use this if you want, say,
   different bars in a chart to have different colors. You simply use
   color_code_to_ag_color(i) as the color of the i'th bar. */
int color_code_to_ag_color(int color_code)
{
  int numcols = sizeof(Amut_cols)/sizeof(int);
  int index = color_code % numcols;
  return Amut_cols[index];
}

void ongr_framed_box( frame *f, onpoint *p1, onpoint *p2 , int col)
{
  onpoint abs1,abs2;
  int old_col = ag_pen_color();
  r2a(f,p1,&abs1);
  r2a(f,p2,&abs2);
  ag_set_pen_color(col);
  ag_box(abs1.x,abs1.y,abs2.x,abs2.y);
  ag_set_pen_color(AG_BLACK);
  ag_rectangle(abs1.x,abs1.y,abs2.x,abs2.y);
  ag_set_pen_color(old_col);
}

void graphcoords_to_relcoords(ongr *on,double x,double y,onpoint *p)
{
  axis *xax = &on->x_axis;
  axis *yax = &on->y_axis;

  p ->x = xax->start.x + 
          (xax->end.x - xax->start.x) * (x - xax->lo) / (xax->hi - xax->lo);
  p ->y = yax->start.y + 
          (yax->end.y - yax->start.y) * (y - yax->lo) / (yax->hi - yax->lo);
}

void ongr_colored_bordered_rectangle(frame *fr,ongr *on,double x1,double y1,
				     double x2,double y2,int ag_col)
{
  onpoint p1[1];
  onpoint p2[1];

  graphcoords_to_relcoords(on,x1,y1,p1);
  graphcoords_to_relcoords(on,x2,y2,p2);


  ongr_framed_box(fr,p1,p2,ag_col);
}

double dym_sum_row(dym *x,int row)
{
  dyv *dv = mk_dyv_from_dym_row(x,row);
  double result = dyv_sum(dv);
  free_dyv(dv);
  return result;
}

/* Draws a histogram with N bars where N = dym_rows(freqs)
   Each bar is compartmented into K subbars, on top of each other, where
   K = dym_cols(freqs) and the j'th subcomponent of the i'th bar has
   height = dym_ref(freqs,i,j) and is on top of the j-1'th subcomponent. 
   The i'th bar is centered on x=xlo + (i+0.5) * (xhi - xlo) / num_bars */
void ongr_plot_hist(frame *fr,ongr *on,dym *freqs,bool ratio)
{
  int i;
  int num_bars = dym_rows(freqs);
  int num_classes = dym_cols(freqs);
  double xlo = on -> x_axis.lo;
  double xhi = on -> x_axis.hi;
  double bar_width = 0.9 * (xhi - xlo) / num_bars;

  for ( i = 0 ; i < num_bars ; i++ )
  {
    double xmid = xlo + (i + 0.5) * (xhi - xlo) / num_bars;
    double x1 = xmid - bar_width/2;
    double x2 = xmid + bar_width/2;
    int j;
    double sumy = 0.0;
    double total_sum_y = real_max(1e-5,dym_sum_row(freqs,i));

    for ( j = 0 ; j < num_classes ; j++ )
    {
      double y1 = sumy;
      double dy = dym_ref(freqs,i,j) / ((ratio) ? (total_sum_y/100.0) : 1.0);
      double y2 = sumy + dy;
      int amut_col = color_code_to_ag_color(j);
      ongr_colored_bordered_rectangle(fr,on,x1,y1,x2,y2,amut_col);
      sumy = y2;
    }
  }
}

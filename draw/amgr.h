/*
   File:            amgr.h
   Author:          Andrew Moore
   Created:         7th Feb 1990

   Description:     Standard andrew header file for graphics
*/

#ifndef AMGR_H
#define AMGR_H

#include "./utils/amut.h"
#include "./draw/amxw.h"

/*
   Example of use:

    All coordinates are in doubles (x,y) with origin at bottom left.
    The screen covers (0,0) to (512,512)

  {
    double x,y;
    int button;

    ag_on("plop.ps"); Will open a 512x512 pixel x-window, and also save
                      results to a postscript file called plop.ps.
                      If you don't want anything saved to a file, call
                      ag_on("");

    ag_circle(256.0,256.0,128.0);  Draws a circle centered on the center
                                   of the window, radius quarter window width.

    ag_print(10.0,10.0,"hi"); printed at bottom left

    button = ag_get_xy(&x,&y);  Waits for the user to click button on window.

    printf("You clicked %s button\n",(button==1)?"left":(button==2)?"middle":
                                     "right"
          );
    ag_dot(x,y);  Draws dot where user clicked button

    ag_off();     The drawing is now saved in the postscript file
  }
*/

#ifdef PC_MVIS_PLATFORM
extern double CharHeight, CharWidth;
#define AG_CHAR_HEIGHT CharHeight
#define AG_CHAR_WIDTH  CharWidth
#else
#define AG_CHAR_HEIGHT 16.0
#define AG_CHAR_WIDTH  7.0
#endif

typedef struct agbox
{
  double xlo;
  double ylo;
  double xhi;
  double yhi;
} agbox;

void ag_on(char *fname);

/*Under Windows, this function will update the graphics window.
  Under other operating systems, it will do nothing.
  This is needed because under Windows, any graphics commands that are
  issued will not be seen until a call to ag_off() or ag_show().*/
void ag_show(void);

void ag_dot(double x, double y);

/* Draw a tiny dot (one pixel in size on ag_ window, 1/512th of width
   of entire drawing area in postscript) */
void ag_pixel(double x, double y);

void ag_line(double u, double v, double x, double y);

void ag_circle(double x, double y, double r);

void ag_box(double x, double y, double u, double v);

void ag_disc(double x, double y, double r);

/* This function is like ag_print except it ignores the current
   scale factor. Thus if you call with x = 10.0, y = 490.0 you are
   guaranteed to draw in the top left of the window no matter what
   scale factor is currently in use. */
void ag_print_absolute_coords(double x, double y, char *s);

void ag_print(double x, double y, char *s);

void ag_pen(double c);

/* Declared but not defined - sir /9/6/2000
void ag_boxed_print(double x, double y, double u, double v, char *s); */

/* Waits for the user to click in the graphics window.
   Stores the coordinates of where the user clicked in
   (*xptr,*yptr).

    Returns 1 if the user clicked the left button
    Returns 3 if right button
    Returns 2 if middle button, but I guess that doesn't
       happen under NT */
int ag_get_xy(double *x_ptr, double *y_ptr);

/*Does the same as ag_get_xy(), but prompts the user*/
int ag_get_xy_prompted(char *prompt,double *x_ptr, double *y_ptr);

void ag_off(void);

void ag_cross(double x, double y);

/* Declared but not defined - sir /9/6/2000
void ag_polystart(void);
void ag_polypoint(double x, double y);
void ag_polyfill(void);
void ag_polyedge(void); */

void ag_arrow(double x, double y, double u, double v);

void ag_rectangle(double x, double y, double u, double v);

void ag_clipped_dot(double x, double y);

void ag_clipped_line(double p, double q, double u, double v);

/* Declared but not defined - sir /9/6/2000
void ag_char_point(double x, double y, char c);
void ag_char_square(double x, double y, char c); */

void ag_uparrow(double x, double y);

void ag_downarrow(double x, double y);

void ag_rounded_box(double x, double y, double u, double v, double r);

void click_on_message(double xlo, double ylo, double xhi, double yhi, char *m);

void wait_for_click(void);

/* Declared but not defined - sir /9/6/2000
void ag_lock(void);
void ag_unlock(void); */

/* By default, graphics, postscript and apicts are created using a
   coordinate system in which bottom left is at (0,0) and top
   right is (512,512). This allows you to change the bottom left
   and top right. */
void set_ag_frame(double xlo,double ylo,double xhi,double yhi);

/* Change the physical size of the window on the screen. The coordinate
   system stretches so that something that would previously have appeared
   in the top right corner still appears in the top right corner. (and
   ditto for bottom left) */
void ag_window_shape(int width_in_pixels,int height_in_pixels);


/* Causes the ag_ window to disappear. It can still reapper if
   ag_on("") is called again */
void ag_disappear();

void ag_circle_transformed_radius(double x, double y, double r);

/* Foreground and background colors */
extern int AG_DISPLAY;
extern int AG_BKGD;
extern int AG_HIGHLIGHT;
extern int AG_HYPERLINK;

/* declarations added by Mary on 2 July 95; part of amxw color code */

void ag_set_pen_color(int amut_col);
int ag_pen_color();
int ag_spectrum_color(double fract);
void ag_rgb(int ag_color, double *r_r, double *r_g, double *r_b);

void ag_print_vertical(double x,double y,char *string);

/* Draws a line segment in the ag graphics window consisting of
   that set of points { (x,y) | ( (x,y) - (px,py) ) . (dx,dy) == 0 }.

   This is a line perpendicular to (dx,dy) that travels through the point
   (px,py).

   The function is clever about clipping.
*/
void ag_cut(double px,double py,double dx,double dy);

/*Now we define all functions that deal with bounding boxes, for
  our purposes called an ag_box.*/
agbox *mk_agbox(double xlo,double ylo,double xhi,double yhi);
void free_agbox(agbox *agb);
agbox *mk_copy_agbox(agbox *agb);
void fprintf_agbox(FILE *s,char *m1,agbox *x,char *m2);
double break_agval(double lo,double hi,double frac);
double vertical_break_agval(agbox *agb,double frac);

void agbox_piebar(agbox *agb,dyv *widths
                  ,string_array *left,string_array *right
                  ,ivec *colors);

void ag_cross(double x, double y);
void ag_smallbox(double x, double y);
void ag_smalltriangle(double x, double y);
void ag_smallcircle(double x, double y);

#define DOT_MARKTYPE 0
#define CROSS_MARKTYPE 1
#define SMALLBOX_MARKTYPE 2
#define SMALLTRIANGLE_MARKTYPE 3
#define SMALLCIRCLE_MARKTYPE 4
#define PIXEL_MARKTYPE 5
#define FANCY_MARKTYPE 6
#define BIG_CROSS_MARKTYPE 7
#define NUM_MARKTYPES 8

void ag_mark(double x,double y,int marktype);

/*
  300secs          400secs       500secs                      700secs
   +------------------------------------------------------------+
   |                |              |                            |
   |                |              |                            |
   +------------------------------------------------------------+

   Suppose parent was the large box above. Suppose you were told that parent
   covers the time period 300secs--700secs. And suppose you wanted to
   draw in a subregion corresponding to 400-500 secs. What should your agbox
   be? The answer is returned by this function.
*/
agbox *mk_agbox_corresponding_to_time(agbox *parent,
				      double parent_start_time,
				      double parent_end_time,
				      double local_start_time,
				      double local_end_time);

double ag_value_corresponding_to_time(agbox *parent,
				      double parent_start_time,
				      double parent_end_time,
				      double local_time);

double get_time_corresponding_to_agbox(agbox *agb,double start,double end,double x);

void agbox_render_string_at_bottom(agbox *agb,char *string);

/* Prints horizontally at top left of box */
void agbox_render_string_horizontal(agbox *agb,char *string);

/* Prints vertically at top left of box */
void agbox_render_string_vertical(agbox *agb,char *string);

/* Prints horizontally at top left of box, unless it can't fit the string
   in horizontally, and it can fit more letters in vertically, in which
   case prints vertically.  It does clip letters.*/
void agbox_render_string(agbox *agb,char *string);

bool enough_room_for_string_in_agbox(agbox *agb,char *string);

/* As above, except only renders if enough room */
void agbox_render_string_if_enough_room(agbox *agb,char *string);

/* As above, except resizes agb to the leftover horizontal space.
  Returns TRUE iff there was enoguh room.*/
bool agbox_render_string_if_enough_room_and_resize(agbox **agb,char *string);

/*Fits each string into the box the best way it knows how.  Returns how many rows/cols it
  needed to fit them all, or 0 if it couldn't do it.*/
int agbox_render_string_array(agbox *agb,string_array *sa);
int render_string_array_horizontal_if_possible(agbox *agb,string_array *sa);
int render_string_array_vertical_if_possible(agbox *agb,string_array *sa);

/* Uses full name if it fits, else uses short name
   if it fits. Shortname may be NULL in which case it's never
   drawn.  */
void agbox_render_string_or_abbrev(agbox *agb,char *name,char *shortname);

agbox *mk_horizontal_stripe_absolute(agbox *parent,double aglo,double aghi);
agbox *mk_horizontal_stripe_relative(agbox *parent,
				     double fraclo,double frachi);

/* Imagine dividing the first arg (parent) up into "numstrips" horizontal
   strips of equal height, exactly covering and partitioning parent.
   What would be the agbox corresponding to the stripnum'th of these stripes
   from the top (where stripnum==0 corresponds to the topmost stripe?) */
agbox *mk_agbox_horizontal_substripe(agbox *parent,int stripnum,int numstrips);

agbox *mk_vertical_stripe_absolute(agbox *parent,double aglo,double aghi);
agbox *mk_vertical_stripe_relative(agbox *parent,double fraclo,double frachi);
agbox *mk_agbox_vertical_substripe(agbox *parent,int stripnum,int numstrips);

/* Draws a line segment in absolute coordinates.
   It DOES clip the line so that only its components INSIDE agb are
   drawn. */
void render_line_segment_absolute(agbox *agb,
				  double x0,double y0,double x1,double y1);

/* Draws a relative-coordinate line segment. Assume that the
   bottom left of the agbox should be identified with the
   relative coordinate (0,0) and the top right should be identified
   with (1,1).

   It DOES clip the line so that only its components INSIDE agb are
   drawn. */
void render_line_segment_relative(agbox *agb,
				  double x0,double y0,double x1,double y1);

/* Draws a virtual-coordinate line segment. Assume that the
   bottom left of the agbox should be identified with the
   virtual-coordinates-point (vx0,vy0) and the top right with
   (vx1,vy1). Then this draws the line segment that begins with
   (x0,y0) in virtual coordinates and ends with (x1,y1).

   It DOES clip the line so that only its components INSIDE agb are
   drawn. */
void render_virtual_line_segment(agbox *agb,
				 double vx0,double vy0,double vx1,double vy1,
				 double x0,double y0,double x1,double y1);

double agbox_width(agbox *agb);
double agbox_height(agbox *agb);
double agbox_thinnest_dimension(agbox *agb);

void agbox_border(agbox *agb);
void agbox_color_border(agbox *agb,int agcolor);
void agbox_fill(agbox *agb);
void agbox_color_fill(agbox *agb,int agcolor);
/*Returns TRUE iff (x,y) is within agb (or on the boundary)*/
bool is_in_agbox(agbox *agb,double x,double y);

/* If you are working in an application that is driven by dispatch_xxx( )
   type functions, it is often worth using special_ag_on instead of ag_on.
   This is because when the comment on sessionname is invoked by the
   user, all graphs will automatically be saved in auto-named .ps files
   and linked within a latex document into a session log. It is harmless
   to use this instead of ag_on if you're not sure. */
void special_ag_on(char *filename);

void ag_cross(double x, double y);
void ag_smallbox(double x, double y);
void ag_smalltriangle(double x, double y);
void ag_smallcircle(double x, double y);

void ag_mark(double x,double y,int marktype);

/* GJG: added the ability to choose how ag_colors show up in saved PostScript
   files.  Default is as before, different ag_colors show up as actual colors
   in the PostScript file.  But, on a grayscale device, it's hard to tell the
   difference between different colors.  So, two options: colors get mapped to
   different gray levels (this option was there already, but I changed the
   mapping so that colors with adjacent indices were no longer similar), or
   colors get mapped to different patterns (e.g., solid or crosshatched or with
   diagonal lines). The following three functions switch between these
   behaviors in the obvious way.

   IMPORTANT: the call to ps_usepats must come before the call to ag_on,
   since ps_prologue decides whether to emit the pattern definitions based
   on the value of PS_MONO. */

void ps_usegray();
void ps_usepats();
void ps_usecolor();


#endif /* #ifndef AMGR_H */

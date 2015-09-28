
/*
   File:            amgr.c
   Author:          Andrew Moore
   Created:         9th Feb 1990
   Modified:        4 Feb 96

   Description:     Andrews Graphics
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "amgr.h"
#include "./utils/command.h"

/* The default foreground, background, highlight and hyperlink colors. */
int AG_DISPLAY = AG_BLUE;
int AG_HYPERLINK = AG_DARKGREEN;
int AG_HIGHLIGHT = AG_RED;
int AG_BKGD = AG_WHITE;


int PS_MONO = 0;  /* GJG: flag to make PostScript use
		     patterns instead of colors */
#define NUM_PSPATS 6
char *pspatnames[NUM_PSPATS] = {"XPat","LSlashPat","RSlashPat","PlusPat","HPat","VPat"};

/* Functions from either amxw.c, apict.c, or amsv.c (Sunview) */

extern void screen_show(void);
extern void screen_update(void);
extern void screen_clear(void);
extern void screen_dot(double x, double y);
void screen_pixel(double x, double y);
extern void screen_line(double x1, double y1, double x2, double y2);
extern void screen_box(double x, double y, double w, double h);
extern void screen_prbox(double u, double v, double x, double y, char *s);
extern void screen_circle(double x, double y, double r);
extern void screen_disc(double x, double y, double r);
extern void screen_print(char *s, double x, double y);
extern void screen_on(void);
extern void screen_off(void);
extern void screen_lock(void);
extern void screen_unlock(void);
void screen_window_shape(int width_in_pixels,int height_in_pixels);
void screen_disappear();

extern void apict_show(void);
extern void apict_clear(void);
extern void apict_dot(double x, double y);
extern void apict_pixel(double x, double y);
extern void apict_line(double x1, double y1, double x2, double y2);
extern void apict_box(double x, double y, double w, double h);
extern void apict_prbox(double u, double v, double x, double y, char *s);
extern void apict_circle(double x, double y, double r);
extern void apict_disc(double x, double y, double r);
extern void apict_print(char *s, double x, double y);
extern void apict_on(void);
extern void apict_off(void);
extern bool apict_on_p();

extern void simple_set_pen_color(int color);
extern int simple_pen_color();
extern int simple_spectrum_color(double fract);
extern int simple_rgb_color(int rgb);

extern void simple_read_keyboard(void);
extern int simple_get_xy(char *prompt, double *rx, double *ry);
extern int simple_get_release_xy(char *prompt, double *rx, double *ry);
extern void simple_set_my_max_xy(double mx, double my);

/************ GRAPHICS.C FOLLOWS *****************/
/* Postscript Graphics Functions */

bool Ag_active = FALSE;
bool Ps_active = FALSE;

FILE *Ps_stream = NULL;

void ps_error(char *s)
{
  if ( Ps_stream != NULL )
    if ( fclose(Ps_stream) )
      printf("\n*** Warning from ps_error: Problems closing Ps_stream\n");
  printf("Postscript Generation error:\n%s\n",s);
  exit(1);
}

void ps_prologue(FILE *stream,double maxx,double maxy)
{
  fprintf(stream,"%%!PS-Adobe-but-minimally-conforming\n");
  fprintf(stream,"%%%%DocumentFonts: Times-Roman\n");
  fprintf(stream,"%% MTPGD-generated postscript (by Andrew W. Moore)\n");
  fprintf(stream,"/dot { newpath 2.2 0 360 arc fill } def\n");
  fprintf(stream,"/circle { newpath 0 360 arc stroke } def\n");
  fprintf(stream,"/greycirc { newpath 0 360 arc\n");
  fprintf(stream,"         gsave 0.5 setgray fill grestore } def\n");
  fprintf(stream,"/disc { newpath 0 360 arc\n");
  fprintf(stream,"         fill } def\n");
  fprintf(stream,"/line { newpath moveto lineto stroke } def\n");
  fprintf(stream,"/box { newpath moveto 1 index 0 rlineto\n");
  fprintf(stream,"               0 1 index rlineto\n");
  fprintf(stream,"               1 index neg 0 rlineto\n");
  fprintf(stream,"               closepath\n");
  fprintf(stream,"               fill pop pop\n");
  fprintf(stream,"     } def\n");
  fprintf(stream,"/px { newpath moveto 1 0 rlineto\n");
  fprintf(stream,"               0 1 rlineto\n");
  fprintf(stream,"               -1 0 rlineto\n");
  fprintf(stream,"               closepath\n");
  fprintf(stream,"               fill\n");
  fprintf(stream,"     } def\n");
  fprintf(stream,"/amprint { /Courier findfont 13 scalefont setfont\n");
  fprintf(stream,"         newpath moveto show } def\n");
  fprintf(stream,"/point { newpath 2 copy 5 0 360 arc\n");
  fprintf(stream,"         gsave 1 setgray fill grestore\n");
  fprintf(stream,"         gsave 0.1 setlinewidth stroke grestore newpath\n");
  fprintf(stream,"         /Helvetica findfont 5 scalefont setfont\n");
  fprintf(stream,"         moveto\n");
  fprintf(stream,"         dup stringwidth pop 2 div neg -1.77 rmoveto\n");
  fprintf(stream,"         show\n");
  fprintf(stream,"       } def\n");
  fprintf(stream,"/q { setgray newpath moveto dup 0 rlineto\n");
  fprintf(stream,"     0 1 index rlineto neg 0 rlineto\n");
  fprintf(stream,"     closepath fill\n");
  fprintf(stream,"   } def\n");
  fprintf(stream,"/prbox { /prboxh exch def /prboxw exch def \n");
  fprintf(stream,"         /prboxy exch def /prboxx exch def\n");
  fprintf(stream,"         newpath /Courier findfont prboxh 0.75 mul \n");
  fprintf(stream,"         scalefont setfont dup stringwidth pop prboxw \n");
  fprintf(stream,"         exch sub 2 div prboxx add prboxh 4 div prboxy\n");
  fprintf(stream,"         add moveto show } def \n");
  fprintf(stream,"/bounder { 0 setgray\n");

//#define NORMAL

#ifdef NORMAL
  fprintf(stream,"           newpath -1 -1 moveto 512 -1 lineto 512 512 lineto\n");
  fprintf(stream,"                                -1 512 lineto closepath stroke\n");
#else
  fprintf(stream,"           newpath -1 -1 moveto %g -1 lineto %g %g lineto\n"
                 "                   -1 %g lineto closepath stroke\n",
	  maxx,maxx,maxy,maxy);
#endif
  fprintf(stream,"         } def\n");
  fprintf(stream,"0.5 setlinewidth 1 setlinecap 1 setlinejoin\n");
#ifdef NORMAL
  fprintf(stream,"%%%%BoundingBox: 30 180 542 692\n");
#else
  fprintf(stream,"%%%%BoundingBox: 30 180 %g %g\n",30+maxx,180+maxy);
#endif
  fprintf(stream,"%%%%EndProlog\n");
  fprintf(stream,"%%%%Page: 1 1\n");
  fprintf(stream,"30 180 translate 1 1 scale\n");
  fprintf(stream,"gsave\n");
  fprintf(stream,"newpath \n");
#ifdef NORMAL
  fprintf(stream,"0 0 moveto\n");
  fprintf(stream,"512 0 lineto\n");
  fprintf(stream,"512 512 lineto\n");
  fprintf(stream,"0 512 lineto\n");
#else
  fprintf(stream,"0 0 moveto\n");
  fprintf(stream,"%g 0 lineto\n",maxx);
  fprintf(stream,"%g %g lineto\n",maxx,maxy);
  fprintf(stream,"0 %g lineto\n",maxy);
#endif
  fprintf(stream,"closepath\n");
  fprintf(stream,"clip\n");
  fprintf(stream,"newpath\n");

  /* add pattern dictionaries -- GJG */
  if (2 == PS_MONO) {
    /* crosshatch */
    fprintf(stream, "<< /PatternType 1\n");
    fprintf(stream, "   /PaintType 1\n");
    fprintf(stream, "   /TilingType 1\n");
    fprintf(stream, "   /BBox [0 0 7 7]\n");
    fprintf(stream, "   /XStep 7 /YStep 7\n");
    fprintf(stream, "   /PaintProc { begin\n");
    fprintf(stream, "       1 1 1 setrgbcolor\n");
    fprintf(stream, "       .7 setgray 0 0 7 7 rectfill 0 setgray\n");
    fprintf(stream, "       0 0 moveto 7 7 lineto stroke\n");
    fprintf(stream, "       7 0 moveto 0 7 lineto stroke\n");
    fprintf(stream, "     end } bind\n");
    fprintf(stream, ">>\n");
    fprintf(stream, "matrix makepattern /XPat exch def\n");

    /* left diagonal lines */
    fprintf(stream, "<< /PatternType 1\n");
    fprintf(stream, "   /PaintType 1\n");
    fprintf(stream, "   /TilingType 1\n");
    fprintf(stream, "   /BBox [0 0 7 7]\n");
    fprintf(stream, "   /XStep 7 /YStep 7\n");
    fprintf(stream, "   /PaintProc { begin\n");
    fprintf(stream, "       1 1 1 setrgbcolor\n");
    fprintf(stream, "       .7 setgray 0 0 7 7 rectfill 0 setgray\n");
    fprintf(stream, "       0 0 moveto 7 7 lineto stroke\n");
    fprintf(stream, "     end } bind\n");
    fprintf(stream, ">>\n");
    fprintf(stream, "matrix makepattern /LSlashPat exch def\n");

    /* right diagonal lines */
    fprintf(stream, "<< /PatternType 1\n");
    fprintf(stream, "   /PaintType 1\n");
    fprintf(stream, "   /TilingType 1\n");
    fprintf(stream, "   /BBox [0 0 7 7]\n");
    fprintf(stream, "   /XStep 7 /YStep 7\n");
    fprintf(stream, "   /PaintProc { begin\n");
    fprintf(stream, "       1 1 1 setrgbcolor\n");
    fprintf(stream, "       .7 setgray 0 0 7 7 rectfill 0 setgray\n");
    fprintf(stream, "       0 7 moveto 7 0 lineto stroke\n");
    fprintf(stream, "     end } bind\n");
    fprintf(stream, ">>\n");
    fprintf(stream, "matrix makepattern /RSlashPat exch def\n");

    /* plusses */
    fprintf(stream, "<< /PatternType 1\n");
    fprintf(stream, "   /PaintType 1\n");
    fprintf(stream, "   /TilingType 1\n");
    fprintf(stream, "   /BBox [0 0 7 7]\n");
    fprintf(stream, "   /XStep 7 /YStep 7\n");
    fprintf(stream, "   /PaintProc { begin\n");
    fprintf(stream, "       1 1 1 setrgbcolor\n");
    fprintf(stream, "       .7 setgray 0 0 7 7 rectfill 0 setgray\n");
    fprintf(stream, "       0 4 moveto 7 4 lineto stroke\n");
    fprintf(stream, "       4 0 moveto 4 7 lineto stroke\n");
    fprintf(stream, "     end } bind\n");
    fprintf(stream, ">>\n");
    fprintf(stream, "matrix makepattern /PlusPat exch def\n");

    /* horiz lines */
    fprintf(stream, "<< /PatternType 1\n");
    fprintf(stream, "   /PaintType 1\n");
    fprintf(stream, "   /TilingType 1\n");
    fprintf(stream, "   /BBox [0 0 5 5]\n");
    fprintf(stream, "   /XStep 5 /YStep 5\n");
    fprintf(stream, "   /PaintProc { begin\n");
    fprintf(stream, "       1 1 1 setrgbcolor\n");
    fprintf(stream, "       .7 setgray 0 0 5 5 rectfill 0 setgray\n");
    fprintf(stream, "       0 2 moveto 5 2 lineto stroke\n");
    fprintf(stream, "     end } bind\n");
    fprintf(stream, ">>\n");
    fprintf(stream, "matrix makepattern /HPat exch def\n");

    /* vert lines */
    fprintf(stream, "<< /PatternType 1\n");
    fprintf(stream, "   /PaintType 1\n");
    fprintf(stream, "   /TilingType 1\n");
    fprintf(stream, "   /BBox [0 0 5 5]\n");
    fprintf(stream, "   /XStep 5 /YStep 5\n");
    fprintf(stream, "   /PaintProc { begin\n");
    fprintf(stream, "       1 1 1 setrgbcolor\n");
    fprintf(stream, "       .7 setgray 0 0 5 5 rectfill 0 setgray\n");
    fprintf(stream, "       2 0 moveto 2 5 lineto stroke\n");
    fprintf(stream, "     end } bind\n");
    fprintf(stream, ">>\n");
    fprintf(stream, "matrix makepattern /VPat exch def\n");

    /* solid */
    fprintf(stream, "<< /PatternType 1\n");
    fprintf(stream, "   /PaintType 1\n");
    fprintf(stream, "   /TilingType 1\n");
    fprintf(stream, "   /BBox [0 0 5 5]\n");
    fprintf(stream, "   /XStep 5 /YStep 5\n");
    fprintf(stream, "   /PaintProc { begin\n");
    fprintf(stream, "       1 1 1 setrgbcolor\n");
    fprintf(stream, "       0 setgray 0 0 5 5 rectfill\n");
    fprintf(stream, "     end } bind\n");
    fprintf(stream, ">>\n");
    fprintf(stream, "matrix makepattern /SolidPat exch def\n");
  }
}

void ps_usegray()
{
  PS_MONO = 1;
}

void ps_usepats()
{
  PS_MONO = 2;
}

void ps_usecolor()
{
  PS_MONO = 0;
}

void ps_on(char *filename,double maxx,double maxy)
{
  if ( Ps_stream != NULL )
    ps_error("Called ps_on() when ps_recording already active\n");

  Ps_stream = fopen(filename,"w");

  if ( Ps_stream == NULL )
  {
    printf("Cannot open file %s for writing\n",filename);
    ps_error("File Problems");
  }

  ps_prologue(Ps_stream,maxx,maxy);
}

void ps_off(void)
{
  fprintf(Ps_stream,"grestore\n");
#ifdef NORMAL
  fprintf(Ps_stream,"bounder\n");
#endif
  fprintf(Ps_stream,"showpage\n%%%%Trailer\n");
  if ( fclose(Ps_stream) )
    ps_error("ps_off: problems closing Ps_stream");

  Ps_stream = NULL;
}

void ps_dot(double x, double y)
{
  fprintf(Ps_stream,"%g %g dot\n",x,y);
}

void ps_pixel(double x, double y)
{
  fprintf(Ps_stream,"%g %g px\n",x,y);
}

void ps_circle(double x, double y, double r)
{
  fprintf(Ps_stream,"%g %g %g circle\n",x,y,r);
}

void ps_disc(double x, double y, double r)
{
  fprintf(Ps_stream,"%g %g %g disc\n",x,y,r);
}

void ps_grey_circle(double x, double y, double r)
{
  fprintf(Ps_stream,"%g %g %g greycirc\n",x,y,r);
}

void ps_line(double u, double v, double x, double y)
{
  fprintf(Ps_stream,"%g %g %g %g line\n",u,v,x,y);
}

void ps_point(double x, double y, int n)
{
  fprintf(Ps_stream,"(%d) %g %g point\n",n,x,y);
}

double Color_now = 0.0;

void ps_color(int amut_col)
{
  /* fprintf(Ps_stream,"%g setgray\n",((double) amut_col)/(NUM_COLORS-1)); */
  if (amut_col == AG_BLACK)
    fprintf(Ps_stream, "0 setgray\n");
  else if (amut_col == AG_WHITE)
    fprintf(Ps_stream, "1 setgray\n");
  else
    fprintf(Ps_stream,"%g setgray\n",
	    ((double) (amut_col * 17 % NUM_COLORS))/(NUM_COLORS+1));
  /* GJG: changed so that adjacent colors would contrast more, and no color
     is too close to white */
}

void ps_rgb(double r,double g,double b)
{
  fprintf(Ps_stream,"%g %g %g setrgbcolor\n",r,g,b);
}

void ps_print(double x, double y, char *s)
{
  fprintf(Ps_stream,"(%s) %G %g amprint\n",s,x,y);
}

void ps_filled_rectangle(double x, double y, double u, double v)
{
  fprintf(Ps_stream,"%g %g %g %g box\n",x-u,y-v,u,v);
}

void ps_prbox(double u, double v, double x, double y, char *s)
{
  double t;
  if ( x < u )
  {
    t = x;
    x = u;
    u = t;
  }
  if ( y < v )
  {
    t = y;
    y = v;
    v = t;
  }
  fprintf(Ps_stream,"(%s) %g %g %g %g prbox\n",s,u,v,x-u,y-v);
}

#define MAX_POINTS 200

static double ps_x_points[MAX_POINTS];
static double ps_y_points[MAX_POINTS];
int ps_points_index = 0;

void ps_polystart(void)
{
  ps_points_index = 0;
}

void ps_polypoint(double x, double y)
{
  if ( ps_points_index >= MAX_POINTS )
    printf("graphics.c: POLYGON TOO LARGE: More than %d points!\n",MAX_POINTS);
  else
  {
    ps_x_points[ps_points_index] = x;
    ps_y_points[ps_points_index] = y;
    ps_points_index += 1;
  }
}

void ps_polyfill(void)
{
  if ( ps_points_index > 2 )
  {
    int i;
    fprintf(Ps_stream,"newpath %g %g moveto\n",ps_x_points[0],ps_y_points[0]);
    for ( i = 1 ; i < ps_points_index ; i++ )
    {
      fprintf(Ps_stream,"%g %g lineto\n",ps_x_points[i],ps_y_points[i]);
    }
    fprintf(Ps_stream,"closepath fill\n");
  }
}

void ps_polyedge(void)
{
  if ( ps_points_index > 2 )
  {
    int i;
    fprintf(Ps_stream,"newpath %g %g moveto\n",ps_x_points[0],ps_y_points[0]);
    for ( i = 1 ; i < ps_points_index ; i++ )
    {
      fprintf(Ps_stream,"%g %g lineto\n",ps_x_points[i],ps_y_points[i]);
    }
    fprintf(Ps_stream,"closepath stroke\n");
  }
}

/************* NOW THE AG STUFF *****************/

double Frame_x_lo = 0.0;
double Frame_y_lo = 0.0;
double Frame_x_hi = 512.0;
double Frame_y_hi = 512.0;

void ag_transform(double *x,double *y)
{
  *x = 512.0 * (*x - Frame_x_lo) / (Frame_x_hi - Frame_x_lo);
  *y = 512.0 * (*y - Frame_y_lo) / (Frame_y_hi - Frame_y_lo);
}

void ag_untransform(double *x,double *y)
{
  *x = Frame_x_lo + (Frame_x_hi - Frame_x_lo) * *x / 512.0;
  *y = Frame_y_lo + (Frame_y_hi - Frame_y_lo) * *y / 512.0;
}

/* By default, graphics, postscript and apicts are created using a
   coordinate system in which bottom left is at (0,0) and top
   right is (512,512). This allows you to change the bottom left
   and top right. */
void set_ag_frame(double xlo,double ylo,double xhi,double yhi)
{
  Frame_x_lo = xlo;
  Frame_y_lo = ylo;
  Frame_x_hi = xhi;
  Frame_y_hi = yhi;
}

void ag_window_shape(int width_in_pixels,int height_in_pixels)
{
  screen_window_shape(width_in_pixels,height_in_pixels);
}

/* Causes the ag_ window to disappear. It can still reapper if
   ag_on("") is called again */
void ag_disappear()
{
  ag_off();
  screen_disappear();
}

/* Pictures/graphs can be built up in apict structures or, on the Unix
   side, may be displayed directly in a window on the screen.

   To start building a picture in an apict structure first call
   apict_on.  To start drawing a picture directly onto the screen --
   currently only implemented on the Unix side -- first call ag_on....
   Calls to ag_line, ag_circle, etc., will now add lines, circles,
   etc., to the picture being created.  When the picture is finished
   call ag_off/apict_off.  The latter returns the apict containing the
   picture....

   In the Unix world, it is fine to call both ag_on and apict_on, so
   that the picture both appears on the screen and is also built up in
   an apict.
*/

bool ag_on_p = FALSE; /* Set to TRUE when ag_on is called.  Set to
                         FALSE when ag_off is called. */

void ag_on_with_non_standard_maxx_maxy(char *fname,double maxx,double maxy)
{
  extern FILE *Ps_stream;

  ag_on_p = TRUE;
  screen_on();
  screen_show();

  if ( !Ps_active ) Ps_active = ( Ps_stream != NULL );

  if ( Ps_active )
  {
    printf("*** ag_on: You haven't turned previous graphic off.\n");
    printf("*** I will do so for you.\n");
    ps_off();
    Ps_active = FALSE;
  }

  if ( fname != NULL && !eq_string(fname,"*") && !eq_string(fname,"") )
  {
    FILE *s = fopen(fname,"w");
    if ( s != NULL )
    {
      fclose(s);
      ps_on(fname,maxx,maxy);
      Ps_active = TRUE;
    }
  }
  Ag_active = TRUE;
}

void ag_on(char *fname)
{
  ag_on_with_non_standard_maxx_maxy(fname,512.0,512.0);
}

void ag_off(void)
{
  screen_show();
  if ( Ps_active ) ps_off();
  Ag_active = FALSE;
  Ps_active = FALSE;
  screen_off();
  ag_on_p = FALSE;
}

void ag_show(void)
{
#ifdef PC_MVIS_PLATFORM
  screen_update();
#endif
}

/* ==================================================================== */
/* The following functions determine whether ag_on and/or apict_on
   have been called more recently than ag_off and/or apict_off.  If so
   they add the picture element to the screen and/or to the apict being
   built up.  */
/* ==================================================================== */

void ag_dot(double x, double y)
{
  ag_transform(&x,&y);
  if (ag_on_p)
  {
    screen_dot(x,y);
    screen_show();
  }
  if (apict_on_p())
  {
    apict_dot(x,y);
  }
  if ( Ps_active ) ps_dot(x,y);
}

void ag_pixel(double x, double y)
{
  ag_transform(&x,&y);
  if (ag_on_p)
  {
    screen_pixel(x,y);
    screen_show();
  }
  if (apict_on_p())
  {
    apict_pixel(x,y);
  }
  if ( Ps_active ) ps_pixel(x,y);
}

//#define THICK_LINES

void ag_untransformed_line(double u, double v, double x, double y)
{
  if (ag_on_p)
  {
    screen_line(u,v,x,y);
#ifdef THICK_LINES
    screen_line(u+1.0,v,x+1.0,y);
    screen_line(u,v+1.0,x,y+1.0);
    screen_line(u+1.0,v+1.0,x+1.0,y+1.0);
#endif
    screen_show();
  }
  if (apict_on_p())
  {
    apict_line(u,v,x,y);
#ifdef THICK_LINES
    apict_line(u+1.0,v,x+1.0,y);
    apict_line(u,v+1.0,x,y+1.0);
    apict_line(u+1.0,v+1.0,x+1.0,y+1.0);
#endif
  }
  if ( Ps_active )
  {
    ps_line(u,v,x,y);
#ifdef THICK_LINES
    ps_line(u+1.0,v,x+1.0,y);
    ps_line(u,v+1.0,x,y+1.0);
    ps_line(u+1.0,v+1.0,x+1.0,y+1.0);
#endif
  }
}

void ag_line(double u, double v, double x, double y)
{
  ag_transform(&u,&v);
  ag_transform(&x,&y);
  ag_untransformed_line(u,v,x,y);
}

void ag_untransformed_circle(double x, double y, double r)
{
  if (ag_on_p)
  {
    screen_circle(x,y,r);
    screen_show();
  }
  if (apict_on_p())
    apict_circle(x,y,r);
  if ( Ps_active ) ps_circle(x,y,r);
}

void ag_untransformed_disc(double x, double y, double r)
{
  if (ag_on_p)
  {
    screen_disc(x,y,r);
    screen_show();
  }
  if (apict_on_p())
    apict_disc(x,y,r);
  if ( Ps_active ) ps_disc(x,y,r);
}

void ag_circle(double x, double y, double r)
{
  ag_transform(&x,&y);
  ag_untransformed_circle(x,y,r);
}

void ag_circle_transformed_radius(double x, double y, double r)
{
  double xedge = x;
  double yedge = y+r;
  ag_transform(&x,&y);
  ag_transform(&xedge,&yedge);
  r = yedge - y;

  if (ag_on_p)
  {
    screen_circle(x,y,r);
    screen_show();
  }
  if (apict_on_p())
    apict_circle(x,y,r);
  if ( Ps_active ) ps_circle(x,y,r);
}

void ag_box(double x, double y, double u, double v)
{
  ag_transform(&x,&y);
  ag_transform(&u,&v);
  if (ag_on_p)
  {
    screen_box(x,y,u-x,v-y);
    screen_show();
  }
  if (apict_on_p())
    apict_box(x,y,u-x,v-y);
  if ( Ps_active ) ps_filled_rectangle(x,y,u,v);
}

void ag_disc(double x, double y, double r)
{
  ag_transform(&x,&y);
  ag_untransformed_disc(x,y,r);
}

void ag_disc_transformed_radius(double x, double y, double r)
{
  double xedge = x;
  double yedge = y+r;
  ag_transform(&x,&y);
  ag_transform(&xedge,&yedge);
  r = yedge - y;

  if (ag_on_p)
  {
    screen_disc(x,y,r);
    screen_show();
  }
  if (apict_on_p())
    apict_disc(x,y,r);
  if ( Ps_active ) ps_disc(x,y,r);
}

void ag_print_absolute_coords(double x, double y, char *s)
{
  if ( s == NULL ) my_error("ag_print: NULL string");
  if (ag_on_p)
  {
    screen_print(s,x,y);
    screen_show();
  }
  if (apict_on_p())
    apict_print(s,x,y);
  if ( Ps_active ) ps_print(x,y,s);
}

void ag_print(double x, double y, char *s)
{
  ag_transform(&x,&y);
  ag_print_absolute_coords(x,y,s);
}

/* ================================================================= */
/* Miscellaneous functions.                                          */
/* ================================================================= */

void ag_cross(double x, double y)
{
  double c = 10.0;
  ag_transform(&x,&y);
  ag_untransformed_line(x-c,y-c,x+c,y+c);
  ag_untransformed_line(x-c,y+c,x+c,y-c);
}

void ag_big_cross(double x, double y)
{
  double c = 20.0;
  int i,j;
  ag_transform(&x,&y);
  for ( i = -1 ; i <= 1 ; i++ )
    for ( j = -1 ; j <= 1 ; j++ )
    {
      ag_untransformed_line(x-c+i,y-c+j,x+c+i,y+c+j);
      ag_untransformed_line(x-c+i,y+c+j,x+c+i,y-c+j);
    }
}

void ag_smallbox(double x, double y)
{
  double c = 10.0;
  ag_transform(&x,&y);
  ag_untransformed_line(x-c,y-c,x+c,y-c);
  ag_untransformed_line(x+c,y+c,x+c,y-c);
  ag_untransformed_line(x+c,y+c,x-c,y+c);
  ag_untransformed_line(x-c,y-c,x-c,y+c);
}

void ag_smalltriangle(double x, double y)
{
  double c = 10.0;
  ag_transform(&x,&y);
  ag_untransformed_line(x-c,y-c,x+c,y-c);
  ag_untransformed_line(x-c,y-c,x,y+c);
  ag_untransformed_line(x+c,y-c,x,y+c);
}

void ag_smallcircle(double x, double y)
{
  double c = 10.0;
  ag_transform(&x,&y);
  ag_untransformed_circle(x,y,c);
}

void ag_fancy_disc(double x, double y)
{
  int col = ag_pen_color();
  ag_transform(&x,&y);

  ag_set_pen_color(AG_BLACK);
  ag_untransformed_disc(x,y,10.0);
  ag_set_pen_color(AG_WHITE);
  ag_untransformed_disc(x,y,8.0);
  ag_set_pen_color(col);
  ag_untransformed_disc(x,y,6.0);
}

void ag_mark(double x,double y,int marktype)
{
  switch ( marktype )
  {
    case DOT_MARKTYPE: ag_dot(x,y); break;
    case CROSS_MARKTYPE: ag_cross(x,y); break;
    case BIG_CROSS_MARKTYPE: ag_big_cross(x,y); break;
    case SMALLBOX_MARKTYPE: ag_smallbox(x,y); break;
    case SMALLTRIANGLE_MARKTYPE: ag_smalltriangle(x,y); break;
    case SMALLCIRCLE_MARKTYPE: ag_smallcircle(x,y); break;
    case PIXEL_MARKTYPE: ag_pixel(x,y); break;
    case FANCY_MARKTYPE: ag_fancy_disc(x,y); break;
    default: my_error("Unknown marktype");
  }
}

#define ARROW_DELTA 4.0

void ag_arrow(double x, double y, double u, double v)
{
  double mag;

  mag = sqrt( (x - u) * (x - u) + (y - v) * (y - v) );
  ag_line(x,y,u,v);

  if ( mag > 1e-4 )
  {
    double alpha = ARROW_DELTA * (u - x) / mag;
    double beta  = ARROW_DELTA * (v - y) / mag;

    ag_line(u - alpha - beta,v - beta + alpha,u,v);
    ag_line(u - alpha + beta,v - beta - alpha,u,v);
/*
    ag_line(u - alpha - beta,v - beta + alpha,
            u - alpha + beta,v - beta - alpha
           ); */
  }
}

void ag_rectangle(double x, double y, double u, double v)
{
  ag_line(x,y,x,v);
  ag_line(x,v,u,v);
  ag_line(u,v,u,y);
  ag_line(u,y,x,y);
}

void ag_clipped_line(double p, double q, double u, double v)
{
  if ( Ps_stream == NULL )
    ag_line(p,q,u,v);
  else
  {
    double minx = 0.0;
    double maxx = 512.0;
    double miny = 0.0;
    double maxy = 512.0;

    if ( p > u )
    {
      double t = p;
      p = u;
      u = t;
      t = q;
      q = v;
      v = t;
    }

    if ( u < minx || p > maxx ) return;

    if ( u - p < 1e-5 )
    {
      if ( q > v )
      {
        double t = q;
        q = v;
        v = t;
      }
      ag_line(p,real_max(miny,q),u,real_min(maxy,v));
    }
    else
    {
      double m = (q - v) / (p - u);
      double c = q - m * p;

      if ( p < minx )
      {
        p = minx;
        q = p * m + c;
      }

      if ( u > maxx )
      {
        u = maxx;
        v = u * m + c;
      }

      if ( q > v )
      {
        double t = p;
        p = u;
        u = t;
        t = q;
        q = v;
        v = t;
      }

      if ( v < miny || q > maxy ) return;

      if ( v - q < 1e-5 )
        ag_line(p,q,u,v);
      else
      {
        m = (p - u) / (q - v);
        c = m * q - p;

        if ( q < miny )
        {
          q = miny;
          p = m * q + c;
        }

        if ( v > maxy )
        {
          v = maxy;
          u = m * v + c;
        }

        ag_line(p,q,u,v);
      }
    }
  }
}

void ag_clipped_dot(double x, double y)
{
  if ( x > 0.0 && x < 512.0 && y > 0.0 && y < 512.0 )
    ag_dot(x,y);
}

void ag_uparrow(double x, double y)
{
  ag_line(x-ARROW_DELTA/2.0,y-ARROW_DELTA/2.0,x,y);
  ag_line(x+ARROW_DELTA/2.0,y-ARROW_DELTA/2.0,x,y);
}

void ag_downarrow(double x, double y)
{
  ag_line(x-ARROW_DELTA/2.0,y+ARROW_DELTA/2.0,x,y);
  ag_line(x+ARROW_DELTA/2.0,y+ARROW_DELTA/2.0,x,y);
}

void ag_rounded_box(double x, double y, double u, double v, double r)
{
  ag_box(x-r,y,u+r,v);
  ag_box(x,y-r,u,v+r);
  ag_disc(x,y,r);
  ag_disc(x,v,r);
  ag_disc(u,y,r);
  ag_disc(u,v,r);
}

/* Backwards compatibility with earlier versions.  If called with 0.0
   it sets the pen to black; if called with 1.0 it sets the pen to
   white; if called with anything else, it sets the pen to purple. */
void ag_pen(double c)
{
  if (c == 0.0)
    simple_set_pen_color(AG_BLACK);
  else if (c == 1.0)
    simple_set_pen_color(AG_WHITE);
  else
    simple_set_pen_color(AG_PURPLE);
  if ( Ps_active )
  {
    if ( c == 0.0 )
      ps_color(AG_BLACK);
    else if (c == 1.0)
      ps_color(AG_WHITE);
    else
      ps_color(AG_PURPLE);
  }
}

/* ===================================================================== */
/* General notes on color handling.                                      */
/*

   Low-level color handling for the amut graphics is platform-dependent,
   but the functional interface is the same across different platforms:

   By default the "pen" will draw in black.  The user can specify
   another color for the pen by calling ag_set_pen_color.  This
   function takes an integer argument from 0 through NUM_COLORS - 1.

   Currently amut graphics only provides a small palette of colors.
   Three methods are provided to aid the user in selecting colors:
   1) The integer constants AG_BLACK, AG_WHITE, AG_RED, AG_YELLOW,
      AG_GREEN, AG_CYAN, AG_BLUE, and AG_PURPLE are hash-defined to
      appropriate shades.
   2) The user can also call ag_spectrum_color with a double in
      the range [0.0, 1.0]; this will then return the integer that is
      closest to the specified fraction of the way along a rainbow
      spectrum running from red (0.0) to purple (1.0).  Note that
      white is excluded from this spectrum.
   3) The user can call ag_rgb_color with a three byte integer.
      The first byte specifies the desired intensity for the red
      component of the color (0 = off; 11111111 = full on; first bit
      most significant).  The second byte specifies the desired
      intensity for the green component.  The third byte specifies the
      desired intensity for the blue component.  ag_rgb_color
      returns the number of the closest color in the amut palette.
      *NOTE*: this function is not implemented in the PC/Windows NT
      world.                                                             */
/* ===================================================================== */

/* Sets the pen color to the specified amut color (legal values being
   the integers from 0 through NUM_COLORS - 1. */
void ag_set_pen_color(int amut_col)
{
  simple_set_pen_color(amut_col);

  if ( Ps_active )
  {
    if (PS_MONO == 2)
      {
	if (amut_col == AG_BLACK) {
	  fprintf(Ps_stream, "0 setgray\n");
	} else if (amut_col == AG_WHITE) {
	  fprintf(Ps_stream, "1 setgray\n");
	} else
	  fprintf(Ps_stream, "/Pattern setcolorspace %s setcolor\n",
		  pspatnames[amut_col % NUM_PSPATS]);
      }
    else if (PS_MONO == 1)
      {
	ps_color(amut_col);
      }
    else
      {
	double r,g,b;

	ag_rgb(amut_col, &r, &g, &b);
	ps_rgb(r, g, b);
      }
  }
}

/* Returns the current pen color */
int ag_pen_color()
{
  return(simple_pen_color());
}

/* Calling spectrum_color with a double in the range [0, 1.0] will
   return the integer that is closest to the specified fraction of the
   way along the amut approximate rainbow spectrum described in
   amxw.c.  0.0 corresponds to the red end and 1.0 to the violet end.
        Note that white is excluded from the sepctrum.  */
int ag_spectrum_color(double fract)
{
  return(simple_spectrum_color(fract));
}

/* The user can call rgb_color with a three byte integer.  The
   first byte specifies the desired intensity for the red component of
   the color (0 = off; 11111111 = full on; first bit most
   significant).  The second byte specifies the desired intensity for
   the green component.  The third byte specifies the desired
   intensity for the blue component.
        rgb_color returns the number of the closest color in the amut
   palette.  Closeness is determined by the maximum difference between
   R, G, and B values in the two colors.

   NOTE: not implemented for PC/Windows NT world.
*/
void ag_rgb(int ag_color, double *r_r, double *r_g, double *r_b)
{
  void simple_rgb(int ag_color, double *r_r, double *r_g, double *r_b);

  simple_rgb(ag_color, r_r, r_g, r_b);
}

/* ===================================================================== */
/* Functions related to mouse input.                                     */
/* ===================================================================== */

int ag_get_xy(double *x_ptr, double *y_ptr)
{
  int button = simple_get_xy("",x_ptr,y_ptr);
  ag_untransform(x_ptr,y_ptr);
  return(button);
}

int ag_get_xy_prompted(char *prompt,double *x_ptr, double *y_ptr)
{
  int button = simple_get_xy(prompt,x_ptr,y_ptr);
  ag_untransform(x_ptr,y_ptr);
  return(button);
}

void click_on_message(double xlo, double ylo, double xhi, double yhi, char *m)
{
  double x = xhi + 10.0;
  double y = yhi + 10.0;
  ag_pen(1.0);
  ag_box(xlo,ylo,xhi,yhi);
  ag_pen(0.0);
  ag_rectangle(xlo,ylo,xhi,yhi);
  ag_print(xlo+5.0,ylo+5.0,m);

  while ( x < xlo || y < ylo || x > xhi || y > yhi )
    ag_get_xy(&x,&y);

  ag_box(xlo,ylo,xhi,yhi);
}

void wait_for_click(void)
{
  click_on_message(20.0,20.0,280.0,60.0,"Click here to continue..");
}

void ag_print_vertical(double x,double y,char *string)
{
  int i;
  for ( i = 0 ; string[i] != '\0' ; i++ )
  {
    char buff[100];
    char c = string[i];
    if ( c == '-' || c == '_' ) c = '|';
    else if ( c == '|' ) c = '-';

    sprintf(buff,"%c",string[i]);
    ag_print(x,y-17.0 * i,buff);
  }
}

void amgr_wait_for_key()
{
  printf("Hit Return.\n");
  while ( getchar() != '\n' )
  {
    /* skip */
  }
}

void test_amgr(int argc,char *argv[])
{
  int i;

  printf("Will turn on graphics after key press...\n");
  amgr_wait_for_key();
  ag_on("");
  printf("Will draw line after key press...\n");
  amgr_wait_for_key();
  ag_line(0.0,0.0,512.0,256.0);
  printf("Will draw red circle after key press...\n");
  amgr_wait_for_key();
  ag_set_pen_color(AG_RED);
  ag_disc(128.0,384.0,50.0);
  printf("Will ag_off() after key...\n");
  amgr_wait_for_key();
  ag_set_pen_color(AG_WHITE);
  ag_box(90.0,340.0,280.0,450.0);
  ag_set_pen_color(AG_GRAY);
  ag_rectangle(90.0,390.0,280.0,450.0);
  ag_set_pen_color(AG_PURPLE);
  ag_print(100.0,400.0,"I'm a");
  ag_print(200.0,420.0,"Jumbo");
  ag_set_pen_color(AG_GREEN);
  ag_print(160.0,410.0,"Jelly");
  ag_off();
  printf("Bye!\n");
  amgr_wait_for_key();

  ag_on("plop.ps");
  for ( i = 195 ; i > 0 ; i-= 5 )
  {
    ag_set_pen_color(ag_spectrum_color(i / 200.0));
    ag_disc(150.0+i/5.0,250.0+i/10.0,(double)i);
  }
  ag_off();
}

/* Is there a point r=(x,y) on the line segment between
   (xa,ya) (xb,yb) for which (r - (xp,yp)).(dx,dy) == 0?

   If so, return TRUE and put result in *r_x and *r_y
*/
bool pivhelp(double xa,double ya,double xb,double yb,
               double xp,double yp,
               double xd,double yd,
               double *r_x,double *r_y)
{
  double num = (xp - xa) * xd + (yp - ya) * yd;
  double den = (xb - xa) * xd + (yb - ya) * yd;
  double lambda = (fabs(den) > 1e-6) ? (num / den) :
                  (den < 0.0) ? (num / -1e-6) : (num / 1e-6);
  bool result = lambda >= 0.0 && lambda <= 1.0;

  *r_x = xa + lambda * (xb - xa);
  *r_y = ya + lambda * (yb - ya);

  return(result);
}

/* Draws a line segment in the ag graphics window consisting of
   that set of points { (x,y) | ( (x,y) - (px,py) ) . (dx,dy) == 0 }.

   This is a line perpendicular to (dx,dy) that travels through the point
   (px,py). */
void ag_cut(double px,double py,double dx,double dy)
{
  double x1 = -77.7;
  double x2 = -77.7;
  double y1 = -77.7;
  double y2 = -77.7;
  double x,y;
  bool defined1 = FALSE;
  bool defined2 = FALSE;

  if ( pivhelp(Frame_x_lo,Frame_y_lo,Frame_x_hi,Frame_y_lo,px,py,dx,dy,&x,&y) )
  { if ( defined1 ) { x2 = x; y2 = y; defined2 = TRUE; }
    else { x1 = x; y1 = y; defined1 = TRUE; }
  }

  if ( pivhelp(Frame_x_hi,Frame_y_hi,Frame_x_hi,Frame_y_lo,px,py,dx,dy,&x,&y) )
  { if ( defined1 ) { x2 = x; y2 = y; defined2 = TRUE; }
    else { x1 = x; y1 = y; defined1 = TRUE; }
  }

  if ( pivhelp(Frame_x_hi,Frame_y_hi,Frame_x_lo,Frame_y_hi,px,py,dx,dy,&x,&y) )
  { if ( defined1 ) { x2 = x; y2 = y; defined2 = TRUE; }
    else { x1 = x; y1 = y; defined1 = TRUE; }
  }

  if ( pivhelp(Frame_x_lo,Frame_y_lo,Frame_x_lo,Frame_y_hi,px,py,dx,dy,&x,&y) )
  { if ( defined1 ) { x2 = x; y2 = y; defined2 = TRUE; }
    else { x1 = x; y1 = y; defined1 = TRUE; }
  }

  if ( defined1 && defined2 )
    ag_line(x1,y1,x2,y2);
}


agbox *mk_agbox(double xlo,double ylo,double xhi,double yhi)
{
  agbox *agb = AM_MALLOC(agbox);
  agb -> xlo = xlo;
  agb -> ylo = ylo;
  agb -> xhi = xhi;
  agb -> yhi = yhi;
  return agb;
}

void free_agbox(agbox *agb)
{
  AM_FREE(agb,agbox);
}

agbox *mk_copy_agbox(agbox *agb)
{
  agbox *newagb = AM_MALLOC(agbox);
  newagb -> xlo = agb -> xlo;
  newagb -> ylo = agb -> ylo;
  newagb -> xhi = agb -> xhi;
  newagb -> yhi = agb -> yhi;
  return newagb;
}

void fprintf_agbox(FILE *s,char *m1,agbox *x,char *m2)
{
  fprintf(s,"%s -> xlo = %g%s\n",m1,x->xlo,m2);
  fprintf(s,"%s -> ylo = %g%s\n",m1,x->ylo,m2);
  fprintf(s,"%s -> xhi = %g%s\n",m1,x->xhi,m2);
  fprintf(s,"%s -> yhi = %g%s\n",m1,x->yhi,m2);
}

double break_agval(double lo,double hi,double frac)
{
  return lo + frac * (hi - lo);
}

double vertical_break_agval(agbox *agb,double frac)
{
  return break_agval(agb->ylo,agb->yhi,frac);
}

double horizontal_break_agval(agbox *agb,double frac)
{
  return break_agval(agb->xlo,agb->xhi,frac);
}

void agbox_piebar(agbox *agb,dyv *widths
                  ,string_array *left,string_array *right
                  ,ivec *colors){
  int num_segs = dyv_size(widths);
  int num_colors = ivec_size(colors);
  agbox *i_agb, *i_agb_right, *i_agb_left;
  double x = 0, dx, xx;
  double width = dyv_sum(widths);
  char *lstr, *rstr;
  int i,count=0;
  bool alternate_colors = num_colors<num_segs;

  for(i=0;i<num_segs;i++){
    if(dyv_ref(widths,i)>0){
      dx = dyv_ref(widths,i);
      lstr = string_array_ref(left,i);
      rstr = string_array_ref(right,i);
      i_agb = mk_agbox_corresponding_to_time(agb,0,width,x,x+dx);

      agbox_color_fill(i_agb,ivec_ref(colors,(count%num_colors)));
      ag_set_pen_color(AG_BLACK);
      x += dx;

      if(((strlen(rstr)+4)*AG_CHAR_WIDTH) < agbox_width(i_agb)){
        xx = i_agb->xhi-(strlen(rstr)+1)*AG_CHAR_WIDTH;
        i_agb_right = mk_vertical_stripe_absolute(i_agb,xx,i_agb->xhi);
        i_agb_left = mk_vertical_stripe_absolute(i_agb,i_agb->xlo,xx);
        agbox_render_string(i_agb_right,rstr);
      } else {
        i_agb_right = NULL;
        i_agb_left = mk_copy_agbox(i_agb);
      }
      agbox_render_string_if_enough_room(i_agb_left,lstr);
      agbox_border(i_agb);

      if(i_agb_right) free_agbox(i_agb_right);
      free_agbox(i_agb_left);
      free_agbox(i_agb);
      count++;
    }
    else if(!alternate_colors) count++;
  }
}

double ag_value_corresponding_to_time(agbox *parent,
				      double parent_start_time,
				      double parent_end_time,
				      double local_time)
{
  double z = (local_time - parent_start_time) /
             (parent_end_time - parent_start_time);
  return vertical_break_agval(parent,z);
}

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
				      double local_end_time)
{
  double d = parent_end_time - parent_start_time;
  double fracstart = (local_start_time - parent_start_time) / d;
  double fracend = (local_end_time - parent_start_time) / d;
  return mk_vertical_stripe_relative(parent,fracstart,fracend);
}

double get_time_corresponding_to_agbox(agbox *agb,double start,double end,double x){
  double agb_x = (x - agb->xlo) / agbox_width(agb);
  return start + agb_x * (end - start);
}

/* Prints horizontally at top left of box */
void agbox_render_string_horizontal(agbox *agb,char *string)
{
  ag_print(agb->xlo+0.4*AG_CHAR_WIDTH,
	   agb->yhi-AG_CHAR_HEIGHT-1.0,string);
}

/* Prints vertically at top left of box */
void agbox_render_string_vertical(agbox *agb,char *string)
{
  ag_print_vertical(agb->xlo+0.5*AG_CHAR_WIDTH,
		    agb->yhi-AG_CHAR_HEIGHT-1.0,string);
}

/* Prints horizontally at top left of box, unless it can't fit the string
   in horizontally, and it can fit more letters in vertically, in which
   case prints vertically.  It does clip letters.*/
void agbox_render_string(agbox *agb,char *string)
{
  int numchars = (int)strlen(string);
  int numchars_h = int_min(numchars,(int)(agbox_height(agb)/(AG_CHAR_HEIGHT+1.0)));
  int numchars_w = int_min(numchars,(int)(agbox_width(agb)/AG_CHAR_WIDTH));
  char *buff = AM_MALLOC_ARRAY(char,numchars+1);

  sprintf(buff,"%s",string);

  if ( numchars_w < numchars && numchars_h > numchars_w ) {
    buff[numchars_h] = '\0';
    agbox_render_string_vertical(agb,buff);
  } else {
    buff[numchars_w] = '\0';
    agbox_render_string_horizontal(agb,buff);
  }

  AM_FREE_ARRAY(buff,char,numchars+1);
}

bool enough_room_for_string_in_agbox(agbox *agb,char *string)
{
  int len = (int)strlen(string);
  double box_width = agbox_width(agb);
  double box_height = agbox_height(agb);
  bool ok = box_width >= AG_CHAR_WIDTH &&
            box_height >= AG_CHAR_HEIGHT;
  if ( ok )
  {
    double real_width = (len+.4) * (double)AG_CHAR_WIDTH;
    double real_height =  (len+.5) * (double)AG_CHAR_HEIGHT;
    ok = real_width < box_width || real_height < box_height;
  }

  return ok;
}

/* As above, except only renders if enough room */
void agbox_render_string_if_enough_room(agbox *agb,char *string)
{
  if ( enough_room_for_string_in_agbox(agb,string) )
    agbox_render_string(agb,string);
}

/* As above, except returns the leftover horizontal space*/
bool agbox_render_string_if_enough_room_and_resize(agbox **agb,char *string)
{
  double width = AG_CHAR_WIDTH*((double) (strlen(string)+1));
  agbox *agb_old;
  bool ret = FALSE;

  if(width<agbox_width(*agb) && AG_CHAR_HEIGHT<agbox_height(*agb)){
    ret = TRUE;
    agbox_render_string_horizontal(*agb,string);
    agb_old = *agb;
    *agb = mk_vertical_stripe_absolute(agb_old,agb_old->xlo+width,agb_old->xhi);
    free_agbox(agb_old);
  }
  return ret;
}

int render_string_array_horizontal_if_possible(agbox *agb,string_array *sa){
  string_array *used_sa = mk_string_array(0);
  double width = agbox_width(agb);
  double size = 0.0;
  char *str = NULL;
  agbox *agb_next;
  int ok = 0;

  while(size<width && string_array_size(sa)!=0){
    str = string_array_ref(sa,0);
    size += (strlen(str)+1)*AG_CHAR_WIDTH;
    if(size<width){
      add_to_string_array(used_sa,str);
      string_array_remove(sa,0);
    }
  }

  if(string_array_size(used_sa)>0 && size>width){
    if(agbox_height(agb)>(AG_CHAR_HEIGHT*2)){
      agb_next = mk_horizontal_stripe_absolute(agb,agb->ylo,agb->yhi-AG_CHAR_HEIGHT);
      ok = render_string_array_horizontal_if_possible(agb_next,sa);
      free_agbox(agb_next);
    } else {
      ok = 0;
    }
  }

  if(size<width || ok){
    str = mk_string_from_string_array(used_sa);
    agbox_render_string_horizontal(agb,str);
    free_string(str);
    free_string_array(used_sa);
    return ok+1;
  }
  free_string_array(used_sa);
  return 0;
}

int render_string_array_vertical_if_possible(agbox *agb,string_array *sa){
  string_array *used_sa = mk_string_array(0);
  double height = agbox_height(agb);
  double size = 0.0;
  char *str;
  agbox *agb_next;
  int ok = 0;

  while(size<height && string_array_size(sa)!=0){
    str = string_array_ref(sa,0);
    size += (strlen(str)+1)*AG_CHAR_HEIGHT;
    if(size<height){
      add_to_string_array(used_sa,str);
      string_array_remove(sa,0);
    }
  }

  if(string_array_size(used_sa)>0 && size>height){
    if(agbox_width(agb)>(AG_CHAR_WIDTH*3)){
      agb_next = mk_vertical_stripe_absolute(agb,agb->xlo+(AG_CHAR_WIDTH*2.0),agb->xhi);
      ok = render_string_array_vertical_if_possible(agb_next,sa);
      free_agbox(agb_next);
    } else {
      ok = 0;
    }
  }

  if(size<height || ok){
    str = mk_string_from_string_array(used_sa);
    agbox_render_string_vertical(agb,str);
    free_string(str);
    free_string_array(used_sa);
    return ok+1;
  }
  free_string_array(used_sa);
  return 0;
}

/*Fits each string into the box the best way it knows how.  Returns how many rows/cols it
  needed to fit them all, or 0 if it couldn't do it.*/
int agbox_render_string_array(agbox *agb,string_array *sa){
  string_array *copy_sa = NULL;
  int ok;

  copy_sa = mk_copy_string_array(sa);
  ok = render_string_array_horizontal_if_possible(agb,copy_sa);
  free_string_array(copy_sa);
  if(!ok){
    copy_sa = mk_copy_string_array(sa);
    ok = render_string_array_vertical_if_possible(agb,copy_sa);
    free_string_array(copy_sa);
  }
  return ok;
}

void agbox_render_string_or_abbrev(agbox *agb,char *name,char *shortname)
{
  if ( enough_room_for_string_in_agbox(agb,name) )
    agbox_render_string(agb,name);
  else if ( shortname != NULL && enough_room_for_string_in_agbox(agb,shortname) )
    agbox_render_string(agb,shortname);
}


void agbox_render_string_at_bottom(agbox *agb,char *string)
{
  agbox *bottom =
    mk_horizontal_stripe_absolute(agb,agb->ylo,
				  agb->ylo + 1.5 * AG_CHAR_HEIGHT);
  agbox_render_string_horizontal(bottom,string);
  free_agbox(bottom);
}

agbox *mk_horizontal_stripe_absolute(agbox *parent,double aglo,double aghi)
{
  return mk_agbox(parent->xlo,aglo,parent->xhi,aghi);
}

agbox *mk_horizontal_stripe_relative(agbox *parent,double fraclo,double frachi)
{
  return mk_horizontal_stripe_absolute(parent,
				       vertical_break_agval(parent,fraclo),
				       vertical_break_agval(parent,frachi));
}

/* Imagine dividing the first arg (parent) up into "numstrips" horizontal
   strips of equal height, exactly covering and partitioning parent.
   What would be the agbox corresponding to the stripnum'th of these stripes
   from the top (where stripnum==0 corresponds to the topmost stripe?) */
agbox *mk_agbox_horizontal_substripe(agbox *parent,int stripnum,int numstrips)
{
  double x = (double) numstrips;
  double y = x - stripnum - 1.0;
  return mk_horizontal_stripe_relative(parent,y/x,(y+1)/x);
}

agbox *mk_vertical_stripe_absolute(agbox *parent,double aglo,double aghi)
{
  return mk_agbox(aglo,parent->ylo,aghi,parent->yhi);
}

agbox *mk_vertical_stripe_relative(agbox *parent,double fraclo,double frachi)
{
  return mk_vertical_stripe_absolute(parent,
				     horizontal_break_agval(parent,fraclo),
				     horizontal_break_agval(parent,frachi));
}

agbox *mk_agbox_vertical_substripe(agbox *parent,int stripnum,int numstrips)
{
  return mk_vertical_stripe_relative(parent,
				     stripnum/(double) numstrips,
				     (stripnum+1)/(double) numstrips);
}

/*m must be an inverted slope, i.e. inv_slope = 1 / slope*/
void clip_point_vertical(double *x,double *y,double c,double m){
  *x = *x + (c - *y) * m;
  *y = c;
}

void clip_point_horizontal(double *x,double *y,double c,double m){
  *y = *y + (c - *x) * m;
  *x = c;
}

/* Draws a line segment in absolute coordinates.
   It DOES clip the line so that only its components INSIDE agb are
   drawn. */
void render_line_segment_absolute(agbox *agb,
				  double x0,double y0,double x1,double y1)
{
  double slope = 0, inv_slope = 0;
  double temp;

  if(x0>x1){
    temp = x0; x0 = x1; x1 = temp;
    temp = y0; y0 = y1; y1 = temp;
  }
  if(y0!=y1){
    inv_slope = (x1 - x0) / (y1 - y0);
    if(y0<agb->ylo) clip_point_vertical(&x0,&y0,agb->ylo,inv_slope);
    else if(y0>agb->yhi) clip_point_vertical(&x0,&y0,agb->yhi,inv_slope);
    if(y1<agb->ylo) clip_point_vertical(&x1,&y1,agb->ylo,inv_slope);
    else if(y1>agb->yhi) clip_point_vertical(&x1,&y1,agb->yhi,inv_slope);
  }
  if(x0!=x1){
    slope = (y1 - y0) / (x1 - x0);
    if(x0<agb->xlo) clip_point_horizontal(&x0,&y0,agb->xlo,slope);
    if(x1>agb->xhi) clip_point_horizontal(&x1,&y1,agb->xhi,slope);
  }

  if((y0!=y1 || (y0>=agb->ylo && y0<=agb->yhi)) && (x0!=x1 || (x0>=agb->xlo && x0<=agb->xhi)))
    ag_line(x0,y0,x1,y1);
}

/* Draws a relative-coordinate line segment. Assume that the
   bottom left of the agbox should be identified with the
   relative coordinate (0,0) and the top right should be identified
   with (1,1).

   It DOES clip the line so that only its components INSIDE agb are
   drawn. */
void render_line_segment_relative(agbox *agb,
				  double x0,double y0,double x1,double y1)
{
  render_line_segment_absolute(agb,
			       horizontal_break_agval(agb,x0),
			       vertical_break_agval(agb,y0),
			       horizontal_break_agval(agb,x1),
			       vertical_break_agval(agb,y1));
}

/* Draws a virtual-coordinate line segment. Assume that the
   bottom left of the agbox should be identified with the
   virtual-coordinates-point (vx0,vy0) and the top right with
   (vx1,vy1). Then this draws the line segment that begins with
   (x0,y0) in virtual coordinates and ends with (x1,y1).

   It DOES clip the line so that only its components INSIDE agb are
   drawn. */
void render_virtual_line_segment(agbox *agb,
				 double vx0,double vy0,double vx1,double vy1,
				 double x0,double y0,double x1,double y1)
{
  render_line_segment_relative(agb,(x0-vx0)/(vx1-vx0),(y0-vy0)/(vy1-vy0),
			           (x1-vx0)/(vx1-vx0),(y1-vy0)/(vy1-vy0));
}

double agbox_width(agbox *agb)
{
  return agb->xhi - agb->xlo;
}

double agbox_height(agbox *agb)
{
  return agb->yhi - agb->ylo;
}

double agbox_thinnest_dimension(agbox *agb)
{
  return real_min(agbox_width(agb),agbox_height(agb));
}

void agbox_border(agbox *agb)
{
  ag_rectangle(agb->xlo,agb->ylo,agb->xhi,agb->yhi);
}

void agbox_color_border(agbox *agb,int agcolor)
{
  ag_set_pen_color(agcolor);
  agbox_border(agb);
}

void agbox_fill(agbox *agb){
  ag_box(agb->xlo,agb->ylo,agb->xhi,agb->yhi);
}

void agbox_color_fill(agbox *agb,int agcolor){
  ag_set_pen_color(agcolor);
  agbox_fill(agb);
}

bool is_in_agbox(agbox *agb,double x,double y){
  return (x>=agb->xlo && x<=agb->xhi && y>=agb->ylo && y<=agb->yhi);
}

/* If you are working in an application that is driven by dispatch_xxx( )
   type functions, it is often worth using special_ag_on instead of ag_on.
   This is because when the comment on sessionname is invoked by the
   user, all graphs will automatically be saved in auto-named .ps files
   and linked within a latex document into a session log. It is harmless
   to use this instead of ag_on if you're not sure. */
void special_ag_on(char *filename)
{
  bool comment_on = comment_is_on();
  if ( !comment_on )
    ag_on(filename);
  else
  {
    extern int Graphic_number; /* from command.c */
    bool save_to_file = filename != NULL && !eq_string(filename,"");
    char *savename;

    if ( save_to_file )
      savename = mk_copy_string(filename);
    else
      savename = mk_printf("%s.%04d.ps",comment_session_name(),Graphic_number);

    ag_on(savename);

    printf("\n----------------------------------------------------------\n");
    printf("\\end{verbatim}\n");
    printf("\\bigfig{%s}\n",savename);
    printf("\\begin{verbatim}\n");

    Graphic_number += 1;

    free_string(savename);
  }
}


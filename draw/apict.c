/*
       File:             apict.c
       Author:           Mary Soon Lee
       Created:          16 Jan 96
       Modified:         3 Nov 97

       Description:      Code supporting apicts, the structures used
                         to store pictures.  Elements of these pictures
                         may be circles, lines, strings, rectangles,
                         and discs (solid circles).

           Andrew added "apict_set" data structure. Contains multiple
           apicts.
*/

#include <stdio.h>
#include <math.h>
#include "./utils/ambs.h"
#include "apict.h"
#include "amgr.h"

/* ==================================================================== */
/* Basic functions to allocate/free/print apicts and their components.  */
/* ==================================================================== */

/* Allocates and returns an empty apict of default size
   (width = height = 512). */
apict *mk_apict()
{
  apict *result = AM_MALLOC(apict);
  result->width = 512;
  result->height = 512;
  result->lines = NULL;
  result->circles = NULL;
  result->strings = NULL;
  result->rects = NULL;
  result->discs = NULL;
  result->pixels = NULL;
  return(result);
}

void am_free_ap_line(ap_line *l)
{
  AM_FREE(l, ap_line);
}

void am_free_ap_circle(ap_circle *c)
{
  AM_FREE(c, ap_circle);
}

void am_free_ap_string(ap_string *s)
{
  AM_FREE(s, ap_string);
}

void am_free_ap_rect(ap_rect *r)
{
  AM_FREE(r, ap_rect);
}

void am_free_ap_disc(ap_disc *d)
{
  AM_FREE(d, ap_disc);
}

void am_free_ap_pixel(ap_pixel *d)
{
  AM_FREE(d, ap_pixel);
}

void am_free_ap_circles(ap_circles *circles_ptr)
{
  while ( circles_ptr != NULL )
  {
    ap_circles *next = circles_ptr->rest_circles;
    am_free_ap_circle(circles_ptr->c1);
    AM_FREE(circles_ptr, ap_circles);
    circles_ptr = next;
  }
}

void am_free_ap_lines(ap_lines *lines_ptr)
{
  while ( lines_ptr != NULL )
  {
    ap_lines *next = lines_ptr->rest_lines;
    am_free_ap_line(lines_ptr->l1);
    AM_FREE(lines_ptr, ap_lines);
    lines_ptr = next;
  }
}

void am_free_ap_strings(ap_strings *strings_ptr)
{
  while ( strings_ptr != NULL )
  {
    ap_strings *next = strings_ptr->rest_strings;
    am_free_ap_string(strings_ptr->s1);
    AM_FREE(strings_ptr, ap_strings);
    strings_ptr = next;
  }
}

void am_free_ap_rects(ap_rects *rects_ptr)
{
  while ( rects_ptr != NULL )
  {
    ap_rects *next = rects_ptr->rest_rects;
    am_free_ap_rect(rects_ptr->r1);
    AM_FREE(rects_ptr, ap_rects);
    rects_ptr = next;
  }
}

void am_free_ap_discs(ap_discs *discs_ptr)
{
  while ( discs_ptr != NULL )
  {
    ap_discs *next = discs_ptr->rest_discs;
    am_free_ap_disc(discs_ptr->d1);
    AM_FREE(discs_ptr, ap_discs);
    discs_ptr = next;
  }
}

void am_free_ap_pixels(ap_pixels *pixels_ptr)
{
  while ( pixels_ptr != NULL )
  {
    ap_pixels *next = pixels_ptr->rest_pixels;
    am_free_ap_pixel(pixels_ptr->p1);
    AM_FREE(pixels_ptr, ap_pixels);
    pixels_ptr = next;
  }
}

/* AWM- N.B. I'm changing the preferred name of this to free_apict, which
   is a SPR/Auton convention. am_free_apict retained for compatibility
*/
void free_apict(apict *g)
{
  am_free_ap_lines(g->lines);
  am_free_ap_circles(g->circles);
  am_free_ap_strings(g->strings);
  am_free_ap_rects(g->rects);
  am_free_ap_discs(g->discs);
  am_free_ap_pixels(g->pixels);
  AM_FREE(g, apict);
}

void am_free_apict(apict *g)
{
  free_apict(g);
}

void print_ap_line(ap_line *l)
{
  printf("Line: (%g, %g) -> (%g, %g)\n", l->x1, l->y1, l->x2, l->y2);
}

void print_ap_circle(ap_circle *c)
{
  printf("Circle: center (%g, %g), radius %g\n", c->x, c->y, c->r);
}

void print_ap_string(ap_string *s)
{
  printf("String: start (%g, %g), text \"%s\"\n", s->x, s->y, s->str);
}

void print_ap_rect(ap_rect *r)
{
  printf("Rect: (%g, %g) -> (%g, %g)\n",
         r->x, r->y, r->x + r->w, r->y + r->h);
}

void print_ap_disc(ap_disc *d)
{
  printf("Disc: center (%g, %g), radius %g\n", d->x, d->y, d->r);
}

void print_ap_pixel(ap_pixel *d)
{
  printf("Pixel: location (%g, %g)\n", d->x, d->y);
}

void print_apict(apict *g)
{
  ap_lines *l = g->lines;
  ap_circles *c = g->circles;
  ap_strings *s = g->strings;
  ap_rects *r = g->rects;
  ap_discs *d = g->discs;
  ap_pixels *p = g->pixels;

  while (l != NULL)
  {
    print_ap_line(l->l1);
    l = l->rest_lines;
  }
  while (c != NULL)
  {
    print_ap_circle(c->c1);
    c = c->rest_circles;
  }
  while (s != NULL)
  {
    print_ap_string(s->s1);
    s = s->rest_strings;
  }
  while (r != NULL)
  {
    print_ap_rect(r->r1);
    r = r->rest_rects;
  }
  while (d != NULL)
  {
    print_ap_disc(d->d1);
    d = d->rest_discs;
  }
  while (p != NULL)
  {
    print_ap_pixel(p->p1);
    p = p->rest_pixels;
  }
}

ap_line *mk_ap_line(double x1, double y1, double x2, double y2, int col)
{
  ap_line *result = AM_MALLOC(ap_line);
  result->x1 = x1;
  result->x2 = x2;
  result->y1 = y1;
  result->y2 = y2;
  result->color = col;
  return(result);
}

ap_circle *mk_ap_circle(double x, double y, double r, int col)
{
  ap_circle *result = AM_MALLOC(ap_circle);
  result->x = x;
  result->y = y;
  result->r = r;
  result->color = col;
  return(result);
}

ap_string *mk_ap_string(double x, double y, char *s, int col)
{
  ap_string *result = AM_MALLOC(ap_string);
  result->x = x;
  result->y = y;
  if (((int) strlen(s)) > MAX_APICT_STR_LEN)
    my_error("mk_ap_string: Too long a string");
  sprintf(result->str, "%s", s);
  result->color = col;
  return(result);
}

ap_rect *mk_ap_rect(double x, double y, double w, double h, int col)
{
  ap_rect *result = AM_MALLOC(ap_rect);
  result->x = x;
  result->y = y;
  result->w = w;
  result->h = h;
  result->color = col;
  return(result);
}

ap_disc *mk_ap_disc(double x, double y, double r, int col)
{
  ap_disc *result = AM_MALLOC(ap_disc);
  result->x = x;
  result->y = y;
  result->r = r;
  result->color = col;
  return(result);
}

ap_pixel *mk_ap_pixel(double x, double y, int col)
{
  ap_pixel *result = AM_MALLOC(ap_pixel);
  result->x = x;
  result->y = y;
  result->color = col;
  return(result);
}

/* ==================================================================== */
/* Basic functions to add lines/circles/strings/rects/discs to an apict */
/* ==================================================================== */

void ap_add_line(apict *g, double x1, double y1, double x2, double y2,
                 int col)
{
  ap_lines *new_lines = AM_MALLOC(ap_lines);
  ap_line *nl = mk_ap_line(x1, y1, x2, y2, col);

  new_lines->l1 = nl;
  if (g->lines == NULL)
    new_lines->rest_lines = NULL;
  else
    new_lines->rest_lines = g->lines;
  g->lines = new_lines;
}

void ap_add_circle(apict *g, double x, double y, double r, int col)
{
  ap_circles *new_circles = AM_MALLOC(ap_circles);
  ap_circle *nc = mk_ap_circle(x, y, r, col);

  new_circles->c1 = nc;
  if (g->circles == NULL)
    new_circles->rest_circles = NULL;
  else
    new_circles->rest_circles = g->circles;
  g->circles = new_circles;
}

void ap_add_string(apict *g, double x, double y, char *s, int col)
{
  ap_strings *new_strings = AM_MALLOC(ap_strings);
  ap_string *ns = mk_ap_string(x, y, s, col);

  new_strings->s1 = ns;
  if (g->strings == NULL)
    new_strings->rest_strings = NULL;
  else
    new_strings->rest_strings = g->strings;
  g->strings = new_strings;
}

/* Take care to add the disc to the *end* of the list, so that
   discs added last will be drawn up last. */
/* March 99 - AWM puts it in wrong order disobeying above to get speed */
void ap_add_disc(apict *a, double x, double y, double r, int col)
{
  ap_discs *old_discs = a->discs;
  ap_discs *new_discs = AM_MALLOC(ap_discs);

  new_discs->d1 = mk_ap_disc(x, y, r, col);
  new_discs->rest_discs = NULL;
  if (old_discs == NULL)
    a->discs = new_discs;
  else
  {
    new_discs -> rest_discs = a->discs;
    a->discs = new_discs;
  }
}

/* Take care to add the pixel to the *end* of the list, so that
   pixels added last will be drawn up last. */
/* March 99 - AWM puts it in wrong order disobeying above to get speed */
void ap_add_pixel(apict *a, double x, double y, int col)
{
  ap_pixels *old_pixels = a->pixels;
  ap_pixels *new_pixels = AM_MALLOC(ap_pixels);

  new_pixels->p1 = mk_ap_pixel(x, y, col);
  new_pixels->rest_pixels = NULL;
  if (old_pixels == NULL)
    a->pixels = new_pixels;
  else
  {
    new_pixels -> rest_pixels = a->pixels;
    a->pixels = new_pixels;
  }
}

/* Take care to add the rectangle to the *end* of the list, so that
   rectangles added last will be drawn up last. */
/* March 99 - AWM puts it in wrong order disobeying above to get speed */
void ap_add_rect(apict *g, double x, double y, double w, double h, int col)
{
  ap_rects *old_rects = g->rects;
  ap_rects *new_rects = AM_MALLOC(ap_rects);

  new_rects->r1 = mk_ap_rect(x, y, w, h, col);
  new_rects->rest_rects = NULL;
  if (old_rects == NULL)
    g->rects = new_rects;
  else
  {
    new_rects -> rest_rects = g->rects;
    g->rects = new_rects;
  }
}

/* ==================================================================== */
/* Basic functions to copy lines/circles/strings/rects/discs/apicts and */
/* to clear an apict.                                                   */
/* ==================================================================== */

/* Clears the apict, freeing its subcomponents */
void ap_clear(apict *a)
{
  am_free_ap_lines(a->lines);
  am_free_ap_circles(a->circles);
  am_free_ap_strings(a->strings);
  am_free_ap_rects(a->rects);
  am_free_ap_discs(a->discs);
  am_free_ap_pixels(a->pixels);
  a->lines = NULL;
  a->circles = NULL;
  a->strings = NULL;
  a->rects = NULL;
  a->discs = NULL;
  a->pixels = NULL;
}

/* Copies the lines in a_from into a_to.  Assumes that a_to's lines
   field is NULL initially (i.e. nothing needs to be freed). */
void ap_copy_lines(apict *a_from, apict *a_to)
{
  ap_lines *ls = a_from->lines;
  ap_line *l;

  while (ls != NULL)
  {
    l = ls->l1;
    ap_add_line(a_to, l->x1, l->y1, l->x2, l->y2, l->color);
    ls = ls->rest_lines;
  }
}

/* Copies the circles in a_from into a_to.  Assumes that a_to's circles
   field is NULL initially (i.e. nothing needs to be freed). */
void ap_copy_circles(apict *a_from, apict *a_to)
{
  ap_circles *cs = a_from->circles;
  ap_circle *c;

  while (cs != NULL)
  {
    c = cs->c1;
    ap_add_circle(a_to, c->x, c->y, c->r, c->color);
    cs = cs->rest_circles;
  }
}

/* Copies the strings in a_from into a_to.  Assumes that a_to's strings
   field is NULL initially (i.e. nothing needs to be freed). */
void ap_copy_strings(apict *a_from, apict *a_to)
{
  ap_strings *ss = a_from->strings;
  ap_string *s;

  while (ss != NULL)
  {
    s = ss->s1;
    ap_add_string(a_to, s->x, s->y, s->str, s->color);
    ss = ss->rest_strings;
  }
}

/* Copies the rects in a_from into a_to.  Assumes that a_to's rects
   field is NULL initially (i.e. nothing needs to be freed).

   4 Dec 95: amended so that the rectangles are in the same order in
   the copied list.  This is an ugly revision, needed so that when
   one box is drawn after another it will really overlay the other
   box.
*/
void ap_copy_rects(apict *a_from, apict *a_to)
{
  ap_rects *rs = a_from->rects;
  ap_rect *r;
  apict *tem = mk_apict();

  while (rs != NULL)
  {
    r = rs->r1;
    ap_add_rect(tem, r->x, r->y, r->w, r->h, r->color);
    rs = rs->rest_rects;
  }
  /* tem->rects now holds a reversed list of the rectangles... */
  rs = tem->rects;
  while (rs != NULL)
  {
    r = rs->r1;
    ap_add_rect(a_to, r->x, r->y, r->w, r->h, r->color);
    rs = rs->rest_rects;
  }
  am_free_apict(tem);
}

/* Copies the discs in a_from into a_to.  Assumes that a_to's discs
   field is NULL initially (i.e. nothing needs to be freed).

   4 Dec 95: amended so that the discs are in the same order in
   the copied list.  This is an ugly revision, needed so that when
   one disc is drawn after another it will really overlay the other
   disc.
*/
void ap_copy_discs(apict *a_from, apict *a_to)
{
  ap_discs *ds = a_from->discs;
  ap_disc *d;
  apict *tem = mk_apict();

  while (ds != NULL)
  {
        d = ds->d1;
    ap_add_disc(tem, d->x, d->y, d->r, d->color);
    ds = ds->rest_discs;
  }
  /* tem->discs now holds a reversed list of the discs... */
  ds = tem->discs;
  while (ds != NULL)
  {
        d = ds->d1;
    ap_add_disc(a_to, d->x, d->y, d->r, d->color);
    ds = ds->rest_discs;
  }
  am_free_apict(tem);
}

/* Copies the pixels in a_from into a_to.  Assumes that a_to's pixels
   field is NULL initially (i.e. nothing needs to be freed).

   4 Dec 95: amended so that the pixels are in the same order in
   the copied list.  This is an ugly revision, needed so that when
   one pixel is drawn after another it will really overlay the other
   pixel.
*/
void ap_copy_pixels(apict *a_from, apict *a_to)
{
  ap_pixels *ps = a_from->pixels;
  ap_pixel *p;
  apict *tem = mk_apict();

  while (ps != NULL)
  {
    p = ps->p1;
    ap_add_pixel(tem, p->x, p->y, p->color);
    ps = ps->rest_pixels;
  }
  /* tem->pixels now holds a reversed list of the pixels... */
  ps = tem->pixels;
  while (ps != NULL)
  {
        p = ps->p1;
    ap_add_pixel(a_to, p->x, p->y, p->color);
    ps = ps->rest_pixels;
  }
  am_free_apict(tem);
}

/* Copies the contents of a_from to a_to, first clearing a_to. */
void ap_copy_from_to(apict *a_from, apict *a_to)
{
  ap_clear(a_to);
  ap_copy_lines(a_from, a_to);
  ap_copy_circles(a_from, a_to);
  ap_copy_strings(a_from, a_to);
  ap_copy_rects(a_from, a_to);
  ap_copy_discs(a_from, a_to);
  ap_copy_pixels(a_from, a_to);
  a_to->width = a_from->width;
  a_to->height = a_from->height;
}

apict *mk_copy_apict(apict *ap)
{
  apict *result = mk_apict();
  ap_copy_from_to(ap,result);
  return(result);
}

/* ==================================================================== */
/* If apict_on has been called more recently than apict_off, then any   */
/* calls to ag_line, ag_circle, etc., should save the picture elements  */
/* to the global variable Current_Apict.  This is achieved by a series  */
/* of functions called apict_line, apict_circle, etc., that mimic       */
/* screen_line, screen_circle, etc.  (See amgr.c/amgr.h.)               */
/*                                                                      */
/* When apict_off is called, the Current_Apict structure is returned,   */
/* and Current_Apict is reset to an empty apict.                        */
/* ==================================================================== */

apict *Current_Apict = NULL;

/* Sets Current_Apict to hold a new, empty apict structure.  Subsequent
   calls to apict_line, apict_circle, etc., will add elements to this
   apict.  When apict_off is next called, it will return this apict,
   and reset Current_Apict to NULL.   This allows the user to incrementally
   build up a structure representing a picture.

   It is an error to call apict_on when Current_Apict is not NULL; hence
   it is an error to call apict_on twice without calling apict_off in
   between.
*/
void apict_on()
{
  if (Current_Apict != NULL)
    my_error("apict_on called when the Current_Apict was non-NULL");
  Current_Apict = mk_apict();
}

/* This returns the structure in Current_Apict, and resets the global
   variable to NULL.  */
apict *apict_off()
{
  apict *res = Current_Apict;

  Current_Apict = NULL;
  return(res);
}

/* ==================================================================== */
/* Principal functions.  These are the functions that the programmer    */
/* should use to manipulate apicts.                                     */
/* ==================================================================== */

/* Clears the current apict */
void apict_clear()
{
  ap_clear(Current_Apict);
}

/* Adds a line from (x1, y1) to (x2, y2) to the current apict */
void apict_line(double x1, double y1, double x2, double y2)
{
  extern int Pen_Color;
  ap_add_line(Current_Apict, x1, y1, x2, y2, Pen_Color);
}

/* Adds a circle, center (x, y), radius r, to the current apict */
void apict_circle(double x, double y, double r)
{
  extern int Pen_Color;
  ap_add_circle(Current_Apict, x, y, r, Pen_Color);
}

/* Adds a string, s, that starts at coordinates (x, y), to the current
   apict */
void apict_print(char *s, double x, double y)
{
  extern int Pen_Color;
  ap_add_string(Current_Apict, x, y, s, Pen_Color);
}

/* Adds a filled rectangle to the current apict.  The bottom left of the
   rectangle is at (x, y); its width and height are w and h.  Assumption:
   w and h are positive. */
void apict_box(double x, double y, double w, double h)
{
  extern int Pen_Color;
  ap_add_rect(Current_Apict, x, y, w, h, Pen_Color);
}

/* Adds a disc, center (x, y), radius r, to the current apict */
void apict_disc(double x, double y, double r)
{
  extern int Pen_Color;
  ap_add_disc(Current_Apict, x, y, r, Pen_Color);
}

/* Draws a pixel at (x, y). */
void apict_pixel(double x, double y)
{
  extern int Pen_Color;
  ap_add_pixel(Current_Apict,x, y,Pen_Color);
}

/* Draws a dot at (x, y). */
void apict_dot(double x, double y)
{
  apict_disc(x, y, 3.0);
}

/* Returns the apict currently being built up.  */
apict *cur_apict()
{
  return(Current_Apict);
}

/* Returns TRUE iff an apict is currently being built up
   (i.e. apict_on has been called more recently than apict_off).  */
bool apict_on_p()
{
  bool res = FALSE;

  if (Current_Apict != NULL)
    res = TRUE;
  return(res);
}

/* ==================================================================== */
/* Miscellaneous functions                                              */
/* ==================================================================== */

/* Sets a to hold colored concentric circles.  */
void set_apict_circles(apict *a)
{
  int i, j;
  int simple_spectrum_color(double fract);

  ap_clear(a);
  for (i = 200; i >= 0; i -= 1)
  {
    for (j = 0; j < 5; j++)
    ap_add_circle(a, 100 + 105 * j, 450.0, 80.0 * ((1 + i) / 201.0),
                  simple_spectrum_color(((double) i)/201.0));
  }
}

/* Sets a to hold a grid of colored boxes.  */
void set_apict_boxes(apict *a)
{
  int i, j;
  int simple_spectrum_color(double fract);

  ap_clear(a);
  for (j = 0; j < 5; j++)
        for (i = 0; i < 3; i++)
          ap_add_rect(a, 20 + j * 70.0, 380 - i * 150.0, 60.0, 140.0,
                  i * 5 + j);
}

/* This function modifies the elements in ap so that all the x
   coordinates get replaced by dx + scale * x, and all the y
   coordinates get replaced by dy + scale * y.
*/
void transform_apict(apict *ap, double dx, double dy, double scale)
{
  ap_lines *lines = ap->lines;
  ap_circles *circles = ap->circles;
  ap_strings *strings = ap->strings;
  ap_rects *rects = ap->rects;
  ap_discs *discs = ap->discs;
  ap_pixels *pixels = ap->pixels;

  ap->width = (int) (ap->width * scale);
  ap->height = (int) (ap->height * scale);

  /* Transform the lines */
  while (lines != NULL)
  {
    ap_line *l = lines->l1;

    l->x1 = dx + scale * l->x1;
    l->x2 = dx + scale * l->x2;
    l->y1 = dy + scale * l->y1;
    l->y2 = dy + scale * l->y2;
    lines = lines->rest_lines;
  }

  /* Transform the circles */
  while (circles != NULL)
  {
    ap_circle *c = circles->c1;

    c->x = dx + scale * c->x;
    c->y = dy + scale * c->y;
    c->r = scale * c->r;
    circles = circles->rest_circles;
  }

  /* Transform the strings */
  while (strings != NULL)
  {
    ap_string *s = strings->s1;

    s->x = dx + scale * s->x;
    s->y = dy + scale * s->y;
    strings = strings->rest_strings;
  }

  /* Transform the rects */
  while (rects != NULL)
  {
    ap_rect *r = rects->r1;

    r->x = dx + scale * r->x;
    r->y = dy + scale * r->y;
    r->w = scale * r->w;
    r->h = scale * r->h;
    rects = rects->rest_rects;
  }

  /* Transform the discs */
  while (discs != NULL)
  {
    ap_disc *d = discs->d1;

    d->x = dx + scale * d->x;
    d->y = dy + scale * d->y;
    d->r = scale * d->r;
    discs = discs->rest_discs;
  }
  /* Transform the pixels */
  while (pixels != NULL)
  {
    ap_pixel *p = pixels->p1;

    p->x = dx + scale * p->x;
    p->y = dy + scale * p->y;
    pixels = pixels->rest_pixels;
  }
}

/* =================================================================== */

/* Code to display an apict in the Unix X-windows or Windows NT world.
   This will only have an effect if ag_on has been called earlier.  */
void render_apict(apict *g)
{
  ap_lines *l = g->lines;
  ap_line *l1;
  ap_circles *c = g->circles;
  ap_circle *c1;
  ap_strings *s = g->strings;
  ap_string *s1;
  ap_rects *r = g->rects;
  ap_rect *r1;
  ap_discs *d = g->discs;
  ap_disc *d1;
  ap_pixels *p = g->pixels;
  ap_pixel *p1;
  int old_color;

  old_color = ag_pen_color();

  while (l != NULL)
  {
    l1 = l->l1;
    if (l1->color != ag_pen_color())
      ag_set_pen_color(l1->color);
    ag_line(l1->x1, l1->y1, l1->x2, l1->y2);
    l = l->rest_lines;
  }
  while (c != NULL)
  {
    c1 = c->c1;
    if (c1->color != ag_pen_color())
      ag_set_pen_color(c1->color);
    ag_circle(c1->x, c1->y, c1->r);
    c = c->rest_circles;
  }
  while (r != NULL)
  {
    r1 = r->r1;
    if (r1->color != ag_pen_color())
      ag_set_pen_color(r1->color);
    ag_box(r1->x, r1->y, r1->x + r1->w, r1->y + r1->h);
    r = r->rest_rects;
  }
  while (d != NULL)
  {
    d1 = d->d1;
    if (d1->color != ag_pen_color())
      ag_set_pen_color(d1->color);
    ag_disc(d1->x, d1->y, d1->r);
    d = d->rest_discs;
  }
  while (p != NULL)
  {
    p1 = p->p1;
    if (p1->color != ag_pen_color())
      ag_set_pen_color(p1->color);
    ag_pixel(p1->x, p1->y);
    p = p->rest_pixels;
  }
  while (s != NULL)
  {
    s1 = s->s1;
    if (s1->color != ag_pen_color())
      ag_set_pen_color(s1->color);
    ag_print(s1->x, s1->y, s1->str);
    s = s->rest_strings;
  }
  ag_set_pen_color(old_color);
}

/***** apict_array time! ****/

#define INITIAL_APICT_SET_SIZE 10

apict_set *mk_empty_apict_set()
{
  apict_set *aps = AM_MALLOC(apict_set);
  aps -> size = 0;
  aps -> array_size = INITIAL_APICT_SET_SIZE;
  aps -> array = AM_MALLOC_ARRAY(apict_ptr,aps->array_size);
  return(aps);
}

void add_to_apict_set(apict_set *aps,apict *this_apict)
/*
     Assume apict_set is previously of size n. After this it is of size
   n+1, and the n+1'th element is a COPY of this_apict.
*/
{
  if ( aps -> size == aps -> array_size )
  {
    int new_size = 2 + 2 * aps->array_size;
    apict **new_array = AM_MALLOC_ARRAY(apict_ptr,new_size);
    int i;
    for ( i = 0 ; i < aps -> array_size ; i++ )
      new_array[i] = aps->array[i];
    AM_FREE_ARRAY(aps->array,apict_ptr,aps->array_size);
    aps -> array = new_array;
    aps -> array_size = new_size;
  }
  aps->array[aps->size] = mk_copy_apict(this_apict);
  aps->size += 1;
}

int apict_set_size(apict_set *aps)
{
  return(aps->size);
}

apict *apict_set_ref(apict_set *aps,int index)
/*
     Returns a pointer (not a copy) to the index'th element stored in
   the apict_set. Error if index < 0 or index >= size
*/
{
  apict *result;
  if ( index < 0 || index >= apict_set_size(aps) )
  {
    result = NULL;
    my_error("apict_set_ref");
  }
  else
    result = aps->array[index];
  return(result);
}

void fprintf_apict_set(FILE *s,char *m1,apict_set *aps,char *m2)
{
  if ( apict_set_size(aps) == 0 )
    fprintf(s,"%s = <apict_set with zero entries>%s",m1,m2);
  else
  {
    int i;
    for ( i = 0 ; i < apict_set_size(aps) ; i++ )
    {
      fprintf(s,"\n%s[%2d] follows (on stdout).....\n\n",m1,i);
      print_apict(apict_set_ref(aps,i));
    }
  }
}

void free_apict_set(apict_set *aps)
{
  int i;
  for ( i = 0 ; i < apict_set_size(aps) ; i++ )
    free_apict(aps->array[i]);
  AM_FREE_ARRAY(aps->array,apict_ptr,aps->array_size);
  AM_FREE(aps,apict_set);
}

apict_set *mk_copy_apict_set(apict_set *aps)
{
  apict_set *new_ar = mk_empty_apict_set();
  int i;

  for ( i = 0 ; i < apict_set_size(aps) ; i++ )
    add_to_apict_set(new_ar,apict_set_ref(aps,i));

  return(new_ar);
}


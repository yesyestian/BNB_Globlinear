/*
       File:             apict.h
       Author:           Mary Soon Lee
       Created:          16 Jan 96
       Modified:         3 Nov 97

       Description: Header file for code supporting apicts, the
                    structures used to store pictures.  Elements of
                    these pictures may be circles, lines, strings,
                    rectangles, and discs (solid circles).
*/

#ifndef APICT_H
#define APICT_H

#include "./utils/ambs.h"
#include "./utils/amma.h"

#define MAX_APICT_STR_LEN 400

/* The structure to represent a line */
typedef struct ap_line
{
  double x1;
  double x2;
  double y1;
  double y2;
  int color;
} ap_line;

typedef struct ap_lines
{
  ap_line *l1;
  struct ap_lines *rest_lines;
} ap_lines;

/* The structure to represent a circle */
typedef struct ap_circle
{
  double x;   /* The x-coordinate of the center */
  double y;   /* The y-coordinate of the center */
  double r;   /* The radius */
  int color;
} ap_circle;

typedef struct ap_circles
{
  ap_circle *c1;
  struct ap_circles *rest_circles;
} ap_circles;

/* The structure to represent a string */
typedef struct ap_string
{
  double x;   /* The x-coordinate of the start */
  double y;   /* The y-coordinate of the start */
  char str[MAX_APICT_STR_LEN + 1];
  int color;
} ap_string;

typedef struct ap_strings
{
  ap_string *s1;
  struct ap_strings *rest_strings;
} ap_strings;

/* The structure to represent a rectangle */
typedef struct ap_rect
{
  double x;   /* The x-coordinate of the bottom left */
  double y;   /* The y-coordinate of the bottom left */
  double w;   /* The width */
  double h;   /* The height */
  int color;
} ap_rect;

typedef struct ap_rects
{
  ap_rect *r1;
  struct ap_rects *rest_rects;
} ap_rects;

/* The structure to represent a disc */
typedef struct ap_disc
{
  double x;   /* The x-coordinate of the center */
  double y;   /* The y-coordinate of the center */
  double r;   /* The radius */
  int color;
} ap_disc;

typedef struct ap_discs
{
  ap_disc *d1;
  struct ap_discs *rest_discs;
} ap_discs;

/* The structure to represent a pixel */
typedef struct ap_pixel
{
  double x;   /* The x-coordinate of the center */
  double y;   /* The y-coordinate of the center */
  int color;
} ap_pixel;

typedef struct ap_pixels
{
  ap_pixel *p1;
  struct ap_pixels *rest_pixels;
} ap_pixels;

typedef struct apict
{
  int width;  /* The width of the apict.  */
  int height; /* The height of the apict.  */
  ap_lines *lines;
  ap_circles *circles;
  ap_strings *strings;
  ap_rects *rects;
  ap_discs *discs;
  ap_pixels *pixels;
} apict,*apict_ptr;

typedef struct apict_set
{
  int size;
  int array_size;
  apict **array;
} apict_set;

apict *mk_apict();
void free_apict(apict *g);
void am_free_apict(apict *g); /* Same as free_apict. retained for compat. */
apict *mk_copy_apict(apict *ap);

void print_apict(apict *g);
void ap_add_line(apict *g, double x1, double y1, double x2, double y2,
                 int col);
void ap_add_circle(apict *g, double x, double y, double r, int col);
void ap_add_string(apict *g, double x, double y, char *s, int col);
void ap_add_rect(apict *g, double x, double y, double w, double h, int col);
void ap_add_disc(apict *g, double x, double y, double r, int col);
void ap_copy_from_to(apict *a_from, apict *a_to);
void set_apict_circles(apict *a);
apict *cur_apict();
void render_apict(apict *g);
void apict_on();
apict *apict_off();

/* This function modifies the elements in ap so that all the x
   coordinates get replaced by dx + scale * x, and all the y
   coordinates get replaced by dy + scale * y.
*/
void transform_apict(apict *ap, double dx, double dy, double scale);

apict_set *mk_empty_apict_set();
void add_to_apict_set(apict_set *aps,apict *this_apict);
int apict_set_size(apict_set *aps);
apict *apict_set_ref(apict_set *aps,int index);
void fprintf_apict_set(FILE *s,char *m1,apict_set *aps,char *m2);
void free_apict_set(apict_set *aps);
apict_set *mk_copy_apict_set(apict_set *aps);

#endif /* #ifndef APICT_H */


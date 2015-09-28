/*
   File:        lingraph.c
   Author:      Andrew W. Moore
   Created:     Wed May  8 19:16:36 EDT 1996
   Description: Graphs with lines on them

   Copyright (C) 1996, A. W. Moore
*/

#include "lingraph.h"
#include "ongr.h"

laxis *mk_empty_laxis()
{
  laxis *x = AM_MALLOC(laxis);
  x -> name = NULL;
  x -> min_defined = FALSE;
  x -> min_value = -777;
  x -> max_defined = FALSE;
  x -> max_value = -777;
  return(x);
}

lingraph *mk_empty_lingraph()
{
  lingraph *x = AM_MALLOC(lingraph);
  x -> lines = mk_empty_dym_array();
  x -> color_codes = mk_ivec(0);
  x -> styles = mk_ivec(0);
  x -> labels = mk_string_array(0);
  x -> title = NULL;
  x -> string_locations = mk_dym(0,2);
  x -> strings = mk_string_array(0);

  x -> x = mk_empty_laxis();
  x -> y = mk_empty_laxis();
  return(x);
}

void free_laxis(laxis *x)
{
  if ( x->name != NULL ) free_string(x->name);
  AM_FREE(x,laxis);
}

void free_lingraph(lingraph *x)
{
  if ( x->lines != NULL ) free_dym_array(x->lines);
  if ( x->color_codes != NULL ) free_ivec(x->color_codes);
  if ( x->styles != NULL ) free_ivec(x->styles);
  if ( x->labels != NULL ) free_string_array(x->labels);
  if ( x->title != NULL ) free_string(x->title);
  if ( x->string_locations != NULL ) free_dym(x->string_locations);
  if ( x->strings != NULL ) free_string_array(x->strings);
  if ( x->x != NULL ) free_laxis(x->x);
  if ( x->y != NULL ) free_laxis(x->y);
  AM_FREE(x,lingraph);
}

laxis *mk_copy_laxis(laxis *x)
{
  laxis *nx = AM_MALLOC(laxis);
  nx->name = (x->name==NULL) ? NULL : mk_copy_string(x->name);
  nx -> min_defined = x -> min_defined;
  nx -> min_value = x -> min_value;
  nx -> max_defined = x -> max_defined;
  nx -> max_value = x -> max_value;
  return(nx);
}

lingraph *mk_copy_lingraph(lingraph *x)
{
  lingraph *nx = AM_MALLOC(lingraph);
  nx->lines = (x->lines==NULL) ? NULL : mk_copy_dym_array(x->lines);
  nx->color_codes = (x->color_codes==NULL) ? NULL : mk_copy_ivec(x->color_codes);
  nx->styles = (x->styles==NULL) ? NULL : mk_copy_ivec(x->styles);
  nx->labels = (x->labels==NULL) ? NULL : mk_copy_string_array(x->labels);
  nx->title = (x->title==NULL) ? NULL : mk_copy_string(x->title);
  nx->string_locations = (x->string_locations==NULL) ? NULL : mk_copy_dym(x->string_locations);
  nx->strings = (x->strings==NULL) ? NULL : mk_copy_string_array(x->strings);
  nx->x = (x->x==NULL) ? NULL : mk_copy_laxis(x->x);
  nx->y = (x->y==NULL) ? NULL : mk_copy_laxis(x->y);
  return(nx);
}

void fprintf_laxis(FILE *s,char *m1,laxis *x,char *m2)
{
  char buff[100];
  sprintf(buff,"%s -> name",m1);
  fprintf_string(s,buff,x->name,m2);
  sprintf(buff,"%s -> min_defined",m1);
  fprintf_bool(s,buff,x->min_defined,m2);
  sprintf(buff,"%s -> min_value",m1);
  fprintf_double(s,buff,x->min_value,m2);
  sprintf(buff,"%s -> max_defined",m1);
  fprintf_bool(s,buff,x->max_defined,m2);
  sprintf(buff,"%s -> max_value",m1);
  fprintf_double(s,buff,x->max_value,m2);
}

void fprintf_lingraph(FILE *s,char *m1,lingraph *x,char *m2)
{
  char buff[100];
  sprintf(buff,"%s -> lines",m1);
  fprintf_dym_array(s,buff,x->lines,m2);
  sprintf(buff,"%s -> color_codes",m1);
  fprintf_ivec(s,buff,x->color_codes,m2);
  sprintf(buff,"%s -> styles",m1);
  fprintf_ivec(s,buff,x->styles,m2);
  sprintf(buff,"%s -> labels",m1);
  fprintf_string_array(s,buff,x->labels,m2);
  sprintf(buff,"%s -> title",m1);
  fprintf_string(s,buff,x->title,m2);
  sprintf(buff,"%s -> string_locations",m1);
  fprintf_dym(s,buff,x->string_locations,m2);
  sprintf(buff,"%s -> strings",m1);
  fprintf_string_array(s,buff,x->strings,m2);
  sprintf(buff,"%s -> x",m1);
  fprintf_laxis(s,buff,x->x,m2);
  sprintf(buff,"%s -> y",m1);
  fprintf_laxis(s,buff,x->y,m2);
}

void plingraph(lingraph *lg)
{
  fprintf_lingraph(stdout,"lg",lg,"\n");
}

int next_color_code(ivec *color_codes)
{
  if ( ivec_size(color_codes)==0 )
    return(0);
  else
    return(1 + ivec_max(color_codes));
}

void maybe_expand_num_lines(lingraph *lg,int line_num)
{
  if ( line_num < 0 || line_num > 1000 )
    my_error("lingraph linenums should be between 0 -- 1000");
  while ( dym_array_size(lg->lines) <= line_num )
  {
    dym *new_line = mk_dym(0,2);
    add_to_dym_array(lg->lines,new_line);
    add_to_ivec(lg->color_codes,next_color_code(lg->color_codes));
    add_to_ivec(lg->styles,JOINED_STYLE);
    add_to_string_array(lg->labels,NULL);
    free_dym(new_line);
  }
}

int lingraph_num_lines(lingraph *lg)
{
  return(dym_array_size(lg->lines));
}

dym *lingraph_line_ref(lingraph *lg,int line_num)
{
  maybe_expand_num_lines(lg,line_num);
  return(dym_array_ref(lg->lines,line_num));
}

void add_to_lingraph(lingraph *lg,int line_num,double x,double y)
{
  dym *line = lingraph_line_ref(lg,line_num);
  add_row(line);
  dym_set(line,dym_rows(line)-1,0,x);
  dym_set(line,dym_rows(line)-1,1,y);
}

void set_lines_same_color(lingraph *lg,int line_num1,int line_num2)
{
  maybe_expand_num_lines(lg,line_num1);
  maybe_expand_num_lines(lg,line_num2);
  if ( line_num1 != line_num2 )
  {
    int col1 = ivec_ref(lg->color_codes,line_num1);
    int col2 = ivec_ref(lg->color_codes,line_num2);
    if ( col1 != col2 )
    {
      int i;
      ivec_set(lg->color_codes,line_num2,col1);
      if ( !is_in_ivec(lg->color_codes,col2) )
      {
        int new_max_col = next_color_code(lg->color_codes)-1;
        if ( new_max_col > col2 )
          for ( i = 0 ; i < ivec_size(lg->color_codes) ; i++ )
            if ( ivec_ref(lg->color_codes,i) == new_max_col )
              ivec_set(lg->color_codes,i,col2);
      }
    }
  }
}

void set_line_style_dotted(lingraph *lg,int line_num)
{
  maybe_expand_num_lines(lg,line_num);
  ivec_set(lg->styles,line_num,DOTTED_STYLE);
}

void set_line_style_pixel(lingraph *lg,int line_num)
{
  maybe_expand_num_lines(lg,line_num);
  ivec_set(lg->styles,line_num,PIXEL_STYLE);
}

void maybe_replace(char **x,char *label)
{
  if ( *x != NULL ) free_string(*x);
  *x = mk_copy_string(label);
}

void set_x_axis_label(lingraph *lg,char *label)
{
  maybe_replace(&(lg->x->name),label);
}

void set_y_axis_label(lingraph *lg,char *label)
{
  maybe_replace(&(lg->y->name),label);
}

void set_lingraph_title(lingraph *lg,char *label)
{
  maybe_replace(&(lg->title),label);
}

void add_string_to_lingraph(lingraph *lg,double x,double y,char *label)
{
  add_row(lg->string_locations);
  dym_set(lg->string_locations,dym_rows(lg->string_locations)-1,0,x);
  dym_set(lg->string_locations,dym_rows(lg->string_locations)-1,1,y);
  add_to_string_array(lg->strings,label);
}

void set_line_label(lingraph *lg,int line_num,char *label)
{
  maybe_expand_num_lines(lg,line_num);
  string_array_set(lg->labels,line_num,label);
}

void set_x_min(lingraph *lg,double xmin)
{
  lg->x->min_defined = TRUE;
  lg->x->min_value = xmin;
}

void set_x_max(lingraph *lg,double xmax)
{
  lg->x->max_defined = TRUE;
  lg->x->max_value = xmax;
}

void set_y_min(lingraph *lg,double ymin)
{
  lg->y->min_defined = TRUE;
  lg->y->min_value = ymin;
}

void set_y_max(lingraph *lg,double ymax)
{
  lg->y->max_defined = TRUE;
  lg->y->max_value = ymax;
}

int make_xarr_yarr_from_lingraph(lingraph *lg,int line_num,
				 double **xarr,double **yarr)
{
  dym *line = lingraph_line_ref(lg,line_num);
  int size = dym_rows(line);
  int i;

  *xarr = AM_MALLOC_ARRAY(double,size);
  *yarr = AM_MALLOC_ARRAY(double,size);

  for ( i = 0 ; i < size ; i++ )
  {
    (*xarr)[i] = dym_ref(line,i,0);
    (*yarr)[i] = dym_ref(line,i,1);
  }

  return(size);
}

#ifdef APICT_FROM_LINGRAPH_IS_ACTUALLY_NEEDED_BY_SOMEONE
apict *mk_apict_from_lingraph(lingraph *lg)
{
  apict_on();
  render_lingraph(lg);
  return apict_off();
}
#endif

void render_lingraph_in_agbox(agbox *agb,lingraph *lg)
{
  ongr on[1];
  frame fr[1];
  int i;

  fr->bottom_left.x = agb -> xlo;
  fr->bottom_left.y = agb -> ylo;
  fr->top_right.x = agb -> xhi;
  fr->top_right.y = agb -> yhi;

  if ( lg->title != NULL ) print_and_shrink(fr,lg->title);
  clear_ongr(on);

  for ( i = 0 ; i < lingraph_num_lines(lg) ; i++ )
  {
    double *xarr,*yarr;
    int size = make_xarr_yarr_from_lingraph(lg,i,&xarr,&yarr);
    if ( i == 0 )
    {
      compute_axis_limits(xarr,size,&on->x_axis);
      compute_axis_limits(yarr,size,&on->y_axis);
    }
    else
    {
      maybe_expand_axis_limits(xarr,size,&on->x_axis);
      maybe_expand_axis_limits(yarr,size,&on->y_axis);
    }
    AM_FREE_ARRAY(xarr,double,size);
    AM_FREE_ARRAY(yarr,double,size);
  }

  if ( lg->x->min_defined )
    on->x_axis.lowest_val = real_min(lg->x->min_value,on->x_axis.lowest_val);
  if ( lg->x->max_defined )
    on->x_axis.highest_val = real_max(lg->x->max_value,on->x_axis.highest_val);
  if ( lg->y->min_defined )
    on->y_axis.lowest_val = real_min(lg->y->min_value,on->y_axis.lowest_val);
  if ( lg->y->max_defined )
    on->y_axis.highest_val = real_max(lg->y->max_value,on->y_axis.highest_val);

  if ( lg->x->name != NULL ) sprintf(on->x_axis.label,"%s",lg->x->name);
  if ( lg->y->name != NULL ) sprintf(on->y_axis.label,"%s",lg->y->name);

  compute_axes_details(on,fr);

  draw_axes(fr,on);

  for ( i = 0 ; i < lingraph_num_lines(lg) ; i++ )
  {
    void plot_graphic_in_frame(frame *fr, 
                               ongr *on, 
                               double *x_arr,
                               double *y_arr,
                               int size,
                               char *mark_code,
                               double bar_width,
                               int amut_col
                              );

    double *xarr,*yarr;
    int size = make_xarr_yarr_from_lingraph(lg,i,&xarr,&yarr);
    int z = ivec_ref(lg->color_codes,i);
    int amut_col = color_code_to_ag_color(z);
    char *style = (ivec_ref(lg->styles,i)==JOINED_STYLE) ? "LN" : 
                  (ivec_ref(lg->styles,i)==PIXEL_STYLE) ? "NP" : "ND";
    char *label = string_array_ref(lg->labels,i);
    int old_col = ag_pen_color();
 
    plot_graphic_in_frame(fr,on,xarr,yarr,size,style,-1.0,amut_col);

    if ( label != NULL )
    {
      ag_set_pen_color(amut_col);
      ag_print(100.0,490.0-20.0 * i,label);
    }

    ag_set_pen_color(old_col);
    AM_FREE_ARRAY(xarr,double,size);
    AM_FREE_ARRAY(yarr,double,size);
  }

  ag_set_pen_color(AG_BLACK);

  for ( i = 0 ; i < dym_rows(lg->string_locations) ; i++ )
  {
    void ongr_print( frame *f, onpoint *p, char *string );
    char *string = string_array_ref(lg->strings,i);
    onpoint pt[1];

    pt -> x = dym_ref(lg->string_locations,i,0);
    pt -> y = dym_ref(lg->string_locations,i,1);

    ongr_print(fr,pt,string);
  }
}

void render_lingraph(lingraph *lg)
{
  agbox *agb = mk_agbox(0.0,0.0,512.0,512.0);
  render_lingraph_in_agbox(agb,lg);
  free_agbox(agb);
}

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
void add_ellipse_to_lingraph( lingraph *lg, int linenum, 
			      double x0, double y0, 
			      double a, double b, double phi, 
			      int resolution )
{
  double step = 2*PI/(double)resolution;
  double t = 0.0;
  int i;
  double sinphi = sin(phi*PI/180.0);
  double cosphi = cos(phi*PI/180.0);
  double x,y;

  for( i=0; i<=resolution; i++ )
    {
      x = a*cos(t);
      y = b*sin(t);
      add_to_lingraph( lg, linenum, x0+x*cosphi-y*sinphi, y0+x*sinphi+y*cosphi );
      t += step;
    }

  return;
}

#ifdef NEVER
void test_lingraph(int argc,char *argv[])
{
  lingraph *lg = mk_empty_lingraph();
  int num_lines = 8;
  int i;

  plingraph(lg);
  wait_for_key();
  set_lines_same_color(lg,1,3);
  plingraph(lg);
  wait_for_key();

  for ( i = 0 ; i < num_lines ; i++ )
  {
    int j;
    for ( j = 0 ; j < 60 ; j++ )
    {
      double x = (double) (i + j/4.0);
      double y = (i + 5.0) * sin(x/3.0) + (i + j/4.0 + 10.0) / 1.5;
      add_to_lingraph(lg,i,x,y);
    }
  }
  set_lines_same_color(lg,3,num_lines-1);
  set_line_style_dotted(lg,4);
  set_x_axis_label(lg,"Hairyness");
  set_y_axis_label(lg,"Aggressiveness");
  set_lingraph_title(lg,"A very sad story");
  set_x_min(lg,-2.0);

  plingraph(lg);
  wait_for_key();
  render_lingraph(lg);
  wait_for_key();
  free_lingraph(lg);
}
#endif /*NEVER*/

void add_dyv_to_lingraph(lingraph *lg,int line,dyv *dv)
{
  int i;
  for ( i = 0 ; i < dyv_size(dv) ; i++ )
    add_to_lingraph(lg,line,(double)i,dyv_ref(dv,i));
}

lingraph *mk_lingraph_from_dyv(char *title,dyv *dv)
{
  lingraph *lg = mk_empty_lingraph();
  set_lingraph_title(lg,title);
  add_dyv_to_lingraph(lg,0,dv);
  return lg;
}

void draw_dyv(char *title,dyv *dv)
{
  lingraph *lg = mk_lingraph_from_dyv(title,dv);
  render_lingraph(lg);
  free_lingraph(lg);
}



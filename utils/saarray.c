/*
File: saarray.cpp
Author: Andrew Moore


*/

#include "saarray.h"

#define INITIAL_sa_array_SIZE 10

sa_array *mk_empty_sa_array()
{
  sa_array *sa_arr = AM_MALLOC(sa_array);
  sa_arr -> size = 0;
  sa_arr -> array_size = INITIAL_sa_array_SIZE;
  sa_arr -> array = AM_MALLOC_ARRAY(string_array_ptr,sa_arr->array_size);
  return(sa_arr);
}

void add_to_sa_array(sa_array *sa_arr,string_array *sa)
/*
     Assume sa_arr is previously of size n. After this it is of size
   n+1, and the n+1'th element is a COPY of sa.
*/
{
  if ( sa_arr -> size == sa_arr -> array_size )
  {
    int new_size = 2 + 2 * sa_arr->array_size;
    string_array **new_array = AM_MALLOC_ARRAY(string_array_ptr,new_size);
    int i;
    for ( i = 0 ; i < sa_arr -> array_size ; i++ )
      new_array[i] = sa_arr->array[i];
    AM_FREE_ARRAY(sa_arr->array,string_array_ptr,sa_arr->array_size);
    sa_arr -> array = new_array;
    sa_arr -> array_size = new_size;
  }
  sa_arr->array[sa_arr->size] = mk_copy_string_array(sa);
  sa_arr->size += 1;
}

int sa_array_size(sa_array *sa_arr)
{
  return(sa_arr->size);
}

string_array *sa_array_ref(sa_array *sa_arr,int index)
/*
     Returns a pointer (not a copy) to the index'th element stored in
   the sa_arr. Error if index < 0 or index >= size
*/
{
  string_array *result;
  if ( index < 0 || index >= sa_array_size(sa_arr) )
  {
    result = NULL;
    my_error("sa_array_ref");
  }
  else
    result = sa_arr->array[index];
  return(result);
}
  
void fprintf_sa_array(FILE *s,char *m1,sa_array *sa_arr,char *m2)
{
  if ( sa_array_size(sa_arr) == 0 )
    fprintf(s,"%s = <sa_arr with zero ensaies>%s",m1,m2);
  else
  {
    int i;
    for ( i = 0 ; i < sa_array_size(sa_arr) ; i++ )
    {
      char buff[100];
      sprintf(buff,"%s[%2d]",m1,i);
      fprintf_string_array(s,buff,sa_array_ref(sa_arr,i),m2);
    }
  }
}

void psa_array(sa_array *saarr)
{
  fprintf_sa_array(stdout,"saarr",saarr,"\n");
}

void free_sa_array(sa_array *sa_arr)
{
  int i;
  for ( i = 0 ; i < sa_array_size(sa_arr) ; i++ )
    free_string_array(sa_arr->array[i]);
  AM_FREE_ARRAY(sa_arr->array,string_array_ptr,sa_arr->array_size);
  AM_FREE(sa_arr,sa_array);
}

sa_array *mk_copy_sa_array(sa_array *sa_arr)
{
  sa_array *new_ar = mk_empty_sa_array();
  int i;

  for ( i = 0 ; i < sa_array_size(sa_arr) ; i++ )
    add_to_sa_array(new_ar,sa_array_ref(sa_arr,i));

  return(new_ar);
}

#ifdef EXTINCT
int tset_num_bar_graphs(tset *ts)
{
  return(dym_array_size(ts->values));
}

int tset_num_bars(tset *ts)
{
  return(string_array_size(ts->bar_names));
}

int tset_num_cols(tset *ts)
{
  return(string_array_size(ts->col_names));
}

int tset_num_lines_in_title(tset *ts,int bar_graph_index)
{
  return(string_array_size(sa_array_ref(ts->titles,bar_graph_index)));
}

char *tset_title_line(tset *ts,int bg_index,int line_num)
{
  return(string_array_ref(sa_array_ref(ts->titles,bg_index),line_num));
}

char *tset_col_name(tset *ts,int col_index)
{
  return(string_array_ref(ts->col_names,col_index));
}

char *tset_bar_name(tset *ts,int bar_index)
{
  return(string_array_ref(ts->bar_names,bar_index));
}

double tset_value(tset *ts,int bg_index,int bar_index,int col_index)
{
  return(dym_ref(dym_array_ref(ts->values,bg_index),bar_index,col_index));
}

double tset_max_value_in_col(tset *ts,int col_index)
{
  double result = 0.0;
  int bg_index;
  int bar_index;
  for ( bg_index = 0 ; bg_index < tset_num_bar_graphs(ts) ; bg_index++ )
    for ( bar_index = 0 ; bar_index < tset_num_bars(ts) ; bar_index++ )
    {
      double value = tset_value(ts,bg_index,bar_index,col_index);
      if ( (bg_index==0 && bar_index==0) || value > result )
        result = value;
    }
  return(result);
}

tset *mk_tset(int num_cols)
{
  tset *ts = AM_MALLOC(tset);
  int col_index;
  ts -> col_names = mk_string_array(0);

  for ( col_index = 0 ; col_index < num_cols ; col_index++ )
  {
    char buff[100];
    sprintf(buff,"col%d",col_index);
    add_to_string_array(ts->col_names,buff);
  }

  ts -> titles = mk_empty_sa_array();
  ts -> values = mk_empty_dym_array();
  ts -> bar_names = mk_string_array(0);

  return(ts);
}

void free_tset(tset *ts)
{
  free_sa_array(ts->titles);
  free_dym_array(ts->values);
  free_string_array(ts->bar_names);
  free_string_array(ts->col_names);
  AM_FREE(ts,tset);
}

void tset_add_bar_graph(tset *ts,char *title)
{
  string_array *new_title = mk_string_array_1(title);
  dym *new_value = mk_zero_dym(tset_num_bars(ts),tset_num_cols(ts));
  add_to_sa_array(ts->titles,new_title);
  add_to_dym_array(ts->values,new_value);
  free_string_array(new_title);
  free_dym(new_value);
}

void tset_add_bar_graph_title(tset *ts,int bg_index,char *new_title_line)
{
  string_array *sa = sa_array_ref(ts->titles,bg_index);
  add_to_string_array(sa,new_title_line);
}

void tset_add_bar(tset *ts,char *name)
{
  int bg_index;
  dyv *zeroes = mk_zero_dyv(tset_num_cols(ts));
  add_to_string_array(ts->bar_names,name);
  for ( bg_index = 0 ; bg_index < tset_num_bar_graphs(ts) ; bg_index++ )
  {
    dym *value = dym_array_ref(ts->values,bg_index);
    add_row(value);
    copy_dyv_to_dym_row(zeroes,value,dym_rows(value)-1);
  }
  free_dyv(zeroes);
}

void tset_set_col_name(tset *ts,int col_index,char *name)
{
  string_array_set(ts->col_names,col_index,name);
}

void tset_set_value(tset *ts,int bg_index,int bar_index,int col_index,double val)
{
  dym_set(dym_array_ref(ts->values,bg_index),bar_index,col_index,val);
}

#define DASHES_STRING "------"

string_matrix *mk_string_matrix_from_tset(tset *ts)
{
  int max_bg = tset_num_bar_graphs(ts);
  int max_bar = tset_num_bars(ts);
  int max_col = tset_num_cols(ts);
  int i,j,k;

  string_matrix *sm = mk_string_matrix(1+(1+max_bar) * max_bg,3+max_col);

  string_matrix_set(sm,0,1,"|");
  for ( i = 0 ; i < max_col ; i++ )
    string_matrix_set(sm,0,3+i,tset_col_name(ts,i));
  for ( i = 0 ; i < max_bg ; i++ )
  {
    int row = 1 + (1 + max_bar) * i;
    for ( k = 0 ; k < 3+max_col ; k++ )
      string_matrix_set(sm,row,k,DASHES_STRING);
    string_matrix_set(sm,row,1,"|");

    for ( j = 0 ; j < tset_num_lines_in_title(ts,i) ; j++ )
      string_matrix_set(sm,row+1+j,0,tset_title_line(ts,i,j));

    for ( j = 0 ; j < max_bar ; j++ )
    {
      string_matrix_set(sm,row+1+j,1,"|");
      string_matrix_set(sm,row+1+j,2,tset_bar_name(ts,j));
      for ( k = 0 ; k < max_col ; k++ )
        string_matrix_real_set(sm,row+1+j,3+k,tset_value(ts,i,j,k));
    }
  }

  return(sm);
}

expo *mk_expo_from_tset(tset *ts)
{
  string_matrix *sm = mk_string_matrix_from_tset(ts);
  expo *ex = mk_expo_from_string_matrix(sm);
  free_string_matrix(sm);
  return(ex);
}

void explain_tset(FILE *s,tset *ts)
{
  expo *ex = mk_expo_from_tset(ts);
  expo_render_plain(s,ex);
  free_expo(ex);
}

char *mk_title_from_sa(tset *ts,string_array *sa)
{
  string_array *tit = mk_string_array(0);
  char *result;
  int i;
  for ( i = 0 ; i < string_array_size(sa) - tset_num_cols(ts) - 2 ; i++ )
    add_to_string_array(tit,string_array_ref(sa,i));
  result = mk_string_from_string_array(tit);
  free_string_array(tit);
  if ( result == NULL ) result = mk_copy_string("");
  return(result);
}

void add_values_from_sa(tset *ts,int bar_index,string_array *sa,char **rem)
{
  if ( *rem == NULL )
  {
    int base = string_array_size(sa) - tset_num_cols(ts);
    int bar_name_index = base-1;
    if ( bar_name_index < 0 )
      *rem = mk_copy_string("Not enough values on the line");
    else
    {
      int col_index;
      char *bar_name = string_array_ref(sa,bar_name_index);
      if ( tset_num_bars(ts) <= bar_index )
        tset_add_bar(ts,bar_name);
      else if ( !caseless_eq_string(bar_name,tset_bar_name(ts,bar_index)) )
        *rem = mk_copy_string("Bar name doesn't match");

      for ( col_index = 0 ; *rem == NULL && col_index < tset_num_cols(ts) ; col_index++ )
      {
        char *vstring = string_array_ref(sa,base+col_index);
        if ( !is_a_number(vstring) )
          *rem = mk_copy_string("Non-numerical value");
        else
        {
          int bg_index = tset_num_bar_graphs(ts)-1;
          tset_set_value(ts,bg_index,bar_index,col_index,atof(vstring));
        }
      }
    }
  }
}

/* Find next line that is non blank and starts with non comment.
   If file ends return NULL 

  */
string_array *mk_interesting_line(FILE *s,int *line)
{
  string_array *sa = NULL;
  bool file_end = FALSE;
  while ( sa == NULL && !file_end )
  {
    sa = mk_string_array_from_line(s);
    *line += 1;
    if ( sa == NULL )
      file_end = TRUE;
    else if ( string_array_size(sa) == 0 ||
              string_array_ref(sa,0)[0] == '#' )
    {
      free_string_array(sa);
      sa = NULL;
    }
  }

  return(sa);
}

/*
   Does nothing if enter with *rem non null.
   If file ends immediately, does nothing, returns FALSE.
   Else returns TRUE.
   Sets *rem non null (with message) if error.
   */
bool can_add_bar_graph_from_stream(FILE *s,tset *ts,char **rem,int *line)
{
  bool result;

  if ( *rem == NULL )
  {
    string_array *sa = mk_interesting_line(s,line);

    if ( sa == NULL )
      result = FALSE;
    else
    {
      result = TRUE;

      if ( *rem == NULL && string_array_size(sa) == 0 )
        *rem = mk_copy_string("Unexpected blank line");

      if ( *rem == NULL )
      {
        int i;
        char *title = mk_title_from_sa(ts,sa);
        bool finished = FALSE;

        tset_add_bar_graph(ts,title);
        add_values_from_sa(ts,0,sa,rem);
        free_string_array(sa);
        free_string(title);
        sa = NULL;
        title = NULL;
 
        for ( i = 1 ; !finished && *rem == NULL ; i++ )
        {
         sa = mk_interesting_line(s,line);

          if ( sa == NULL || (string_array_size(sa) > 0 && eq_string(string_array_ref(sa,0),DASHES_STRING)) )
            finished = TRUE;
          else
          {
            title = mk_title_from_sa(ts,sa);
            if ( title[0] != '\0' ) tset_add_bar_graph_title(ts,tset_num_bar_graphs(ts)-1,title);
            add_values_from_sa(ts,i,sa,rem);
            free_string(title);
            title = NULL;
          }

          if ( sa != NULL ) free_string_array(sa);
          sa = NULL;
        }
      }
    }
  }
  return(result);
}

/*
   Does nothing if enter with *rem non null.
   Sets *rem non null (with message) if error.
   */
tset *mk_tset_from_stream(FILE *s,char **rem,int *line)
{
  tset *ts = NULL;
  if ( *rem == NULL )
  {
    string_array *sa = mk_interesting_line(s,line);
    bool finished = FALSE;

    if ( sa == NULL || string_array_size(sa) < 2 )
      *rem = mk_copy_string("First line should be | <name> [<more names...>]");
    else
    {
      int i;
      ts = mk_tset(string_array_size(sa) - 1);
      for ( i = 0 ; i < string_array_size(sa) - 1 ; i++ )
        tset_set_col_name(ts,i,string_array_ref(sa,i+1));
    }

    if ( sa != NULL ) free_string_array(sa);
    sa = NULL;

    if ( *rem == NULL )
    {
      sa = mk_interesting_line(s,line);
      if ( sa == NULL || string_array_size(sa) == 0 || !eq_string(string_array_ref(sa,0),DASHES_STRING) )
        *rem = mk_copy_string("Second line of bar graph should be dashes");
      if ( sa != NULL ) free_string_array(sa);
      sa = NULL;
    }

    while ( *rem == NULL && !finished )
      finished = !can_add_bar_graph_from_stream(s,ts,rem,line);
  }

  if ( *rem != NULL && ts != NULL )
  {
    free_tset(ts);
    ts = NULL;
  }

  return(ts);
}

tset *mk_tset_from_filename(char *filename)
{
  char *em = NULL;
  int line = 0;
  FILE *s = safe_fopen(filename,"r");
  tset *ts = mk_tset_from_stream(s,&em,&line);
  if ( ts == NULL )
  {
    if ( em == NULL ) my_error("One of ts and em should be non null");
    fprintf(stderr,"Error reading tset from file %s line %d:\n%s\n",filename,line,em);
    my_error("mk_tset_from_filename");
  }
  else if ( em != NULL ) my_error("One of ts and em should be null");

  fclose(s);

  return(ts);
}

void load_main(int argc,char *argv[])
{
  char *filename = (argc < 3) ? "tset.txt" : argv[2];
  tset *ts = mk_tset_from_filename(filename);
  explain_tset(stdout,ts);
  free_tset(ts);
}

#define ytop 512.0
#define xright 512.0
#define string_indent 5.0

static int bar_colors[] = {
    AG_BLUE,
    AG_RED,
    AG_BLACK,
    AG_DARKGREEN,
    AG_MAGENTA,
    AG_BLUEGREEN,
    AG_CYAN,
    AG_PURPLE,
    AG_GREEN,
    AG_DARKBLUE,
    AG_OLIVE,
    AG_YELLOW,
    AG_DARKRED};

#define bar_colors_size (sizeof(bar_colors) / sizeof(int))

/* This draws the title for the bar_graph_num-th bargraph in ts
   into the rectangle bounded by xlow, ylow, xhigh, yhigh.

   xhigh is ignored by the current implementation, but could
   be used if we had the ability to scale the fonts to fit
   the available space.
*/
void draw_bg_title(double xlow, double ylow, double xhigh, double yhigh,
		   tset *ts, int bar_graph_num)
{
  int line_num;
  int num_lines = tset_num_lines_in_title(ts, bar_graph_num);
  double cury, height_per_line;

  if (num_lines > 0)
  {
    height_per_line = (yhigh - ylow) / ((double) num_lines);
    cury = yhigh - height_per_line / 2.0;

    for (line_num = 0; line_num < num_lines; line_num++)
    {
      ag_print(xlow + string_indent, cury, 
	       tset_title_line(ts, bar_graph_num, line_num));
      cury -= height_per_line;
    }
  }
}

/* This draws the bar-names, bar lengths, and column 2 entries
   for the bar_graph_num-th bar graph in ts, and does so at
   the specified position.
*/
void draw_bg_bar_info(double barname_xleft, double histogram_xleft, 
		      double col2_xleft, double ylow, double yhigh,
		      tset *ts, int bar_graph_num)
{
  int bar_num;
  int num_bars = tset_num_bars(ts);
  double tset_max_barvalue = tset_max_value_in_col(ts, 0);
  double cur_bartop, cur_barbottom, cur_barmiddle, height_per_bar;
  int old_pen_color = ag_pen_color();
  double width_for_col1_num = 50.0;
  double width_for_bars = col2_xleft - (histogram_xleft + width_for_col1_num);
  double cur_bar_len, cur_bar_value;
  char buff[100];
  double reduce_bar_height_by_twice_this = 2.0;

  if (num_bars > 0)
  {
    if (num_bars > bar_colors_size)
      my_error("draw_bg_bar_info: would run out of colors");

    height_per_bar = (yhigh - ylow) / ((double) num_bars);
    cur_bartop = yhigh;
    cur_barbottom = yhigh - height_per_bar;

    for (bar_num = 0; bar_num < num_bars; bar_num++)
    {
      ag_set_pen_color(bar_colors[bar_num]);
      cur_barmiddle = (cur_bartop + cur_barbottom) / 2.0;

      /* Draw the bar name */
      ag_print(barname_xleft + string_indent, 
	       cur_barmiddle,
	       tset_bar_name(ts, bar_num));

      /* Draw the bar length. */
      cur_bar_value = tset_value(ts, bar_graph_num, bar_num, 0);
      if (cur_bar_value < 0)
	my_error("draw_bg_bar_info: bar values must be >= 0.");
      cur_bar_len = width_for_bars * (cur_bar_value / tset_max_barvalue);
      ag_box(histogram_xleft, 
	     cur_barbottom + reduce_bar_height_by_twice_this,
	     histogram_xleft + cur_bar_len, 
	     cur_bartop - reduce_bar_height_by_twice_this);
      
      /* Draw the bar value. */
      sprintf(buff, "%5.2f", cur_bar_value);
      ag_print(histogram_xleft + cur_bar_len + string_indent,
	       cur_barmiddle,
	       buff);	       

      /* Draw the column2 value. */
      sprintf(buff, "%5.2f", tset_value(ts, bar_graph_num, bar_num, 1));
      ag_print(col2_xleft + string_indent,
	       cur_barmiddle,
	       buff);	       

      cur_bartop -= height_per_bar;
      cur_barbottom -= height_per_bar;
    }
    ag_set_pen_color(old_pen_color);
  }
}

/* This draws the bar_graph_num-th bar graph in tset at the
   specified position. */
void draw_bg(double title_xleft, double barname_xleft,
	     double histogram_xleft, double col2_xleft,
	     double ylow, double yhigh,
	     tset *ts,
	     int bar_graph_num)
{
  if (bar_graph_num < 0 || bar_graph_num >= tset_num_bar_graphs(ts))
    my_error("draw_bg_title: bar_graph_num out of range.");

  draw_bg_title(title_xleft, ylow, barname_xleft, yhigh, 
		ts, bar_graph_num);
  draw_bg_bar_info(barname_xleft, histogram_xleft, 
		   col2_xleft, ylow, yhigh,
		   ts, bar_graph_num);  
}

/* 
   This use the amgr utilities to draw histograms of the
   first column of data in ts for each of the bar graphs.  
   The data for each bar graph is shown in turn (one above
   the other).  For a particular bargraph we see:

Title for | bar name 1 | bar the length of col 1 entry for bar 1 | col2 value
bar graph | bar name 2 | bar the length of col 1 entry for bar 2 | col2 value
          |    ...     |   ...                                   |  ...

   Preconditions:
     tset must have exactly two columns.
     drawing will be crappy if tset has nore than about 4 bar graphs.

   This function added by Mary on 19 Jan 97.
*/
void draw_tset(tset *ts, char *ps_filename)
{
  double column_title_height = 30.0;
  double ybottom_column_title = ytop - column_title_height;
  double width_for_graph_titles = 150.0;
  double width_for_bar_titles = 70.0;
  double width_for_col2_values = 50.0;
  double xleft_histogram = width_for_graph_titles + width_for_bar_titles;
  double xleft_col2 = xright - width_for_col2_values; 
  int bg;
  int num_bar_graphs = tset_num_bar_graphs(ts);
  int num_bars_per_graph = tset_num_bars(ts);
  double height_per_bar_graph;
  double cury;

  /* Error checks */
  if (tset_num_cols(ts) != 2)
    my_error("draw_tset: needs ts to have exactly two columns");
  if (num_bar_graphs == 0 || num_bars_per_graph == 0)
    my_error("draw_tset: num_bar_graphs and num_bars_per_graph must both be greater than 0");

  ag_on(ps_filename);
 
  /* Draw the column titles, and a line under them */
  ag_set_pen_color(AG_BLUE);
  ag_print(xleft_histogram + string_indent, ytop - column_title_height / 2.0,
           tset_col_name(ts, 0));  
  ag_print(xleft_col2 + string_indent, ytop - column_title_height / 2.0,
           tset_col_name(ts, 1));  
  ag_line(0.0, ybottom_column_title, xright, ybottom_column_title);

  /* Draw vertical lines dividing the various areas in the graph */
  ag_line(width_for_graph_titles, 0.0,
	  width_for_graph_titles, ybottom_column_title);
  ag_line(xleft_histogram, 0.0, xleft_histogram, ytop);
  ag_line(xleft_col2, 0.0, xleft_col2, ytop);

  /* Draw the information for each bar graph in turn, plus the horizontal lines
     between bar graphs. */
  height_per_bar_graph = (ytop - column_title_height) / (double) num_bar_graphs;
  cury = ybottom_column_title;
  for (bg = 0; bg < num_bar_graphs; bg++)
  {
    draw_bg(0.0, width_for_graph_titles, xleft_histogram, xleft_col2,
            cury - height_per_bar_graph, cury,
	    ts,
	    bg);
    cury -= height_per_bar_graph;
    if (bg != num_bar_graphs - 1)
      ag_line(0.0, cury, xright, cury);
  }

  printf("Hit any key to stop\n");
  wait_for_key();
  ag_off();
}

#endif

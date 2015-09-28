/* 
   File:         motutils.cpp
   Author:       Andrew W. Moore
   Created:      Jan 13th, 1998
   Description:  Simple utilities for the motcon project that should
                 be in the standard Auton library but aren't.
*/

#include "motutils.h"

/* Simply prints "<s1>-><s2> = <x>\n" to file s. Conforms to the 
   AM-FPRINT-STANDARD described at the top of motutils.h */
void fprint_int(FILE *s,char *s1,char *s2,int x)
{
  char prefix[1000];
  sprintf(prefix,"%s->%s",s1,s2);
  fprintf_int(s,prefix,x,"\n");
}

/* Simply prints "<s1>-><s2> = <x>\n" to file s. Conforms to the 
   AM-FPRINT-STANDARD described at the top of motutils.h */
void fprint_bool(FILE *s,char *s1,char *s2,bool x)
{
  char prefix[1000];
  sprintf(prefix,"%s->%s",s1,s2);
  fprintf_bool(s,prefix,x,"\n");
}

/* Simply prints "<s1>-><s2> = <x>\n" to file s. Conforms to the 
   AM-FPRINT-STANDARD described at the top of motutils.h */
void fprint_double(FILE *s,char *s1,char *s2,double x)
{
  char prefix[1000];
  sprintf(prefix,"%s->%s",s1,s2);
  fprintf_double(s,prefix,x,"\n");
}

/* Simply prints "<s1>-><s2> = <x>\n" to file s. Conforms to the 
   AM-FPRINT-STANDARD described at the top of motutils.h */
void fprint_string(FILE *s,char *s1,char *s2,char *x)
{
  char prefix[1000];
  sprintf(prefix,"%s->%s",s1,s2);
  fprintf_string(s,prefix,x,"\n");
}

/* Simply prints "<s1>-><s2> = <x>\n" to file s. Conforms to the 
   AM-FPRINT-STANDARD described at the top of motutils.h */
void fprint_dyv(FILE *s,char *s1,char *s2,dyv *x)
{
  char prefix[1000];
  sprintf(prefix,"%s->%s",s1,s2);
  fprintf_dyv(s,prefix,x,"\n");
}

/* Simply prints "<s1>-><s2> = <x>\n" to file s. Conforms to the 
   AM-FPRINT-STANDARD described at the top of motutils.h */
void fprint_dym(FILE *s,char *s1,char *s2,dym *x)
{
  char prefix[1000];
  sprintf(prefix,"%s->%s",s1,s2);
  fprintf_dym(s,prefix,x,"\n");
}

/* Simply prints "<s1>-><s2> = <x>\n" to file s. Conforms to the 
   AM-FPRINT-STANDARD described at the top of motutils.h */
void fprint_ivec(FILE *s,char *s1,char *s2,ivec *x)
{
  char prefix[1000];
  sprintf(prefix,"%s->%s",s1,s2);
  fprintf_ivec(s,prefix,x,"\n");
}

/* Simply prints "<s1>-><s2> = <x>\n" to file s. Conforms to the 
   AM-FPRINT-STANDARD described at the top of motutils.h */
void fprint_string_array(FILE *s,char *s1,char *s2,string_array *x)
{
  char prefix[1000];
  sprintf(prefix,"%s->%s",s1,s2);
  fprintf_string_array(s,prefix,x,"\n");
}

/* Parses the string, which should be the ASCII form of a set
   of space separated numbers, eg "0.8 99 -5e4 876" (which would
   generate a dyv of size 4). Raises a my_error if cannot be
   parsed in this way. String must be null (\0) terminated. */
dyv *mk_dyv_from_string_simple(char *s)
{
  dyv *x = mk_dyv_from_string(s,NULL);
  if ( x == NULL )
  {
    printf("Can't parse this as a vector: %s\n",s);
    my_error("mk_dyv_from_string_simple");
  }
  return x;
}

/* Raises a my_error if it is ever called with the argument value FALSE.
   Efficient versions will implement am_assert as a macro that expands
   to nothing. */
void am_assert(bool assertion)
{
  if ( !assertion )
    my_error("am_assert: assertion failed");
}

/* Increases the number of rows in dm by 1 and copies the contents
   of dv into that new bottom row. Note this means that 
     dyv_size(dv) must equal dym_cols(dm)
   in order that dv will fit (if this is violated there'll be a my_error */
void add_dyv_to_dym(dym *dm,dyv *dv)
{
  add_row(dm);
  copy_dyv_to_dym_row(dv,dm,dym_rows(dm)-1);
}


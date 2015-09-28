/*
   File:         motutils.cpp
   Author:       Andrew W. Moore
   Created:      Jan 13th, 1998
   Description:  Simple utilities for the motcon project that should
                 be in the standard Auton library but aren't.
*/

#ifndef MOTUTILS_H
#define MOTUTILS_H

#include "./utils/amdmex.h"

/* The AM-FPRINT-STANDARD defined.

   Most of the basic types, and all the user defined types in the
   motcon project will have a print function associated with them.
   This print function is always of the form

     fprint_<type>(FILE *s,char *string1,char *string2,<type> *x)

   If the thing is a simple atomic object (like an int or a bool)
   this will simply print

       <string1>-><string2> = <value>

   For instance

      fprint_int(stdout,"matrix","size",7)

   would send

      matrix->size = 7

   to stdout.

  If the thing a compound structure, it will send all fields
  of the structure to the stream in turn.

   The use of two prefix strings (string1 and string2) may look
   weird at first. But it makes it pretty easy to write fprint_
   functions for new data structures easily. See motcon/task.c
   functions fprint_hyper_rect and fprint_task for examples.
*/

void fprint_int(FILE *s,char *s1,char *s2,int x);
void fprint_bool(FILE *s,char *s1,char *s2,bool x);
void fprint_double(FILE *s,char *s1,char *s2,double x);
void fprint_string(FILE *s,char *s1,char *s2,char *x);
void fprint_dyv(FILE *s,char *s1,char *s2,dyv *x);
void fprint_dym(FILE *s,char *s1,char *s2,dym *x);
void fprint_ivec(FILE *s,char *s1,char *s2,ivec *x);
void fprint_string_array(FILE *s,char *s1,char *s2,string_array *x);


/* Parses the string, which should be the ASCII form of a set
   of space separated numbers, eg "0.8 99 -5e4 876" (which would
   generate a dyv of size 4). Raises a my_error if cannot be
   parsed in this way. String must be null (\0) terminated. */
dyv *mk_dyv_from_string_simple(char *s);

/*Parses a space-separated string into its components, and puts them into an
  integer vector.*/
ivec *mk_ivec_from_string(char *s);


/* Raises a my_error if it is ever called with the argument value FALSE.
   Efficient versions will implement am_assert as a macro that expands
   to nothing. */
void am_assert(bool assertion);

/* Increases the number of rows in dm by 1 and copies the contents
   of dv into that new bottom row. Note this means that
     dyv_size(dv) must equal dym_cols(dm)
   in order that dv will fit (if this is violated there'll be a my_error */
void add_dyv_to_dym(dym *dm,dyv *dv);


#endif /* MOTUTILS_H */











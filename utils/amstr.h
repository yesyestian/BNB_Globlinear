
/*
   File:         amstr.h
   Author:       Andrew Moore
   Created:      June 25 1997 from amma.h by Frank Dellaert
   Description:  String related functions.

   Copyright (C) Andrew W. Moore, 1992
*/

#ifndef AMSTR_H
#define AMSTR_H

#include "amma.h"

extern char *make_copy_string(const char *s);

extern char *mk_copy_string(char *s);
char *mk_copy_quoted_string(char *s);
char *mk_copy_to_quoted_string(char *s);
char *mk_string_extension(char *s, char *ext);

char *mk_input_string(char *message);

char *mk_downcase_string(char *s);
char *mk_upcase_string(char *s);
bool caseless_eq_string(char *s1,char *s2);

char *mk_downcase_string_with_length(char *s,int n);
char *mk_upcase_string_with_length(char *s,int n);
bool caseless_eq_string_with_length(char *s1,char *s2,int n);

char *mk_printf(char *fmt, ...);

int count_occurences_of_char_in_string(char *s,char c);
int find_char_in_string(char *s,char c);

/* Works out the length of the non null string s, and am_free's it */
/* The following functions are synonymous (i.e. do the same thing).
   free_string is the preferred choice because it is consistent with naming
   conventions of all other free functions in the Auton libraries */
void free_string(char *s);
void am_free_string(char *s);
char *mk_concat_strings(char *s1, char *s2);

/* pattern is a string with *s * can represent 0 or more characters.
   returns true if pattern is in the string

  example
  pattern   string  return
  *x        hex     TRUE
  *xy       hex     FALSE
  an*w      andrew  TRUE
  *n*w      andrew  TRUE
*/

bool string_pattern_matches(char *pattern, char *string);

#endif /* AMSTR_H */


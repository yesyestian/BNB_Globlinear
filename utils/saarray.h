/*
File: tset.h
Author: Andrew Moore


*/

#ifndef TSET_H
#define TSET_H

#include "amdmex.h"

typedef struct sa_array
{
  int size;
  int array_size;
  string_array **array;
} sa_array;

sa_array *mk_empty_sa_array();

void add_to_sa_array(sa_array *sa_arr,string_array *sa);
int sa_array_size(sa_array *sa_arr);
string_array *sa_array_ref(sa_array *sa_arr,int index);
void fprintf_sa_array(FILE *s,char *m1,sa_array *sa_arr,char *m2);
void psa_array(sa_array *saarr);
void free_sa_array(sa_array *sa_arr);
sa_array *mk_copy_sa_array(sa_array *sa_arr);

#endif

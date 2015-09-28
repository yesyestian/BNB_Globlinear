/*
   File:        lingraphs.h
   Author:      Andrew W. Moore
   Created:     Sat Jun 26 03:47:02 EDT 1999
   Description: Header for An array of lingraphs

   Copyright 1999, Schenley Park Research
*/


#ifndef lingraphS_H
#define lingraphS_H
#include "./utils/genarray.h"
#include "lingraph.h"

typedef struct lingraphs
{
  generic_array *genarr;
} lingraphs;

lingraphs *mk_empty_lingraphs();
void add_to_lingraphs(lingraphs *array,lingraph *element);

/* Does NOT copy in the element, merely its pointer, so after calling, the
   user must NOT access or directly free element ever again */
void add_pointer_to_lingraphs(lingraphs *array,lingraph *element);

int lingraphs_size(lingraphs *array);
lingraph *lingraphs_ref(lingraphs *array,int index);
void fprintf_lingraphs(FILE *s,char *m1,lingraphs *array,char *m2);
void plingraphs(lingraphs *array);
void free_lingraphs(lingraphs *array);
lingraphs *mk_copy_lingraphs(lingraphs *array);

void render_lingraphs(lingraphs *lgs);

#endif /* #ifndef lingraphS_H */

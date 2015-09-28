/*
   File:        lingraphs.c
   Author:      Andrew W. Moore
   Created:     Sat Jun 26 03:47:02 EDT 1999
   Description: An array of lingraphs

   Copyright 1999, Schenley Park Research
*/

#include "lingraphs.h"

void void_free_lingraph(void *data)
{
  free_lingraph((lingraph *)data);
}

void *void_mk_copy_lingraph(void *data)
{
  return (void *) mk_copy_lingraph((lingraph *)data);
}

void void_fprintf_lingraph(FILE *s,char *m1,void *data,char *m2)
{
  fprintf_lingraph(s,m1,(lingraph *) data,m2);
}

lingraphs *mk_empty_lingraphs()
{
  lingraphs *array = AM_MALLOC(lingraphs);
  array -> genarr = mk_empty_generic_array(void_free_lingraph,
					   void_mk_copy_lingraph,
					   void_fprintf_lingraph);
  return array;
}

void add_to_lingraphs(lingraphs *array,lingraph *element)
{
  add_to_generic_array(array->genarr,(void *)element);
}

void add_pointer_to_lingraphs(lingraphs *array,lingraph *element)
{
  add_pointer_to_generic_array(array->genarr,(void *)element);
}

int lingraphs_size(lingraphs *array)
{
  return(generic_array_size(array->genarr));
}

lingraph *lingraphs_ref(lingraphs *array,int index)
{
  return (lingraph *) generic_array_ref(array->genarr,index);
}
  
void fprintf_lingraphs(FILE *s,char *m1,lingraphs *array,char *m2)
{
  fprintf_generic_array(s,m1,array->genarr,m2);
}

void plingraphs(lingraphs *array)
{
  fprintf_lingraphs(stdout,"lingraphs",array,"\n");
}

void free_lingraphs(lingraphs *array)
{
  free_generic_array(array->genarr);
  AM_FREE(array,lingraphs);
}

lingraphs *mk_copy_lingraphs(lingraphs *array)
{
  lingraphs *new_lingraphs = AM_MALLOC(lingraphs);
  new_lingraphs -> genarr = mk_copy_generic_array(array->genarr);
  return new_lingraphs;
}

void render_lingraphs(lingraphs *lgs)
{
  agbox *agb = mk_agbox(0.0,0.0,512.0,512.0);
  int i;
  for ( i = 0 ; i < lingraphs_size(lgs) ; i++ )
  {
    lingraph *lg = lingraphs_ref(lgs,i);
    agbox *sub = 
      mk_agbox_horizontal_substripe(agb,i,lingraphs_size(lgs));
    render_lingraph_in_agbox(sub,lg);
    free_agbox(sub);
  }
  free_agbox(agb);
}



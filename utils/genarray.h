/*
   File:        genarray.h
   Author:      Andrew W. Moore
   Created:     Sat Jun 26 03:22:44 EDT 1999
   Description: Header for Generic Array of Objects

   Copyright 1999, Schenley Park Research
*/

#ifndef GENARRAY_H
#define GENARRAY_H
#include "amut.h"

typedef struct generic_array
{
  int size;
  int array_size;
  void **array;

  void (*free_data)(void *data);
  void * (*mk_copy_data)(void *data);
  void (*fprintf_data)(FILE *s,char *m1,void *data,char *m2); /* May be NULL */
} generic_array;

#define INITIAL_GENERIC_ARRAY_SIZE 10

generic_array *mk_empty_generic_array(
  void (*free_data)(void *data),
  void * (*mk_copy_data)(void *data),
  void (*fprintf_data)(FILE *s,char *m1,void *data,char *m2));

/*
   Assume generic_array is previously of size n. After this it is of size
   n+1, and the n+1'th element is a COPY of thing.
*/
void add_to_generic_array(generic_array *ga,void *thing);

/*
   Assume generic_array is previously of size n. After this it is of size
   n+1, and the n+1'th element is the same pointer to the same thing. After
   calling this, thing should never be directly referenced or 
   directly freed again.
*/
void add_pointer_to_generic_array(generic_array *ga,void *thing);

int generic_array_size(generic_array *ga);

/*
   Returns a pointer (not a copy); to the index'th element stored in
   the generic_array. Error if index < 0 or index >= size
*/
void *generic_array_ref(generic_array *ga,int index);

/* Sets the index'th element to be a COPY of 'element', and frees what
   was there before*/
void generic_array_set(generic_array *ga,int index,void *element);
  
void fprintf_generic_array(FILE *s,char *m1,generic_array *ga,char *m2);

void free_generic_array(generic_array *ga);

generic_array *mk_copy_generic_array(generic_array *ga);

/* Increases the length of ga by 1, puts a COPY of element in as the
   i'th element, and any elements to the right of the i'th element (i.e.
   any elements with index >= i) are moved to having one higher index.
*/
void insert_in_generic_array(generic_array *ga,int i,void *element);

/* Decreases the length of ga by 1, removes the i'th element
   and any elements to the right of the i'th element (i.e.
   any elements with index >= i) are moved to having one lower index.
*/
void generic_array_remove(generic_array *ga,int i);

/* The entry contents of ga[i] are AMFREED. The pointer to sp
   goes in its place with NO COPYING. After this call sp should never again
   be referenced directly or freed directly. */
void generic_array_set_pointer(generic_array *ga,int i,void *item);

/* The pointer to the contents of ga[i] are overwritten with NULL.
   The previous contents are NOT freed, so someone else must have taken
   care of (or be about to take care of) that */
void generic_array_null_without_freeing(generic_array *genarr,int i);

#endif /* #ifndef GENARRAY_H */

#ifndef BOOLVEC_H
#define BOOLVEC_H

#include "amma.h"

typedef struct boolvec
{
  int array_size; /* Size of array in bytes */
  int size;       /* Number of bits in array */
  char *array;
} boolvec, *boolvec_ptr;

bool safe_boolvec_ref(boolvec *bv,int i);

extern int boolvec_mask[8];

#ifdef AMFAST
#define boolvec_ref(bv,i) (bv->array[(i)>>3]&boolvec_mask[(i)&0x07])
#else
#define boolvec_ref(bv,i) safe_boolvec_ref(bv,i)
#endif

boolvec *mk_boolvec(int size);
boolvec *mk_boolvec_init(int size, bool init);
void free_boolvec(boolvec *bv);
void boolvec_set(boolvec *bv,int i,bool value);
int boolvec_size(boolvec *bv);
int find_index_in_boolvec(boolvec *bv, bool value);
bool is_in_boolvec(boolvec *bv, bool value);
void fprintf_boolvec(FILE *s,char *m1,boolvec *bv,char *m2);
void copy_boolvec(boolvec *src,boolvec *dst);
boolvec *mk_copy_boolvec(boolvec *bv);
void add_to_boolvec(boolvec *bv,bool value);
void boolvec_remove(boolvec *bv,int index);
void boolvec_remove_last_element(boolvec *bv);

#endif

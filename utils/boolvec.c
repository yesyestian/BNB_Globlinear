#include "boolvec.h"

/* Boolvecs are adjustable-length vectors of bools.
   They are implemented efficiently so that the memory for a vector
   of length N is N (+ small constant) bits.

   These may be useful enough to eventually deserve porting to damut */

/* BOOLVECS ASSUME THAT CHARS ARE 8 BITS */

int boolvec_mask[8] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};

boolvec *mk_boolvec_init(int size, bool init)
{
  boolvec *bv;
  int i;
  char val = init ? 0xFF : 0x00;
  bv = mk_boolvec(size);
  for (i = 0; i < bv->array_size; i++)
    bv->array[i] = val;
  return (bv);
}

boolvec *mk_boolvec(int size)
{
  boolvec *bv = AM_MALLOC(boolvec);
  bv -> array_size = ((size & 0x07)==0) ? size>>3 : (size>>3)+1;
  bv -> size = size;
  if ( size < 0 ) my_error("mk_boolvec: bad size");

  bv -> array = AM_MALLOC_ARRAY(char,bv->array_size);
  return(bv);
}

void free_boolvec(boolvec *bv)
{
  AM_FREE_ARRAY(bv->array,char,bv->array_size);
  AM_FREE(bv,boolvec);
}

void check_boolvec_access(boolvec *bv,int i,char *mess)
{
  if ( bv==NULL || i >= bv->size || i < 0 )
  {
    printf("bad boolvec access in %s. index = %d, bv->size = %d\n",
           mess,i,bv->size);
    my_error("check_boolvec_access");
  }
}

bool safe_boolvec_ref(boolvec *bv,int i)
{
  int char_index = i >> 3;
  int bit = i & 0x7;

  check_boolvec_access(bv,i,"boolvec_ref");

  return((boolvec_mask[bit] & bv->array[char_index]) != 0);
}

void boolvec_set(boolvec *bv,int i,bool value)
{
  int char_index = i >> 3;
  int bit = i & 0x7;
  int boolvec_mask = 1 << bit;

#ifndef AMFAST
  check_boolvec_access(bv,i,"boolvec_ref");
#endif

  if ( value )
    bv->array[char_index] |= boolvec_mask;  /* a |= b means a := a bitwise-or b */
  else
  {
    int antimask = boolvec_mask ^ 0xFF;  /* a ^ b is a bitwise-exclusive-or b */
    bv->array[char_index] &= antimask;
  }
}

int boolvec_size(boolvec *bv)
{
  if ( bv==NULL )
    my_error("boolvec_size: NULL bv");
  return(bv->size);
}


int find_index_in_boolvec(boolvec *bv, bool value)
{
  int result = -1;
  int i;

  for ( i = 0 ; i < boolvec_size(bv) && result < 0 ; i++ )
    if (value == boolvec_ref(bv,i)) 
      result = i;
  return(result);
}


bool is_in_boolvec(boolvec *bv, bool value)
{
  return(find_index_in_boolvec(bv,value) >= 0);
}

void fprintf_boolvec(FILE *s,char *m1,boolvec *bv,char *m2)
{
  int i;
  fprintf(s,"%s = ",m1);
  if ( bv==NULL )
    printf("(boolvec *)NULL");
  else if ( boolvec_size(bv) == 0 )
    printf("<zero-length boolvec>");
  else
    for ( i = 0 ; i < boolvec_size(bv) ; i++ )
      fprintf(s,"%c",(boolvec_ref(bv,i)) ? '1' : '0');
  fprintf(s,"%s",m2);
}

void pboolvec(boolvec *x)
{
  fprintf_boolvec(stdout,"boolvec",x,"\n");
}

void copy_boolvec(boolvec *src,boolvec *dst)
{
  int i;
  if ( boolvec_size(src) != boolvec_size(dst) )
    my_error("copy_boolvec: wrong size");
  for ( i = 0 ; i < src->array_size ; i++ )
    dst->array[i] = src->array[i];
}

boolvec *mk_copy_boolvec(boolvec *bv)
{
  boolvec *newbv = mk_boolvec(boolvec_size(bv));
  copy_boolvec(bv,newbv);
  return(newbv);
}

void add_to_boolvec(boolvec *bv,bool value)
{
  int bits_in_array = bv->array_size << 3;
  if ( bv->size > bits_in_array ) 
    my_error(",zmnxocwndpcnloikasx");
  else if ( bv->size == bits_in_array )
  {
    int new_array_size = 2 + (int) (1.5 * bv->array_size);
    char *new_array = AM_MALLOC_ARRAY(char,new_array_size);
    int i;
    for ( i = 0 ; i < bv->array_size ; i++ )
      new_array[i] = bv->array[i];
    AM_FREE_ARRAY(bv->array,char,bv->array_size);
    bv->array = new_array;
    bv->array_size = new_array_size;
  }
  bv->size += 1;
  boolvec_set(bv,bv->size-1,value);
}


/**** Removal functions on boolvecs ****/

/* Reduces the size of bv by one.
   boolvec_ref(bv,index) dibvppears.
   Everything to the right of boolvec_ref(bv,index) is copied one to the left.

   Formally: Let bvold be the boolvec value beore calling this function
             Let bvnew be the boolvec value after calling this function.

PRE: boolvec_size(bvold) > 0
     0 <= index < boolvec_size(bvold)

POST: boolvec_size(bvnew) = boolvec_size(bvold)-1
      for j = 0 , 1, 2 ... index-1  : 
         boolvec_ref(bvnew,j) == boolvec_ref(bvold,j)

      for j = index , index+1 , ... boolvec_size(bvnew)-1:
         boolvec_ref(bvnew,j) == boolvec_ref(bvold,j+1)
*/
void boolvec_remove(boolvec *bv,int index)
{
  int i;
  if ( boolvec_size(bv) <= 0 ) my_error("boolvec_remove: empty boolvec");
  if ( index < 0 || index >= boolvec_size(bv) )
    my_error("boolvec_remove: bad index");

  for ( i = index ; i < boolvec_size(bv) - 1; i++ )
    boolvec_set(bv,i,boolvec_ref(bv,i+1));
  bv -> size -= 1;
}

/* Shrinks bv by one element by removing the rightmost element. 
   Example:

  Before: bv == 00100101
    boolvec_remove_last_element(bv)
  After:  bv == 0010010
*/
void boolvec_remove_last_element(boolvec *bv)
{
  boolvec_remove(bv,boolvec_size(bv)-1);
}

void boolvec_main(int argc,char *argv[])
{
  boolvec *bv = mk_boolvec(0);
  int c = '0';

  while ( c != 'q' )
  {
    printf("Hit 0 1 d or q> ");
    c = getchar();
    if ( c == '0' )
      add_to_boolvec(bv,FALSE);
    if ( c == '1' )
      add_to_boolvec(bv,TRUE);
    if ( c == 'd' )
      boolvec_remove_last_element(bv);
    if ( c == 'r' )
      boolvec_remove(bv,input_int("Remove which element? "));
    pboolvec(bv);
  }
  free_boolvec(bv);
}

/* END OF BOOLVECS */

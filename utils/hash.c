/*
   File:        hash.c
   Author:      Jeff Schneider
   Created:     Fri Feb 12 15:11:16 EST 1999
   Description: Hash table implementation

   Copyright (C) 1999, Jeff Schneider
*/

#include "hash.h"

/* This function defaults to a character key type.  See set_hash_table_key_type
   to change this setting. */
hash_table *mk_sized_empty_hash_table(int size,
				      void *(*data_copy_fn)(void *data),
				      void (*data_free_fn)(void *data),
				      int key_type)
{
  int i;
  hash_table *ht = AM_MALLOC(hash_table);
  ht->num_bins = size;
  ht->num_elements = 0;
  ht->key_type = key_type;
  ht->data_copy_fn = data_copy_fn;
  ht->data_free_fn = data_free_fn;
  ht->array = (hlist **)am_malloc(ht->num_bins * sizeof(hlist *));
  for (i=0;i<ht->num_bins;i++) ht->array[i] = NULL;
  return ht;
}

/* key_type should be one of the types #defined in hash.h.  
   WARNING:  it is NOT ok to call this function after you've already begun
   storing data in the hash table.  The mapping to the data you've already
   stored will be lost.
*/
void set_hash_table_key_type(hash_table *ht, int key_type)
{
  if (ht->num_elements > 0)
    my_error("set_hash_table_key_type: don't change key_type after some data is stored in the hash table!\n");
  ht->key_type = key_type;
}

hash_table *mk_empty_hash_table(void *(*data_copy_fn)(void *data),
				void (*data_free_fn)(void *data))
{
  return mk_sized_empty_hash_table(INITIAL_HASH_BINS,
				   data_copy_fn,data_free_fn,CHAR_KEY_TYPE);
}

/* This function copied from htab.c in planit.
   It had no comments there either.
 */
int hash_code_from_string(char *s, int num_bins)
{
  unsigned int h = 0;
  unsigned int g;
  int i;
  int odd_number = (num_bins % 2 == 0) ? num_bins - 1 : num_bins;

  for (i=0;s[i];i++)
  {
    h = (h << 4) + s[i];
    g = h & 0xf0000000;
    if (g)
    {
      h = h ^ (g >> 24);
      h = h ^ g;
    }
  }
  return (h % odd_number);
}

int hash_code_from_ivec(ivec *iv, int num_bins)
{
  unsigned int h = 0;
  unsigned int g;
  int i;
  int odd_number = (num_bins % 2 == 0) ? num_bins - 1 : num_bins;

  for (i=0;i<ivec_size(iv);i++)
  {
    h = (h << 4) + ivec_ref(iv,i);
    g = h & 0xf0000000;
    if (g)
    {
      h = h ^ (g >> 24);
      h = h ^ g;
    }
  }
  return (h % odd_number);
}

hlist *lookup_list_entry_node(hash_table *ht, hlist *hl, char *key)
{
  while(hl)
  {
    if (ht->key_type == CHAR_KEY_TYPE)
    {
      if (eq_string(key,hl->key)) return hl;
    }
    else if (ht->key_type == IVEC_KEY_TYPE)
    {
      if (ivec_equal((ivec *)key,(ivec *)hl->key)) return hl;
    }
    hl = hl->next;
  }
  return NULL;
}

hlist *lookup_hash_entry_node(hash_table *ht, char *key)
{
  int index = 0;

  if (ht->key_type == CHAR_KEY_TYPE)
    index = hash_code_from_string(key,ht->num_bins);
  else if (ht->key_type == IVEC_KEY_TYPE)
    index = hash_code_from_ivec((ivec *)key,ht->num_bins);
  else
    my_error("lookup_hash_entry_node: invalid hash table key_type\n");
  return lookup_list_entry_node(ht,ht->array[index],key);
}

void *lookup_hash_entry(hash_table *ht, char *key)
{
  hlist *node = lookup_hash_entry_node(ht,key);
  if (node) return node->data;
  else      return NULL;
}

/* This function increases the number of bins in a hash table by taking the
   following steps:
   1. Create a new, empty, larger hash table
   2. Go through the elements of the old table and re-enter them into the new
   3. Swap the guts of the two structures so the new one lives in the old shell
   4. Free the new structure (which holds the old guts)
 */
void increase_hash_size(hash_table *ht, int multiple)
{
  hash_table *new_ht = mk_sized_empty_hash_table(ht->num_bins*multiple,
						 ht->data_copy_fn,
						 ht->data_free_fn,
						 ht->key_type);
  int i,temp;
  hlist **tempa;

  for (i=0;i<ht->num_bins;i++)
  {
    hlist *hl = ht->array[i];
    while(hl)
    {
      add_hash_entry(new_ht,hl->key,hl->data);
      hl=hl->next;
    }
  }
  temp=ht->num_bins;ht->num_bins=new_ht->num_bins;new_ht->num_bins=temp;
  temp=ht->num_elements;ht->num_elements=new_ht->num_elements;new_ht->num_elements=temp;
  tempa=ht->array;ht->array=new_ht->array;new_ht->array=tempa;
  free_hash_table(new_ht);
}

/* The data is copied in */
void add_hash_entry(hash_table *ht, char *key, void *data)
{
  int index = -7777;
  hlist *hl;

  if (ht->num_elements > (3 * ht->num_bins)) increase_hash_size(ht,10);

  if (ht->key_type == CHAR_KEY_TYPE)
    index = hash_code_from_string(key,ht->num_bins);
  else if (ht->key_type == IVEC_KEY_TYPE)
    index = hash_code_from_ivec((ivec *)key,ht->num_bins);
  else
    my_error("add_hash_entry: invalid hash table key_type\n");

  hl = AM_MALLOC(hlist);
  if (ht->key_type == CHAR_KEY_TYPE) 
    hl->key = mk_copy_string(key);
  else if (ht->key_type == IVEC_KEY_TYPE)
    hl->key = (char *)mk_copy_ivec((ivec *)key);
  hl->data = ht->data_copy_fn(data);
  hl->next = ht->array[index];
  ht->array[index] = hl;
  ht->num_elements++;
}

void *update_hash_entry(hash_table *ht, char *key, void *data)
{
  hlist *hl = lookup_hash_entry_node(ht,key);
  if (hl)
  {
    ht->data_free_fn(hl->data);
    hl->data = ht->data_copy_fn(data);
    return data;
  }
  else return NULL;
}

hlist *mk_copy_hlist(hash_table *ht, hlist *hl,
		     void *(*data_copy_fn)(void *data))
{
  if (hl)
  {
    hlist *new_hl = AM_MALLOC(hlist);
    if (ht->key_type == CHAR_KEY_TYPE)
      new_hl->key = mk_copy_string(hl->key);
    else if (ht->key_type == IVEC_KEY_TYPE)
      new_hl->key = (char *)mk_copy_ivec((ivec *)hl->key);
    else
      my_error("mk_copy_hlist: invalid hash table key_type\n");
    new_hl->data = data_copy_fn(hl->data);
    new_hl->next = mk_copy_hlist(ht,hl->next,data_copy_fn);
    return new_hl;
  }
  else return NULL;
}

hash_table *mk_copy_hash_table(hash_table *ht)
{
  int i;
  hash_table *new_ht = AM_MALLOC(hash_table);
  new_ht->num_bins = ht->num_bins;
  new_ht->num_elements = ht->num_elements;
  new_ht->key_type = ht->key_type;
  new_ht->data_copy_fn = ht->data_copy_fn;
  new_ht->data_free_fn = ht->data_free_fn;
  new_ht->array = (hlist **)am_malloc(ht->num_bins * sizeof(hlist *));
  for (i=0;i<new_ht->num_bins;i++)
    new_ht->array[i] = mk_copy_hlist(ht,ht->array[i],new_ht->data_copy_fn);
  return new_ht;
}

void free_hlist_element(hash_table *ht, hlist *hl,
			void (*data_free_fn)(void *data))
{
  data_free_fn(hl->data);
  if (ht->key_type == CHAR_KEY_TYPE) free_string(hl->key);
  else if (ht->key_type == IVEC_KEY_TYPE) free_ivec((ivec *)hl->key);
  else my_error("free_hlist_element: invalid hash table key_type\n");
  AM_FREE(hl,hlist);
}

void free_hlist(hash_table *ht,hlist *hl,void (*data_free_fn)(void *data))
{
  if (hl->next) free_hlist(ht,hl->next,data_free_fn);
  free_hlist_element(ht,hl,data_free_fn);
}

void free_hash_table(hash_table *ht)
{
  int i;
  for (i=0;i<ht->num_bins;i++) 
    if (ht->array[i]) free_hlist(ht,ht->array[i],ht->data_free_fn);
  am_free((char *)ht->array,ht->num_bins*sizeof(hlist *));
  AM_FREE(ht,hash_table);
}

void delete_hash_entry(hash_table *ht, char *key)
{
  int index = 0;   
  hlist *prev = NULL;
  hlist *cur = NULL;

  if (ht->key_type == CHAR_KEY_TYPE)
    index = hash_code_from_string(key,ht->num_bins);
  else if (ht->key_type == IVEC_KEY_TYPE)
    index = hash_code_from_ivec((ivec *)key,ht->num_bins);
  else
    my_error("delete_hash_entry: invalid hash table key_type\n");

  cur = ht->array[index];

  while (cur && strcmp(key,cur->key))
  {
    prev = cur;
    cur = cur->next;
  }
  if (cur) /* that means we found it */
  {
    if (prev) prev->next = cur->next;
    else      ht->array[index] = cur->next;
    free_hlist_element(ht,cur,ht->data_free_fn);
    ht->num_elements--;
  }
}

/* A debugging function for those who are concerned about how well the
   hash function and hash table are working.
*/
#define PHTS_HIST_SIZE 6
void print_hash_table_stats(hash_table *ht)
{
  int i;
  ivec *counts = mk_zero_ivec(PHTS_HIST_SIZE);
  printf("%d bins and %d elements => %.2f full\n",
	 ht->num_bins,ht->num_elements,
	 (double)ht->num_elements/(double)ht->num_bins);
  printf("Note: it is reasonable for the fullness to be above or below one,\nbut probably not more than an order of magnitude either way.\n");

  for (i=0;i<ht->num_bins;i++)
  {
    int count = 0;
    hlist *hl = ht->array[i];
    while (hl) {count++; hl=hl->next;}
    if (count < PHTS_HIST_SIZE) ivec_increment(counts,count,1);
    else                        ivec_increment(counts,PHTS_HIST_SIZE-1,1);
  }
  printf("Elements/Bin Percentage_of_bins_with Number_of_bins_with\n");
  for (i=0;i<PHTS_HIST_SIZE;i++)
  {
    if (i<PHTS_HIST_SIZE-1) printf("%8d ",i);
    else                    printf("%7d+ ",i);
    printf("%10.2f%%  ",100.0*(double)ivec_ref(counts,i)/(double)ht->num_bins);
    printf("   %10d\n",ivec_ref(counts,i));
  }
  free_ivec(counts);
}

/* Some basic copy and free functions */

void *vint_copy(void *data)
{
  int *a = (int *)am_malloc(sizeof(int));
  *a = *((int *)data);
  return (void *)a;
}

void vint_free(void *data)
{
  am_free((char*) data,sizeof(int));
}

/* String-Hash is a simple set-of-strings in which membership
   test and addition of a string and deletion is a constant-time
   operation. Implemented using the hash structure defined above. */

stringhash *mk_empty_stringhash()
{
  stringhash *sh = AM_MALLOC(stringhash);
  sh -> ht = mk_empty_hash_table(vint_copy,vint_free);
  return sh;
}

void free_stringhash(stringhash *sh)
{
  free_hash_table(sh->ht);
  AM_FREE(sh,stringhash);
}

void add_to_stringhash(stringhash *sh,char *string)
{
  if ( !is_in_stringhash(sh,string) )
  {
    int fake_data[1];
    fake_data[0] = 1;
    add_hash_entry(sh->ht,string,(void *) fake_data);
  }
}

void remove_from_stringhash(stringhash *sh,char *string)
{
  delete_hash_entry(sh->ht,string);
}

bool is_in_stringhash(stringhash *sh,char *string)
{
  void *x = lookup_hash_entry(sh->ht,string);
  return x != NULL;
}

string_array *mk_string_array_from_stringhash(stringhash *sh)
{
  string_array *sa = mk_string_array(0);
  int i;
  for ( i = 0 ; i < sh->ht->num_bins ; i++ )
  {
    hlist *hl;
    for ( hl = sh->ht->array[i] ; hl != NULL ; hl = hl -> next )
      add_to_string_array(sa,hl->key);
  }
  return sa;
}

string_array *mk_sorted_string_array_from_stringhash(stringhash *sh)
{
  string_array *sa = mk_string_array_from_stringhash(sh);
  sort_string_array(sa);
  return sa;
}

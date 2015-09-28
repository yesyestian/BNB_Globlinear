/*
   File:        hash.h
   Author:      Jeff Schneider
   Created:     Fri Feb 12 15:11:16 EST 1999
   Description: Header for hash table implementation

   Copyright (C) 1999, Jeff Schneider
*/


#ifndef HASH_H
#define HASH_H
#include "amut.h"

/* This is a very simple implementation of hash tables.  You begin by
   calling mk_empty_hash_table().  You may also want to call
   set_hash_table_key_type if you want to use a key different than the
   default (character string).  Each entry in the hash table will
   have a data element (defined by you) of type (void *).  Therefore, you
   must provide two functions to the creation routine.  One is a copy
   function that takes a pointer to your data objects, and allocates a new
   copy of the data object.  The other is a free function, which takes a
   pointer to one of your data objects and frees it.  Having done that
   single call, you're now ready to use the hash table with the following
   functions:
   lookup_hash_entry - you provide a key (a character string) and it returns
                       the corresponding data pointer (or NULL) if that key
                       is not stored in the table.
   add_hash_entry    - you provide a key and a pointer to a data object and
                       the new data is added.  Note that this does NOT check
                       to see whether the key already has been placed in the
                       table.  You must do that yourself first using the
                       lookup function.
   delete_hash_entry - you provide the key, and it deletes the entry if it
                       finds it, otherwise does nothing.
   update_hash_entry - you give it a new data pointer to be stored with the
                       associated key.  if the key was found and changed, it
		       returns the given data pointer, otherwise NULL.
   free_hash_table   - call this when you're done unless you're writing a
                       Microsoft operating system.

   The hash table automatically grows itself based on the number of new
   entries you give it.  It does not automatically shrink again if you
   delete most of the entries.

   There is a test function for this code in main.c.  You can look at it to
   see a simple example of using this code to store integers with keys.
*/

#define INITIAL_HASH_BINS 10

#define CHAR_KEY_TYPE 0
#define IVEC_KEY_TYPE 1

typedef struct hlist
{
  char *key;
  void *data;
  struct hlist *next;
} hlist;

typedef struct hash_table
{
  int num_bins;
  int num_elements;
  int key_type;
  void *(*data_copy_fn)(void *data);
  void (*data_free_fn)(void *data);
  hlist **array;
} hash_table;

hash_table *mk_empty_hash_table(void *(*data_copy_fn)(void *data),
				void (*data_free_fn)(void *data));

void set_hash_table_key_type(hash_table *ht, int key_type);

void *lookup_hash_entry(hash_table *ht, char *key);

void add_hash_entry(hash_table *ht, char *key, void *data);

void delete_hash_entry(hash_table *ht, char *key);

void *update_hash_entry(hash_table *ht, char *key, void *data);

hash_table *mk_copy_hash_table(hash_table *ht);

void free_hash_table(hash_table *ht);

/* A debugging function for those who are concerned about how well the
   hash function and hash table are working.
*/
void print_hash_table_stats(hash_table *ht);

/* Some basic copy and free functions */

void *vint_copy(void *data);
void vint_free(void *data);

/* String-Hash is a simple set-of-strings in which membership
   test and addition of a string and deletion is a constant-time
   operation. Implemented using the hash structure defined above. */

typedef struct stringhash
{
  hash_table *ht;
} stringhash;

/* String-Hash is a simple set-of-strings in which membership
   test and addition of a string and deletion is a constant-time
   operation. Implemented using the hash structure defined above. */

stringhash *mk_empty_stringhash();
void free_stringhash(stringhash *sh);
void add_to_stringhash(stringhash *sh,char *string);
void remove_from_stringhash(stringhash *sh,char *string);
bool is_in_stringhash(stringhash *sh,char *string);
string_array *mk_string_array_from_stringhash(stringhash *sh);
string_array *mk_sorted_string_array_from_stringhash(stringhash *sh);

#endif /* #ifndef HASH_H */


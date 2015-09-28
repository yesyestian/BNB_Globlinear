/*
   File:         amma.c
   Author:       Andrew W. Moore
   Created:      4th September 1992
   Description:  My own memory allocation and deallocation.

   Copyright (C) Andrew W. Moore, 1992
*/

#define TRACE_MEMORY 0x00986998

#ifdef WIN32
#include <windows.h>
CRITICAL_SECTION mutex;
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "amma.h"
#include "amstr.h"

int Total_mallocked = 0;
int Total_freed = 0;
int Max_mallocked = 0;
int Current_large_mallocked = 0;
char Dummy_memory;
int Num_zeroes_allocated = 0;

char *am_basic_malloc(int size)
{
  char *result;
  if ( size < 0 )
  {
    result = NULL;
    fprintf(stderr,"am_basic_malloc(%d)?\n",size);
    my_error("no way dude");
  }
  else if ( size == 0 )
    result = NULL;
  else
  {
    result = (char *) malloc((unsigned int) size);
    if ( result == NULL )
    {
      basic_am_malloc_report();
      fprintf(stderr,"Oh no!!! We're out of memory! In total you mallocked %d bytes\n",
	      Total_mallocked
	      );
      fprintf(stderr,"However, I returned %d bytes via free()\n",Total_freed);
      fprintf(stderr,"So we have currently got %d bytes in use.\n",
              Total_mallocked - Total_freed
	      );
      fprintf(stderr,"We got back NULL from a call to O/S malloc when you asked\n");
      fprintf(stderr,"for something of size %d\n",size);     
      wait_for_key();
      my_error("out of memory");
    }
    Total_mallocked += size;
    Max_mallocked = int_max(Max_mallocked,Total_mallocked - Total_freed);
  }

  return(result);
}

#ifndef AMFAST
typedef struct amnode
{
  int count;
  int size;
  struct amnode *next;
} amnode;

amnode *Free_amnodes = NULL;

amnode *alloc_amnode()
{
  amnode *am;
  if ( Free_amnodes == NULL )
  {
    am = (amnode *) am_basic_malloc(sizeof(amnode));
    am -> count = -7777;
    am -> size = -7777;
    if ( am == NULL ) my_error("out of memory");
  }
  else
  {
    am = Free_amnodes;
    Free_amnodes = am->next;
  }
  return am;
}

void free_amnode(amnode *amn)
{
  amn->next = Free_amnodes;
  amn->count = -777;
  amn->size = -777;
  Free_amnodes = amn;
}

amnode *add_to_amnode(amnode *am,int count,int size)
{
  amnode *result = alloc_amnode();
  result -> count = count;
  result -> size = size;
  result -> next = am;
  return result;
}

amnode *remove_from_amnode(amnode *am,int count,int size,bool *r_found)
{
  amnode *result = NULL;
  amnode *removed = NULL;
  if ( am == NULL )
    *r_found = FALSE;
  else if ( am->count == count )
  {
    *r_found = TRUE;
    result = am->next;
    removed = am;
  }
  else
  {
    result = am;
    while ( am->next != NULL && am->next->count != count )
      am = am->next;
    if ( am->next == NULL )
    {
      result = NULL;
      *r_found = FALSE;
    }
    else
    {
      amnode *new_next = am->next->next;
      *r_found = TRUE;
      removed = am->next;
      am->next = new_next;
    }
  }

  if ( removed != NULL )
  {
    if ( size != removed->size )
    {
      my_error(mk_printf("You're freeing memory you say has size %d but I say it was am_mallocked with size %d",
                         size,removed->size));
    }
    free_amnode(removed);
  }

  return result;
}

#define NUM_COUNT_BUCKETS 7777
amnode *Amnode_buckets[NUM_COUNT_BUCKETS];
bool Amnode_buckets_started = FALSE;

void add_count_to_hash(int count,int size)
{
  int hash = count % NUM_COUNT_BUCKETS;

  if ( !Amnode_buckets_started )
  {
    int i;
    for ( i = 0 ; i < NUM_COUNT_BUCKETS ; i++ )
      Amnode_buckets[i] = NULL;
    Amnode_buckets_started = TRUE;
  }
  Amnode_buckets[hash] = add_to_amnode(Amnode_buckets[hash],count,size);
}

bool can_remove_count_from_hash(int count,int size)
{
  int hash = count % NUM_COUNT_BUCKETS;
  bool found;
  amnode *newlist;

  if ( !Amnode_buckets_started ) my_error("no way");
  
  if ( Amnode_buckets[hash] == NULL )
    my_error("Serious error freeing something...it looks legal but in fact nothing of that size is allocated");

  newlist = remove_from_amnode(Amnode_buckets[hash],count,size,&found);
  if ( found )
    Amnode_buckets[hash] = newlist;

  return found;
}
    
int lowest_count_in_hash()
{
  int result = -1;
  int i;
  for ( i = 0 ; i < NUM_COUNT_BUCKETS ; i++ )
  {
    amnode *am;
    for ( am = Amnode_buckets[i] ; am != NULL ; am = am->next )
      if ( result < 0 || result > am->count )
        result = am->count;
  }
  return result;
}
#endif

typedef struct free_list_struct
{
  char *memory;
  struct free_list_struct *next;
} free_list;

free_list *Frees[MAX_AM_MALLOC_SIZE];
int Num_allocated[MAX_AM_MALLOC_SIZE];

bool Started = FALSE;

free_list *Free_frees = NULL;
int Frees_allocated = 0;

free_list *alloc_free_node(void)
{
  free_list *result;
  if ( Free_frees != NULL )
  {
    result = Free_frees;
    Free_frees = Free_frees -> next;
  }
  else
    result = (free_list *) am_basic_malloc(sizeof(free_list));

  Frees_allocated += 1;

  return(result);
}

void free_free_node(free_list *fr)
{
  Frees_allocated -= 1;
  fr -> next = Free_frees;
  Free_frees = fr;
}

int Num_previous_am_mallocs = 0;
int Stop_at_nth_malloc = -1;
int Stop_at_defined_by_argc_argv = FALSE;

/* Call this function with am_malloc_number set to N if you are trying
   to debug a memory leak and if am_malloc_report() has warned you that
   on am_malloc call number N there was a leak. */
void memory_leak_stop(int am_malloc_number)
{
#ifdef AMFAST
  printf("In order to do memory-leak debugging you must have\n");
  printf("ALL the sources and auton libraries compiled in non-AMFAST mode.\n");
  my_error("You cannot use memory_leak_stop() if running in AMFAST mode.\n");
#endif
  if ( am_malloc_number < Num_previous_am_mallocs )
  {
    printf("It is useless to call memory_leak_stop(%d) at this point\n",
           am_malloc_number);
    printf("in the code. You have already made %d am_malloc calls.\n",
           Num_previous_am_mallocs);
    my_error("memory_leak_stop()");
  }
  Stop_at_nth_malloc = am_malloc_number;
  Stop_at_defined_by_argc_argv = FALSE;
  printf("Will stop the program and let you enter debugger at am_malloc number %d\n",
         am_malloc_number);
}
/* Call this function at the start of main(int argc,char *argv[]) thus:

     memory_leak_check_args(argc,argv)

   Then, if ever am_malloc_report tells you that you have a memory leak on
   the <N>th call the am_malloc, simply include 
      memleak <N>
   on the command line.
*/
void memory_leak_check_args(int argc,char *argv[])
{
  extern int Ignore_next_n;

  Ignore_next_n = int_from_args("num_waitforkey_skips",argc,argv,0);

  if ( index_of_arg("memleak",argc,argv) > 0 )
  {
    int am_malloc_number = int_from_args("memleak",argc,argv,-1);
    if ( am_malloc_number < 0 )
      printf("Ignoring memleak %d on command line (because -ve value)\n",
             am_malloc_number);
    else
    {
#ifdef AMFAST
      printf("In order to do memory-leak debugging you must have\n");
      printf("ALL the sources and auton libraries compiled in non-AMFAST mode.\n");
      printf("You cannot have memleak %d on the command line\n"
              "if running in AMFAST mode.\n",am_malloc_number);
      my_error("memory_leak_check_args()");
#endif
      memory_leak_stop(am_malloc_number);
      Stop_at_defined_by_argc_argv = TRUE;
    }
  }
}

void memory_leak_stop_message()
{
  printf("Entered memory_leak_stop_message(). Set a breakpoint for\n");
  printf("this function if you want you stop when the memory-leaking\n");
  printf("am_malloc() is called.\n");
}
 
#define BYTES_PER_TAG 8

/* If not AMFAST...

      Returns a block of memory of size "base_size + BYTES_PER_TAG bytes".
      In the first four bytes is a count of the number of
      previous calls to am_extended_malloc that had occured
      before this one. This is used in hunting memory leaks.

      At the moment Bytes_per_tag is eight. What do we do with the remaining
      four bytes? Nothing. They are there to help prevent doubles array
      alignment problems.

   If AMFAST

      Returns a block of memory of size "size".
*/
char *am_extended_malloc(int base_size)
{
#ifdef AMFAST
  int ext_size = base_size;
#else
  int ext_size = base_size + BYTES_PER_TAG;
#endif

  char *memory;

  int nhp2_size = next_highest_power_of_two(ext_size);
  /* We will actually malloc a block of memory of this size,
     UNLESS nhp2_size > MAX_AM_MALLOC_SIZE, in which case
     we'll only malloc a block of size "size" */


#ifndef AMFAST
  if ( Num_previous_am_mallocs == Stop_at_nth_malloc )
  {
    int mnum = Num_previous_am_mallocs;
    printf("I am stopping the program because I am about to\n");
    printf("perform am_malloc number %d, and because you\n",mnum);
    printf("had requested me to stop here by including\n\n");
    if ( Stop_at_defined_by_argc_argv )
      printf("    memleak %d\n\non the command line.",mnum);
    else
      printf("    memory_leak_stop(%d)\n\nin your C-code.\n",mnum);
    printf("Assuming you are hunting for memory leaks, you should\n");
    printf("now enter the debugger (^C in gdb), (select the BREAK\n");
    printf("menu item in Visual C++) and inspect the C-program calling\n");
    printf("stack. Your memory-leaking call will be on the stack.\n");
    printf("OR set a break point at memory_leak_stop_message() and hit return.\n");
    really_wait_for_key();
    memory_leak_stop_message();
    printf("Will continue the program when you next hit return.\n");
    really_wait_for_key();
  }
#endif

  if ( !Started )
  {    
    int i;
#ifdef WIN32
    InitializeCriticalSection(&mutex);
#endif
    for ( i = 0 ; i < MAX_AM_MALLOC_SIZE ; i++ )
    {
      Frees[i] = NULL;
      Num_allocated[i] = 0;
    }
    Started = TRUE;
  }
#ifdef WIN32
  EnterCriticalSection(&mutex);
#endif

  if ( base_size < 0 )
  {
    fprintf(stderr,"am_malloc(%d) illegal\n",base_size);
    memory = NULL;
    my_error("amma.c, am_malloc()");
  }
  else if ( ext_size == 0 )
  {
    memory = &Dummy_memory;
    Num_zeroes_allocated += 1;
  }
  else if ( nhp2_size >= MAX_AM_MALLOC_SIZE )
  {
    /* Note, malloc size is "size" not nhp2_size for these large memory sizes*/
    memory = am_basic_malloc(ext_size);
    Current_large_mallocked += ext_size;
  }
  else if ( Frees[nhp2_size] == NULL )
  {
    Num_allocated[ext_size] += 1;
    memory = am_basic_malloc(nhp2_size);
  }
  else
  {
    free_list *old = Frees[nhp2_size];
    memory = old -> memory;
    Frees[nhp2_size] = old -> next;
    old -> next = NULL;
    free_free_node(old);
    Num_allocated[ext_size] += 1;
  }

#ifdef WIN32
  LeaveCriticalSection(&mutex);
#endif

#ifndef AMFAST
  memory[0] = ((Num_previous_am_mallocs & 0xFF000000) >> 24);
  memory[1] = ((Num_previous_am_mallocs & 0x00FF0000) >> 16);
  memory[2] = ((Num_previous_am_mallocs & 0x0000FF00) >> 8);
  memory[3] = ((Num_previous_am_mallocs & 0x000000FF));

  add_count_to_hash(Num_previous_am_mallocs,base_size);

  Num_previous_am_mallocs += 1;
#endif

  return(memory);
}

/* Returns a pointer to a block of memory of
   size "size".

   If AMFAST... that's the end of the story (macro in amma.h)

   If not AMFAST,
       The BYTES_PER_TAG bytes before this block of memory
       are also allocated, a contain a count of the number
       of counts to am_malloc that preceded this one.
*/
#ifndef AMFAST
#ifndef USE_OS_MEMORY_MANAGEMENT
char *am_malloc(int size)
{
  char *mem_minus_BYTES_PER_TAG = am_extended_malloc(size);
  char *result = mem_minus_BYTES_PER_TAG + BYTES_PER_TAG;

  if ( TRACE_MEMORY == (long)result )
    printf("am_malloc(%d) returns %lX (count = %d)\n",size,(long)result,Num_previous_am_mallocs-1);

  return result;
}
#endif
#endif


bool in_free_list(free_list *fl,char *memory)
{
  bool result = FALSE;
  for ( ; fl != NULL && !result ; fl = fl -> next )
    result = ((long)(fl->memory)) == ((long)memory);
  return(result);
}

bool free_list_has_duplicates(free_list *fl)
{
  bool result = FALSE;
  fprintf(stderr,"   [does free list have duplicates?] ......\n");
  for ( ; fl != NULL && !result ; fl = fl -> next )
  {
    result = in_free_list(fl->next,fl->memory);
    if (result) fprintf(stderr,"Memory %p is duplicated\n",fl->memory);
  }
  fprintf(stderr,"   ......finished\n");
  return(result);
}

/*
   If AMFAST, frees memory block of given size
   If not AMFAST, frees the memory block of size "base_size+BYTES_PER_TAG"
*/
void am_extended_free(char *memory, int base_size)
{
#ifdef AMFAST
  int ext_size = base_size;
#else
  int ext_size = base_size + BYTES_PER_TAG;
#endif
  int nhp2_size;
#ifdef WIN32
  EnterCriticalSection(&mutex);
#endif

#ifndef AMFAST
  {
    int count = (memory[0] & 0xFF) << 24 |
                (memory[1] & 0xFF) << 16 |
                (memory[2] & 0xFF) <<  8 |
                (memory[3] & 0xFF);
    bool removed = can_remove_count_from_hash(count,base_size);
    if ( !removed )
    {
      printf("SERIOUS PROBLEM. YOU ARE FREEING SOMETHING THAT YOU "
             "PROBABLY NEVER ALLOCATED (OR MAYBE ALREADY FREED)\n");
      printf("(cryptic internal count code is %d)\n",count); 
      my_error("am_free");
    }
  }
#endif

  nhp2_size = next_highest_power_of_two(ext_size);

  if ( memory == NULL )
  {
    fprintf(stderr,"am_free(NULL,%d)\n",base_size);
    my_error("You can't free NULL");
  }
  if ( !Started )
    my_error("amma.c, am_free(), Not started");

  if ( base_size < 0 )
  {
    fprintf(stderr,"am_free(memory,%d) illegal\n",base_size);
    my_error("amma.c, am_free()");
  }
  else if ( nhp2_size >= MAX_AM_MALLOC_SIZE )
  {
    Total_freed += ext_size;
    Current_large_mallocked -= ext_size;
    free(memory);
  }
  else if ( ext_size == 0 )
  {
    if ( Num_zeroes_allocated <= 0 )
      my_error("free size 0? There's nothing of that size mallocked\n");
    Num_zeroes_allocated -= 1;
  }
  else if ( Num_allocated[ext_size] <= 0 )
  {
    basic_am_malloc_report();
    fprintf(stderr,
            "Free something size %d? Nothing am_malloced of that size\n",base_size
	    );
    fprintf(stderr,"Memory location of this thing is %8lX\n",(long) memory);

    if ( Frees[nhp2_size] == NULL )
      fprintf(stderr,
              "In fact, I never called am_malloc(%d) during entire program\n",
              base_size
	      );
    else
    {
      fprintf(stderr,
              "(Note: you may have previously mallocked & freed things that size)\n"
	      );
      if ( in_free_list(Frees[nhp2_size],memory) )
        fprintf(stderr," ---You've freed that exact piece of memory before\n");
      else
      {
        fprintf(stderr," ---But you've never freed that piece of memory\n");
        fprintf(stderr,"    (calling it size %d while freeing) before\n",base_size);
      }

      if ( free_list_has_duplicates(Frees[nhp2_size]) )
      {
        fprintf(stderr,"There are duplicate entries in the free list for\n");
        fprintf(stderr,"things of that size. You must have earlier on\n");
        fprintf(stderr,"freed the same piece of memory twice. That\n");
        fprintf(stderr,"confused my counters which is why you're getting\n");
        fprintf(stderr,
                "errors now instead of when you did the duplicate free.\n");
      }
      else
      {
        fprintf(stderr,"The free list for things this size has no\n");
        fprintf(stderr,"duplicates. That means this very am_free call\n");
        fprintf(stderr,"must be freeing a never-allocated thing. Or\n");
        fprintf(stderr,"something am_mallocked with a different size.\n");
        fprintf(stderr,"Or, earlier, you freed something else, wrongly.\n");
      }
    }
  
    my_error("amma.c, am_free()");
  }
  else
  {
    free_list *newv = alloc_free_node();
    newv -> memory = memory;
    newv -> next = Frees[nhp2_size];

#ifdef VERY_CAREFUL_CHECKING
    /* added by awm Sep 1 95 */
    if ( in_free_list(Frees[nhp2_size],memory) )
    {
      fprintf(stderr,"am_free(memory=%lx,size=%d)\n",(long)memory,base_size);
      my_error("am_free: freeing something already freed");
    }
#endif

    Frees[nhp2_size] = newv;
    Num_allocated[ext_size] -= 1;
  }

#ifdef WIN32
  LeaveCriticalSection(&mutex);
#endif
}

/* If AMFAST frees the memory (macro in amma.h)
   If not AMFAST, frees the memory, and the four bytes before.
*/
#ifndef AMFAST
#ifndef USE_OS_MEMORY_MANAGEMENT
void am_free(char *memory, int size)
{
  char *ptr = NULL;
  int size_two = 2 * size;

  if ( TRACE_MEMORY == (long)memory )
    printf("am_free(%lX,%d)\n",(long)memory,size);

  if ( 0 ) 
    printf("am_free(...,%d) 2*size = %d\n",size,size_two);
    

  ptr = memory - BYTES_PER_TAG;
  am_extended_free(ptr,size);
}
#endif
#endif

int free_length(free_list *fr)
{
  int result = 0;
  for ( ; fr != NULL ; fr = fr -> next )
    result += 1;
  return(result);
}

/* Give as much memory as we can find back to the operating sys via free().
 * This routine is intended for use from am_exit when our code has been
 * dynamically linked and we won't be doing a real exit() which would give
 * all the memory back
 */
void am_free_to_os()
{
  int i;
  free_list *cur,*tmp;

  cur = Free_frees;
  while(cur)
  {
    tmp = cur;
    cur = cur->next;
    free(tmp->memory);
    free(tmp);
    Total_freed += sizeof(free_list)*2;
  }
  for (i=0;i<MAX_AM_MALLOC_SIZE;i++)
    if (Frees[i])
    {
      cur = Frees[i];
      while (cur)
      {
	tmp = cur;
	cur = cur->next;
	free(tmp->memory);
	free(tmp);
	Total_freed += sizeof(free_list) + i;
      }
    }
}

void reset_malloc_state()
{
  am_free_to_os(); /* first give everything back we can */
  Started = FALSE;
  Free_frees = NULL;
  Frees_allocated = 0;
  Total_mallocked = 0;
  Total_freed = 0;
  Max_mallocked = 0;
  Current_large_mallocked = 0;
  Num_zeroes_allocated = 0;
#ifdef WIN32
  DeleteCriticalSection(&mutex);
#endif
}

void basic_am_malloc_report()
{
  int size;
  bool boring = TRUE;
  bool free_sizes_exist = FALSE;
  bool all_freed = TRUE;
  int free_bytes = 0;

  for ( size = 0 ; size < MAX_AM_MALLOC_SIZE ; size++ )
    free_bytes += ( sizeof(free_list) + size ) * free_length(Frees[size]);

  free_bytes += 2 * sizeof(free_list) * free_length(Free_frees);

  if ( Current_large_mallocked > 0 ) all_freed = FALSE;

  for ( size = 0 ; all_freed && size < MAX_AM_MALLOC_SIZE ; size++ )
    all_freed = Num_allocated[size] == 0;

  fprintf(stdout,"\n#\n# am_malloc report:                       %s is am_free()'d\n#\n",
         (all_freed) ? "EVERYTHING" : "NOT everything"
	 );

  if ( !all_freed )
  {
#ifdef AMFAST
    printf("#   Run the Debug version (AMFAST undefined during compiling)\n");
    printf("#   to find out how to track the leak.\n#\n");
#else
    printf("#   On am_malloc call number %d you got a piece of memory\n",
           lowest_count_in_hash());
    printf("#   that you never freed. To find this memory leak, EITHER\n");
    printf("#   add a call to memory_leak_stop(%d) at the start of your\n",
           lowest_count_in_hash());
    printf("#    main() program. OR (for general purpose use) make sure\n");
    printf("#    there's a call to memory_leak_check_args(argc,argv) at\n");
    printf("#    the start of the main() program, and then include\n");
    printf("#         memleak %d\n",lowest_count_in_hash());
    printf("#    on the command line when you run the program.\n#\n");
#endif
  }

  fprintf(stdout,"# Max bytes in use during this program = %d\n",Max_mallocked);
  fprintf(stdout,"# Total bytes currently in use         = %d\n",Total_mallocked-Total_freed);
  fprintf(stdout,"# Bytes currently on DAMUT free lists  = %d\n",free_bytes);
  fprintf(stdout,"#      Free nodes: %3d allocated; %3d are free themselves\n",
         Frees_allocated,
         free_length(Free_frees)
	 );

  if ( Num_zeroes_allocated > 0 )
    fprintf(stderr,
            "# You currently have %d zero-length byte%s \"allocated\"\n",
            Num_zeroes_allocated,
            (Num_zeroes_allocated == 1) ? "" : "s"
	    );

  for ( size = 0 ; size < MAX_AM_MALLOC_SIZE ; size++ )
  {
    if ( Started && Num_allocated[size] == 0 && Frees[size] != 0 )
      free_sizes_exist = TRUE;

    if ( Started && Num_allocated[size] > 0 )
    {
      boring = FALSE;
      fprintf(stdout,"# Size %4d bytes: %4d allocated; %4d on free list\n",
             size,
             Num_allocated[size],
             free_length(Frees[size])
	     );
    }
  }

  if ( free_sizes_exist )
  {
    int num_shown = 0;
    int n = 9;
    fprintf(stdout,"#\n# Following sizes have nothing allocated right now,\n");
    fprintf(stdout,"# but have free lists (length shown in parentheses)\n#\n");

    for ( size = 0 ; size < MAX_AM_MALLOC_SIZE ; size++ )
    {

      if ( Started && Frees[size] != NULL && Num_allocated[size] == 0 )
      {
        num_shown += 1;
        fprintf(stdout,"%s%d(%d)%s",
               (num_shown % n == 1) ? "# " : "",
               size,free_length(Frees[size]),
               (num_shown % n == 0) ? "\n" : " "
	       );
      }
    }

    if ( num_shown % n != 0 ) fprintf(stdout,"\n");
  }
  else if ( boring )
    fprintf(stdout,"# Nothing allocated, nor on any free list\n");
  
  fprintf(stdout,"#\n");
}

/* A simplified memory report for the PC.  

   Note that we must call close_statistics; otherwise memory will
   still be allocated for the statistics data structure. 
*/
char *pc_memory_report()
{
  int size;
  bool all_freed = TRUE;
  int free_bytes = 0;
  char s[150];
  char *report;
  void close_statistics();

  close_statistics(); 
  for ( size = 0 ; size < MAX_AM_MALLOC_SIZE ; size++ )
    free_bytes += ( sizeof(free_list) + size ) * free_length(Frees[size]);

  free_bytes += 2 * sizeof(free_list) * free_length(Free_frees);

  if ( Current_large_mallocked > 0 ) all_freed = FALSE;

  for ( size = 0 ; all_freed && size < MAX_AM_MALLOC_SIZE ; size++ )
    all_freed = Num_allocated[size] == 0;

  sprintf(s, "All Freed = %s, Total_mallocked (%d) - free_bytes (%d) - Total_freed (%d) = %d", 
	  (all_freed == TRUE) ? "True" : "False", Total_mallocked,
	  free_bytes, Total_freed, Total_mallocked-(free_bytes + Total_freed));
  report = mk_copy_string(s);
  return(report);
}

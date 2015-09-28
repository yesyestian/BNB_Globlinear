/* *
   File:         amma.h
   Author:       Andrew W. Moore
   Created:      4th September 1992
   Description:  My own memory allocation and deallocation.

   Copyright (C) Andrew W. Moore, 1992
*/

#ifndef AMMA_H

#define AMMA_H

#include "ambs.h"
#include "amstr.h"
#include <stdarg.h>

/*This block will ignore all the am_malloc stuff and handle all memory normally*/
#ifdef USE_OS_MEMORY_MANAGEMENT
#define AM_MALLOC(t) ((t *) malloc(sizeof(t)))
#define AM_MALLOC_ARRAY(t,len) ((t *) malloc((len) * sizeof(t)))

#define AM_FREE(thing,type) free(thing)
#define AM_FREE_ARRAY(thing,type,len) free(thing)

#define am_malloc(len) ((char*) malloc(len))
#define am_free(thing,len) free(thing)

#else

#define AM_MALLOC(t) ((t *) am_malloc(sizeof(t)))
#define AM_MALLOC_ARRAY(t,len) ((t *) am_malloc((len) * sizeof(t)))

#define AM_FREE(thing,type) (am_free((char *)(thing),sizeof(type)))
#define AM_FREE_ARRAY(thing,type,len) \
    ( am_free((char *)(thing) , sizeof(type) * (len)) )

#ifdef AMFAST
char *am_extended_malloc(int size);
void am_extended_free(char *memory, int size);
#define am_free(memory,size) am_extended_free(memory,size) 
#define am_malloc(size) am_extended_malloc(size)
#else
char *am_malloc(int size);
void am_free(char *memory, int size);
#endif

#endif /* ifdef USE_OS_MEMORY_MANAGEMENT */

#define MAX_AM_MALLOC_SIZE (1 << 12)

/* Call this function with am_malloc_number set to N if you are trying
   to debug a memory leak and if am_malloc_report() has warned you that
   on am_malloc call number N there was a leak. */
void memory_leak_stop(int am_malloc_number);

/* Call this function at the start of main(int argc,char *argv[]) thus:

     memory_leak_check_args(argc,argv)

   Then, if ever am_malloc_report tells you that you have a memory leak on
   the <N>th call the am_malloc, simply include 
      memleak <N>
   on the command line.
*/
void memory_leak_check_args(int argc,char *argv[]);

void am_free_to_os();

void reset_malloc_state();
extern void basic_am_malloc_report(void);
   /* If possible don't use this. Use am_malloc_report() from ammarep.c
      instead
   */

/* Declared but not defined - sir 8/6/2000
extern void am_malloc_test(void); */

/* #ifdef PC_MVIS_PLATFORM */
/* #include <stdarg.h> */
/* void eprintf ( FILE * file , char * lpszFormat, ...); */
/* #define fprintf eprintf */
/* #endif */

#endif /* #ifndef AMMA_H */

/*
   File:        ambs.h
   Author:      Andrew W. Moore
   Created:     Sat Sep 12 15:22:40 EDT 1992
   Updated:     June 25 97
   Description: Basic

   Copyright 1996, Schenley Park Research
*/

/*
ambs.c and ambs.h contain our set of basic utility functions--- the kind
of set of utilities that most C programming projects have. If you're working
on the Auton project feel free but not obliged to use them. Try to use
them in preference to your own version of the same thing if possible. If
something's missing, please feel free to add it, and then email awm+lab@cs
to tell us.

This file defines bool to be the boolean datatype, (with a typedef to int)
and #defines TRUE and FALSE to 1 and 0 respectively. It does a couple of
other convenient #defines too.
*/

#ifndef AMBS_H
#define AMBS_H

#include "standard.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef BOOL_DEFINED
#define bool int
#define BOOL_DEFINED
#endif

#ifndef DOUBLE_PTR_DEFINED
typedef double *double_ptr;
#define DOUBLE_PTR_DEFINED
#endif

#ifndef CHAR_PTR_DEFINED
typedef char *char_ptr;
#define CHAR_PTR_DEFINED
#endif

#ifndef CHAR_PTR_PTR_DEFINED
typedef char_ptr *char_ptr_ptr;
#define CHAR_PTR_PTR_DEFINED
#endif

#ifndef PI
#define PI 3.1415926535897932384626
#endif
#define ROOT_TWO 1.4142135623730952

#define EQ_PTR(p1,p2) (((long) p1) == ((long) p2))

extern double Verbosity;

extern bool htmlp;
/*
Use this if you wish to print out diagnostics in development mode, but not
in normal use. Conditionally print out if Verbosity is greater than a certain
value. Informally, the user should expect that with Verbosity 0.0, only
useful stuff appears, Verbosity = 1.0 adds a bit more, Verbosity = 10.0
prints a fair amount and occasionally even stops and waits for a key,
Verbosity = 100.0 prints all kinds of rubbish.
*/

extern bool SuppressMessages;

/* Output a message, and then wait for a key
*/
void my_info(char *string);

/* Signal a warning with the message, and then wait for a key
*/
void my_warning(char *string);

/* Signal an *** Auton Error: with the message, and then stop. Actually
    waits for a key before exitting so that in the debugger you get hit
    cntrl-c and examine the program stack
*/
void my_error(char *string);

/* like my_error, but accepts printf-like params and instead of 
   never returning a void, it never returns an int */
int my_errorf(const char *format, ...);

/* Waits for user to hit return before continuing. Also gives user
   the option of doing some other things (including turning off
   future wait_for_keys(). In a windows environment, MIGHT do a
   dialog box instead of a console print */
void wait_for_key(void);

/* Prints the message where the user will see it. In a windows
   environment might write this into a status window instead of
   the regular console. */
void status_report(char *message);

char *safe_malloc(unsigned size);
/*
  Very very basic malloc. Please don't use this. Use am_malloc defined
  in amma.h
*/

FILE *safe_fopen(char *fname, char *access);
/*
   Just like fopen, except signals a my_error if fopen fails
*/

bool eq_string(char *s1, char *s2);
bool eq_string_with_length(char *s1, char *s2, int n);
/*
   Returns TRUE if and only if the strings contain the same characters.
   The first reads up to a NULL character.  The 2nd reads the 1st n chars
*/

/* Uses the time to set a random seed for the random() functions above */
void am_randomize(void);

/* Following seeds the random number generator */
void am_srand(int seed);

void push_current_am_srand_state();
void pop_current_am_srand_state();

/*  Returns a random integer uniformly >= 0 and < n. */
int int_random(int n);

double range_random(double lo, double hi);
/*
   Returns a random double x uniformly s.t. lo <= x <= hi
*/

/* Generate a gaussian random number mean = 0, variance
   (= sdev) = 1.0 */
double gen_gauss();

/*
   Following functions have obvious properties. In some cases standard macro
   implementations exists. You can use these if you fear strange behavior
   of like to have things typechecked properly, or you fear they might
   not be implemented in all c standard libraries
*/

double real_min(double x, double y);
double real_max(double x, double y);
double real_abs(double x);
double real_square(double x);
double real_cube(double x);

int int_min(int x, int y);
int int_max(int x, int y);
int int_abs(int x);
int int_square(int x);
int int_cube(int x);

void wait_for_key(void);
void really_wait_for_key(void);
/*
  Use this if you want to pause the program to let the use see what's
  going on. Offers the user some other alternatives to hitting return,
  including the ability to ignore future wait_for_key()s
*/

char *input_string(char *mess, char *string, int string_size);
double input_realnum(char *mess);
int input_int(char *mess);
bool input_bool(char *mess);


/* The following functions are very useful from extracting user options
   from the command-line argc and argv.

   The first simply tells you whether a string appeared on the command
   line existed, and if so what its index in argv[] is. 

    index_of_arg(string,argc,argv) returns -1 if arg doesn't exist, returns
    index otherwise.

    string_from_args(key,argc,argv,default_string)
      searches for key on the command line. If it finds it, and
      theres another string to the right of it on the command line,
      returns that string (a pointer to the string...allocates no memory).
      If it doesn't find it, returns the default.

     {int,double,bool}_from_args are similar.

    EXAMPLE:

     {
       int size = int_from_args("size",argc,argv,20);
       double weight = double_from_args("weight",argc,argv,12.4);
     }

     if the program was ran with

       prog weight 404

     then after that code segment, size would be 20 and weight would be 404.0

  Details:
     These functions permit a leading dash in front of an argument, e.g.

       prog -weight 404

     would have had the same effect.

  
     If the string arghelp appears on the command line, then all these functions
     tell the user what key they expect, what type they are looking for, and
     what the default value is.
*/ 
    
int index_of_arg(char *opt_string, int argc, char *argv[]);
char *mk_string_from_args(char *key,int argc,char *argv[],char *default_value);
char *string_from_args(char *key, int argc, char *argv[], char *default_value);
double double_from_args(char *key, int argc, char *argv[], double default_value);
int int_from_args(char *key, int argc, char *argv[], int default_value);
bool bool_from_args(char *key, int argc, char *argv[], bool default_value);


int my_irint(double x);
/*
   rounds a double and turns it into an int.
   irint() the library function doesn't exists everywhere, so use this instead
*/


/*
   The following functions print m1 then the argumnet and m2 to the requested
   stream. What use are they? They can be useful in semi-automated contruction
   of fprintf functions for big user-defined structures. 
*/
void fprintf_int(FILE *s,char *m1,int x,char *m2);
void fprintf_realnum(FILE *s,char *m1,double x,char *m2);
void fprintf_float(FILE *s,char *m1,double x,char *m2);
void fprintf_double(FILE *s,char *m1,double x,char *m2);
void fprintf_bool(FILE *s,char *m1,bool x,char *m2);
void fprintf_string(FILE *s,char *m1,char *x,char *m2);

/*
   Trivial string functions
*/

bool char_is_in(char *s,char c);
int num_of_char_in_string(char *s, char c);
int index_of_char(char *s,char c);
int global_time();
FILE *am_fopen(char *filename,char *mode);

/**** sleep doesn't exist on the Windows NT platform.  We have written
      am_sleep to replace it (this works correctly on Unix
      platforms, but just waits for a key to pressed udner Windows
      NT).  Programmers should call am_sleep instead of sleep.
****/
void am_sleep(int secs);
void am_usleep(int microsecs);

bool is_a_number(char *string);
bool bool_from_string(char *s,bool *r_ok);
int int_from_string(char *s,bool *r_ok);
double double_from_string(char *s,bool *r_ok);

extern void sensible_limits(
    double xlo,
    double xhi,
    double *res_lo,
    double *res_hi,
    double *res_delta
  );

int next_highest_power_of_two(int n);

double roundest_number_between(double lo,double hi);

unsigned sleep(unsigned seconds);
int unlink(const char *filename);

void my_breakpoint();

/* Returns TRUE if and only if x is NaN or Inf (i.e. returns FALSE
   if and only if x is a completely legal number) */         
bool is_ill_defined(double x);

#endif /* AMBS_H */

/* **
   File:            ambs.c
   Author:          Andrew William Moore
   Created:         7th Feb 1990
   Updated:         8 Dec 96
   Description:     Standard andrew tools. These are trivial and tedious.

   Copyright 1996, Schenley Park Research.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <setjmp.h>
#include <stdarg.h>
#include "ambs.h"
#include "adgui.h"

#ifdef PC_MVIS_PLATFORM
extern char GetcBuf;
extern HWND hEditIn;
extern HANDLE KeyPress;
extern CRITICAL_SECTION critsect;
#endif

/* The VC++ version on the PCs uses a different (dialog/menu-based)
   interface, so does not include this code.... */
//#ifndef PC_MVIS_PLATFORM

int Ignore_next_n = 0;

void my_breakpoint()
{
#ifndef AMFAST
  wait_for_key();
  printf("I suggest you set a breakpoint at the function my_breakpoint()\n");
#endif
  /* skip */
}

bool SuppressMessages = FALSE;

/* Updated by AWM so that if comment_is_on it also prints to stdout
   in addition to popping up a message. Non-windows unchanged. */
void my_info(char *string)
{
  bool printed_to_stdout = FALSE;
#ifdef PC_MVIS_PLATFORM
  if(!SuppressMessages) MessageBox(hWnd,string,"Note",MB_ICONINFORMATION|MB_OK);
  else 
  {
    printf("%s\n",string);
    printed_to_stdout = TRUE;
  }
#else
  printf("%s\n",string);
  printed_to_stdout = TRUE;
#endif
  if ( !printed_to_stdout )
  {
    bool comment_is_on();
    if ( comment_is_on() )
      printf("NOTE: %s\n",string);
  }
}

/* Updated by AWM so that if comment_is_on it also prints to stdout
   in addition to popping up a message. Non-windows unchanged. */
void my_warning(char *string)
{
  bool printed_to_stdout = FALSE;
#ifdef PC_MVIS_PLATFORM
  if(!SuppressMessages) MessageBox(hWnd,string,"Warning",MB_ICONEXCLAMATION|MB_OK);
  else 
  {
    printf("Warning! %s\n",string);
    printed_to_stdout = TRUE;
  }
#else
  printf("Warning! %s\n",string);
  printed_to_stdout = TRUE;
  wait_for_key();
#endif
  if ( !printed_to_stdout )
  {
    bool comment_is_on();
    if ( comment_is_on() )
      printf("WARNING: %s\n",string);
  }
}

void my_error(char *string)
{
  /* AWM: I put this back so that I can see the printed output in
     the gray window before
     jumping into the debugger. */
  printf("***** SPR software error:\n***** %s\n",string);
  my_breakpoint();
#ifdef PC_MVIS_PLATFORM
  printf("Will now exit...\n");
  MessageBox(hWnd,string,"SPR software error",MB_ICONSTOP|MB_OK);
#else
  printf("Will now exit...\n");
  really_wait_for_key();
#endif
  exit(-1);
}

int my_errorf(const char *format, ...)
{
  /* AWM: I put this back so that I can see the printed output in
     the gray window before
     jumping into the debugger. */
  va_list ap;
  va_start(ap, format);
  printf("***** SPR software error:\n***** ");
  vfprintf(stdout, format, ap);
  putchar('\n');
  va_end(ap);
  my_breakpoint();
#ifdef PC_MVIS_PLATFORM
  printf("Will now exit...\n");
  MessageBox(hWnd,"Implemented this if you want","SPR software error",MB_ICONSTOP|MB_OK);
#else
  printf("Will now exit...\n");
  really_wait_for_key();
#endif
  exit(-1);
  return 0;						/* shut compiler up */
}

void wait_for_key(void)
{
  int c = ' ';
  int n = 0;

  if ( Ignore_next_n > 0 )
    Ignore_next_n -= 1;
  else if ( Verbosity >= 0.0 ) 
  {
    gui_prompt("[Hit Return] (or h for other options)",CLI_COMMAND_INPUT,FALSE);

    while ( c != '\n')
    {
      c = gui_getchar();

      if ( c >= '0' && c <= '9' )
        n = 10 * n + c - '0';
      if ( c == 'q' || c == 'Q' )
      {
        printf("You hit Q at a wait point. Exiting program\n");
#ifdef PC_MVIS_PLATFORM
        gui_unprompt();
        gui_quit();
#else
        exit(0);
#endif
      }
      if ( c == 'b' || c == 'B' )
      {
	my_breakpoint();
      }
      if ( c == 'h' || c == 'H' )
      {
        printf("Enter:\n");
        printf("   <Return>          To continue past this wait point.\n");
        printf("   <number> <Return> To continue past this and the next\n");
        printf("                     <number> wait points.\n");
        printf("   q (or Q)          To exit the program entirely.\n");
        printf("   b (or B)          Call my_breakpoint\n");
        printf("   h (or H)          To read this help menu.\n");
        printf("   v <number>        To change the Verbosity level. The\n");
        printf("                     higher the verbosity level the more\n");
        printf("                     information will be printed on your\n");
        printf("                     terminal (and the more wait points).\n");

#ifdef PC_MVIS_PLATFORM
        SetWindowText(hEditIn,NULL);
#else
        while ( c != '\n') c = getchar();
#endif
        c = ' ';
      }
      if ( c == 'v' || c == 'V' ) 
      {
        char buff[100];
        int i = 0;
        c = gui_getchar();
        while ( i < 99 && c != '\n')
        {
          buff[i] = c;
          buff[i+1] = '\0';
          i += 1;
          c = gui_getchar();
        }

        printf("Old Verbosity: %g\n",Verbosity);
        Verbosity = atof(buff);
        printf("New Verbosity: %g\n",Verbosity);
      }
    }
  
    Ignore_next_n = n;
    gui_unprompt();
  }
}

void status_report(char *message)
{
  printf("%s\n",message);
}

double Verbosity = 0.0;

char *safe_malloc(unsigned size)
{
  char *res;

  if ( size == 0 )
    res = NULL;      /* Doesn't matter what we return, but lets be tidy */
  else
  {
    res = (char *) malloc(size);
    if ( res == NULL )
      my_error("Safe malloc wasn't.\n");
  }
  return(res);
}

FILE *safe_fopen(char *fname, char *access)
{
  FILE *s = fopen(fname,access);
  if ( s == NULL )
  {
    printf("Could not fopen('%s','%s')\n",fname,access);
    my_error("safe_fopen() failed");
  }
  return(s);
}

bool eq_string(char *s1, char *s2)
{
  return(strcmp(s1,s2) == 0 );
}

bool eq_string_with_length(char *s1, char *s2, int n)
{
  return(strncmp(s1,s2,n) == 0 );
}

/*
   New Random Functions
   Author:      Darrell Kindred + a few mods by Justin Boyan
                 + tiny final adjustments by Andrew Moore
   Created:     Sun. Mar  3, 1996 16:22
   Description: the random number generator recommended by Persi Diaconis
                in his lecture at CMU, February 1996.

		The numbers (x_n) in the sequence are generated by
		  x_n = x_{n-PD_A} * x_{n-PD_B}   (mod 2^32)
		where PD_A=55, PD_B=24.  This gives a pseudorandom
		sequence of ODD integers in the interval [1, 2^32 - 1].

   Modified:    Mon Feb 17 13:42:19 EST 1997
                Scott identified hideous correlations in the low-order bit.
		As a temporary hack, I've killed the lowest 3 bits of
		each number generated.  This seems to help, but we should
		look into getting a more trustworthy RNG.  persi_diaconis--
		The old version is backed up in xdamut/backup/97-02-17.ambs.c
*/


#define PD_A 55
#define PD_B 24

typedef struct pd_state
{
  int cursor;
  unsigned int val[PD_A];
} pd_state;

pd_state Pd_state;
bool Pd_state_defined = FALSE;

bool Currently_saving_state = FALSE;
pd_state Saved_pd_state;

void copy_pd_state(pd_state *src,pd_state *dst)
{
  int i;
  dst -> cursor = src -> cursor;
  for ( i = 0 ; i < PD_A ; i++ )
    dst -> val[i] = src -> val[i];
}

void push_current_am_srand_state()
{
  if ( Currently_saving_state )
    my_error("save current am_srand_state: already being pushed. Implement a real stack.");
  copy_pd_state(&Pd_state,&Saved_pd_state);
  Currently_saving_state = TRUE;
}

void pop_current_am_srand_state()
{
  if ( !Currently_saving_state )
    my_error("pop_current am_srand_state: nothing pushed.");
  copy_pd_state(&Saved_pd_state,&Pd_state);
  Currently_saving_state = FALSE;
}

unsigned int pd_random_29bit() 
{
    int c = Pd_state.cursor;
    if (!Pd_state_defined) 
    {
      /* not seeded */
      am_srand(5297);
    }
    /* this assumes 32-bit unsigned int */
    Pd_state.val[c] *= Pd_state.val[(c + PD_A - PD_B) % PD_A];

    Pd_state.cursor = ((c + 1) % PD_A);

    return Pd_state.val[c] >> 3;
}

unsigned int pd_random_32bit() 
{
    unsigned int x = pd_random_29bit(), y = pd_random_29bit();
    return x ^ (y<<3);
}


/**********************************************************************
 * interface for the outside world...  (jab)
 **********************************************************************/


/* Pass in a single seed, which can be any unsigned int */
void am_srand(int seed)
{
  unsigned int us_seed = (unsigned int) seed & 0x3FFFFFFF;
  unsigned int i;
    /* this code requires 32-bit ints */
  if (sizeof(unsigned int) != 4) 
    my_error("this code requires 32-bit integers");

    /* make it odd */
  us_seed = us_seed * 2 + 1;
    /* add +0,+2,+4,... */

  for (i=0; i<PD_A; i++) {
    Pd_state.val[i] = us_seed + 2 * i;
  }
  Pd_state.cursor = 0;
  Pd_state_defined = TRUE;

    /* Throw out the first few numbers (enough to be sure that 
     * all values in the state have "wrapped around" at least twice). */
    for (i=0; i<400+(us_seed % 100); i++) {
	(void) pd_random_29bit();
    }
}

void am_randomize()
{
    long my_time = time(NULL);
    char *time_string = ctime(&my_time);
    int i = 0;
    int seed = 0;
    
    while ( time_string[i] != '\0' )  {
	seed += i * (int) (time_string[i]);
	i += 1;
    }
    am_srand(seed);
}

/* Synonmym for am_randomize. For backwards compat. Dont use */
void andrews_randomize()
{
  am_randomize();
}

double pd_random_double() 
{
    /* returns a double chosen uniformly from the interval [0,1) */
    return ((double) pd_random_29bit()) / ((unsigned int) (1<<29));
}


double basic_range_random(double lo, double hi)
{
    return lo + (hi-lo)*pd_random_double();
}

int rrcalls = 0;

double range_random(double lo, double hi)
{
  double result = basic_range_random(lo,hi);
//  if ( rrcalls < 10 )
//    printf("rrcall %d: range_random(%g,%g) = %g\n",rrcalls,lo,hi,result);
  rrcalls += 1;
  return result;
}

int int_random(int n)
{
  int result = -1;
  if ( n < 1 ) 
  {
    printf("int_random(%d) ?\n",n);
    my_error("You called int_random(<=0)");
    return(0);
  }
  else
    result = pd_random_29bit() % n;

  return(result);
}

int GG_iset = 0;
double GG_gset = -77.77;

static int Gcalls = 0;

double basic_gen_gauss()
{
  double fac = 0.0;
  double r = 0.0;
  double v1 = 0.0;
  double v2 = 0.0;

  double result = -777.777;

  if  (GG_iset == 0) 
  {
    bool ok = FALSE;
    while ( !ok )
    {
      v1=2.0*range_random(0.0,1.0)-1.0;
      v2=2.0*range_random(0.0,1.0)-1.0;
      r=v1*v1+v2*v2;
      ok = r < 1.0;
    }
    fac=sqrt(-2.0*log(r)/r);
    GG_gset=v1*fac;
    GG_iset=1;
    result = v2*fac;

    if (  result < -10000.0 ) 
      my_error("gen_gauss nuts");
  } 
  else 
  {
    GG_iset=0;
    return GG_gset;
  }

  return result;
}

double gen_gauss()
{
  double result = basic_gen_gauss();

  Gcalls += 1;
  return result;
}

/* Below function is for backwards compatibility only */
int very_random_int(int n)
{
  return(int_random(n));
}

double real_min(double x, double y)
{
  return( (x < y) ? x : y );
}

double real_max(double x, double y)
{
  return( (x < y) ? y : x );
}

double real_abs(double x)
{
  return( (x < 0.0) ? -x : x );
}

double real_square(double x)
{
  return( x * x );
}

double real_cube(double x)
{
  return( x * x * x );
}

int int_min(int x, int y)
{
  return( (x < y) ? x : y );
}

int int_max(int x, int y)
{
  return( (x < y) ? y : x );
}

int int_abs(int x)
{
  return( (x < 0.0) ? -x : x );
}

int int_square(int x)
{
  return( x * x );
}

int int_cube(int x)
{
  return( x * x * x );
}


void really_wait_for_key(void)
{
  int c = ' ';
#ifdef PC_MVIS_PLATFORM
  gui_prompt("[Return to continue or q to quit]",CLI_KEY_INPUT,FALSE);
  while ( c != '\n' )
  {
    c = gui_getchar();
#else
  printf("really_wait_for_key(): [Return to continue or q to quit]");
  while ( c != '\n' )
  {
    c = getchar();
#endif
    if ( c == 'q' || c == 'Q' )
    {
      printf("You hit Q at a wait point. Exiting program\n");
#ifdef PC_MVIS_PLATFORM
      gui_unprompt();
      gui_quit();
#else
      exit(0);
#endif
    }
  }
#ifdef PC_MVIS_PLATFORM
  gui_unprompt();
#endif
}

char *input_string(char *mess, char *string, int string_size)
{
  int i;
  int c = '\n';

  while ( c == ' ' || c == '\n' )
  {
    if ( c == '\n' )
      printf("%s",mess);
    c = getchar();
  }

  for ( i = 0 ; i < string_size-1 && c != ' ' && c != '\n' && c != EOF ; i++ )
  {
    string[i] = (char) c;
    c = getchar();
  }

  string[i] = '\0';
  return(string);
}

double input_realnum(char *mess)
{
  char double_name[100];
  bool happy = FALSE;
  double result;
  
  while ( !happy )
  {
    (void) input_string(mess,double_name,100);
    happy = 1 == sscanf(double_name,"%lf",&result);
  }

  return(result);
}

int input_int(char *mess)
{
  double x = input_realnum(mess);
  int result = my_irint(x);
  if ( real_abs(x-result) > 0.001 )
    printf("input_int(): Warning.. truncated your doubleing point number\n");
  return(result);
}

bool input_bool(char *mess)
{
  int c = '\n';
  int cig;
  while ( c != 'n' && c != 'N' && c != 'y' && c != 'Y' )
  {
    c = '\n';
    while ( c == ' ' || c == '\n' )
    {
      if ( c == '\n' )
        printf("%s",mess);
      c = getchar();
    }
    cig = c;
    while ( cig != '\n' && cig != ' ' )
      cig = getchar();
  }

  return(c == 'y' || c == 'Y');
}

int very_basic_index_of_arg(char *opt_string, int argc, char *argv[])
/* Searches for the string among the argv[]
   arguments. If it finds it, returns the index,
   in the argv[] array, of that following argument.
   Otherwise it returns -1.

   argv[] may contain NULL.
*/
{
  int i;
  int result = -1;
  for ( i = 0 ; result < 0 && i < argc ; i++ )
    if ( argv[i] != NULL && strcmp(argv[i],opt_string) == 0 )
      result = i;

  if ( Verbosity > 350.0 )
  {
    if ( result >= 0 )
      printf("index_of_arg(): Key `%s' from in argument number %d\n",
             argv[result],result
            );
    else
      printf("index_of_arg(): Key `%s' not found\n",opt_string);
  }

  return(result);
}

#define MAX_ARG_CHARS 200

int basic_index_of_arg(char *opt_string, int argc, char *argv[])
/* Searches for the string among the argv[]
   arguments. If it finds it, returns the index,
   in the argv[] array, of that following argument.
   Otherwise it returns -1.

   argv[] may contain NULL.

   It is INSENSITIVE TO THE INITIAL - sign in front of a key
*/
{
  int result;

  if ( Verbosity > 2000.0 )
  {
    int i;
    fprintf(stderr,"Welcome to basic_index_of_arg\n");
    fprintf(stderr,"opt_string = %s\n",
            (opt_string==NULL)?"NULL" : opt_string);
    fprintf(stderr,"argc = %d\n",argc);
    for ( i = 0 ; i < argc ; i++ )
      fprintf(stderr,"argv[%d] = %s\n",i,
              (argv[i]==NULL) ? "NULL" : argv[i]);
  }

  result = very_basic_index_of_arg(opt_string,argc,argv);
  if ( result < 0 )
  {
    int i;
    char other_key[MAX_ARG_CHARS+2];
    if ( opt_string[0] == '-' )
    {
      for ( i = 0 ; i < MAX_ARG_CHARS && opt_string[i+1] != '\0' ; i++ )
        other_key[i] = opt_string[i+1];
      other_key[i] = '\0';
    }
    else
    {
      other_key[0] = '-';
      for ( i = 0 ; i < MAX_ARG_CHARS && opt_string[i] != '\0' ; i++ )
        other_key[i+1] = opt_string[i];
      other_key[i+1] = '\0';
    }

    result = very_basic_index_of_arg(other_key,argc,argv);
  }

  return(result);
}
    
char *basic_string_from_args(char *key, int argc, char *argv[], char *default_value)
{
  int index = basic_index_of_arg(key,argc,argv);
  char *result;

  if ( index >= 0 && index < argc-1 )
    result = argv[index+1];
  else
    result = default_value;

  return(result);
}

/**** Now for public, talkative, ones ****/

bool arghelp_wanted(int argc, char *argv[])
{
  bool result = basic_index_of_arg("arghelp",argc,argv) >= 0 ||
                basic_index_of_arg("-arghelp",argc,argv) >= 0;
  return(result);
}

int index_of_arg(char *opt_string, int argc, char *argv[])
/* Searches for the string among the argv[]
   arguments. If it finds it, returns the index,
   in the argv[] array, of that following argument.
   Otherwise it returns -1.

   argv[] may contain NULL.
*/
{
  int result = basic_index_of_arg(opt_string,argc,argv);
  if ( arghelp_wanted(argc,argv) )
    printf("ARGV %21s %8s         %8s %s on commandline\n",
           opt_string,"","",
           (result >= 0) ? "   " : "not"
          );

  return(result);
}

char *mk_string_from_args(char *key,int argc,char *argv[],char *default_value)
{
  extern char *mk_copy_string(char *);
  return (mk_copy_string(string_from_args(key,argc,argv,default_value)));
}

char *string_from_args(char *key, int argc, char *argv[], char *default_value)
{
  char *result = basic_string_from_args(key,argc,argv,default_value);
  if ( arghelp_wanted(argc,argv) )
  {
    printf("ARGV %21s %8s default %8s ",
           key,"<string>",
           (default_value == NULL) ? "NULL" : default_value
          );
    if ( (default_value == NULL && result == NULL) || 
         (default_value != NULL && eq_string(default_value,result))
       )
      printf("\n");
    else
      printf("commandline %s\n",result);
  }

  return(result);
}

double double_from_args(char *key, int argc, char *argv[], double default_value)
{
  char *string = basic_string_from_args(key,argc,argv,(char *)NULL);
  double result = (string == NULL) ? default_value : atof(string);

  if ( arghelp_wanted(argc,argv) )
  {
    printf("ARGV %21s %8s default %8g ",
           key,"<double>",default_value
          );
    if ( string == NULL )
      printf("\n");
    else
      printf("commandline %g\n",result);
  }

  return(result);
}

int int_from_args(char *key, int argc, char *argv[], int default_value)
{
  char *string = basic_string_from_args(key,argc,argv,(char *)NULL);
  int result = (string == NULL) ? default_value : atoi(string);

  if ( arghelp_wanted(argc,argv) )
  {
    printf("ARGV %21s %8s default %8d ",
           key,"<int>",default_value
          );
    if ( string == NULL )
      printf("\n");
    else
      printf("commandline %d\n",result);
  }

  return(result);
}

bool bool_from_args(char *key, int argc, char *argv[], bool default_value)
/*
   Searches for two adjacent lexical items in argv in which the first 
   item equals the string 'key'.

   Returns TRUE if and only if the second item is a string beginning
   with any of the following characters: t T y Y 1
   (this is meant to accomodate any reasonable user idea of words for
    specifying boolean Truthhood)


   If doesn't find anything matching key then returns default_value.
*/
{
  char *string = basic_string_from_args(key,argc,argv,(char *)NULL);
  bool result = (string == NULL) ? default_value : 
                (string[0] == 't' || string[0] == 'T' ||
                 string[0] == 'y' || string[0] == 'Y' ||
                 string[0] == '1'
                );

  if ( arghelp_wanted(argc,argv) )
  {
    printf("ARGV %21s %8s default %8s ",
           key,"<bool>",(default_value) ? "TRUE" : "FALSE"
          );
    if ( string == NULL )
      printf("\n");
    else
      printf("commandline %s\n",(result)?"TRUE":"FALSE");
  }

  return(result);
}

int my_irint(double x)
/* Returns the closest integer to x */
{
  return((int) floor(x + 0.5));
}

void fprintf_int(FILE *s,char *m1,int x,char *m2)
{
  fprintf(s,"%s = %d%s",m1,x,m2);
}

void fprintf_realnum(FILE *s,char *m1,double x,char *m2)
{
  fprintf(s,"%s = %g%s",m1,x,m2);
}

void fprintf_float(FILE *s,char *m1,double x,char *m2)
{
  fprintf_realnum(s,m1,x,m2);
}

void fprintf_double(FILE *s,char *m1,double x,char *m2)
{
  fprintf_realnum(s,m1,x,m2);
}

void fprintf_bool(FILE *s,char *m1,bool x,char *m2)
{
  fprintf(s,"%s = %s%s",m1,(x)?"True":"False",m2);
}

void fprintf_string(FILE *s,char *m1,char *x,char *m2)
{
  fprintf(s,"%s = \"%s\"%s",m1,x,m2);
}

int index_of_char(char *s,char c)
/*
   Returns the least index i such that s[i] == c.
   If no such index exists in string, returns -1
*/
{
  int result = -1;

  if ( s != NULL )
  {
    int i;
    for ( i = 0 ; result < 0 && s[i] != '\0' ; i++ )
      if ( s[i] == c ) result = i;
  }

  return(result);
}

bool char_is_in(char *s,char c)
{
  bool result = index_of_char(s,c) >= 0;
  return(result);
}

/* This returns the number of occurences of character c in string s. 
   It is legal for s to be NULL (in which case the result will be 0).
   
   Added by Mary on 8 Dec 96.
*/
int num_of_char_in_string(char *s, char c)
{
  int res = 0;

  if (s != NULL)
  {
    int i;

    for (i = 0; i < (int) strlen(s); i++)
      if (s[i] == c) res++;
  }
  return(res);
}

int global_time()
{
  int result = (int) time((time_t *) NULL);
  return(result);
}


/* curr_date() and curr_time() added by jab, 3 Aug 95 */

char datestr[10], timestr[10];

char *curr_date(void)
{
  struct tm *thetime;
  time_t secs = time(0);
  thetime = localtime(&secs);
  sprintf(datestr, "%2d-%02d-%02d",
	  thetime->tm_year, 1+thetime->tm_mon, thetime->tm_mday);
  return datestr;
}

char *curr_time(void)
{
  struct tm *thetime;
  time_t secs = time(0);
  thetime = localtime(&secs);
  sprintf(timestr, "%02d:%02d:%02d",
	  thetime->tm_hour, thetime->tm_min, thetime->tm_sec);
  return timestr;
}

FILE *am_fopen(char *filename,char *mode)
{
  FILE *s = NULL;
  bool all_whitespace = TRUE;
  int i;

  if ( filename == NULL )
	  my_error("am_fopen: Called with filename == NULL");
  
  for ( i = 0 ; filename[i] != '\0' && all_whitespace ; i++ )
	all_whitespace = filename[i] == ' ';

  if ( !all_whitespace )
	s = fopen(filename,mode);
  return(s);
}

#ifdef PC_PLATFORM
#include <windows.h>
void am_sleep(int secs)
{
  int millisecs = 1000 * secs;
  Sleep(millisecs);
}
#else
#include <unistd.h> /* for sleep */
void am_sleep(int secs)
{
  sleep(secs);
}
#endif /* not UNIX_PLATFORM */

#ifdef PC_PLATFORM
#include <windows.h>
void am_usleep(int microsecs)
{
  int millisecs = microsecs/1000;
  Sleep(millisecs);
}
#else
#include <unistd.h> /* for usleep */
void am_usleep(int secs)
{
  usleep(secs);
}
#endif /* not UNIX_PLATFORM */


bool is_a_number(char *string)
/*
   uses a finite state machine with 8 states.
   These are the symbols in the FSM alphabet:
    S     + or -
    P     .
    D     0 or 1 or .. or 9
    E     e or E
    X     anything else

  In the following descriptions of states, any unmentioned
  symbol goes to an absorbing, non-accepting, error state
       State 1: S --> 2 , D --> 4 , P --> 3
       State 2: P --> 3 , D --> 4 
       State 3: D --> 5
    *  State 4: D --> 4 , E --> 6 , P --> 5
    *  State 5: D --> 5 , E --> 6
       State 6: D --> 8 , S --> 7
       State 7: D --> 8
    *  State 8: D --> 8
  The starred states are acceptors.
*/
{
  int state = 1;
  int i = 0;
  bool err_state = FALSE;
  bool result;

  while ( string[i] != '\0' && !err_state )
  {
    char c = string[i];
    char symbol = ( c == '.' ) ? 'P' :
                  ( c == 'e' || c == 'E' ) ? 'E' :
                  ( c == '+' || c == '-' ) ? 'S' :
                  ( c >= '0' && c <= '9' ) ? 'D' : 'X';
/*
    printf("i = %d , c = %c , symbol = %c , state = %d\n",i,c,symbol,state);
*/
    if ( state == 1 )
    {
      if ( symbol == 'S' ) state = 2;
      else if ( symbol == 'D' ) state = 4;
      else if ( symbol == 'P' ) state = 3;
      else err_state = TRUE;
    }
    else if ( state == 2 )
    {
      if      ( symbol == 'P' ) state = 3 ;
      else if ( symbol == 'D' ) state = 4 ;
      else err_state = TRUE;
    }
    else if ( state == 3 )
    {
      if      ( symbol == 'D' ) state = 5;
      else err_state = TRUE;
    }
    else if ( state == 4 )
    {
      if      ( symbol == 'D' ) state = 4 ;
      else if ( symbol == 'E' ) state = 6 ;
      else if ( symbol == 'P' ) state = 5;
      else err_state = TRUE;
    }
    else if ( state == 5 )
    {
      if      ( symbol == 'D' ) state = 5 ;
      else if ( symbol == 'E' ) state = 6;
      else err_state = TRUE;
    }
    else if ( state == 6 )
    {
      if      ( symbol == 'D' ) state = 8 ;
      else if ( symbol == 'S' ) state = 7;
      else err_state = TRUE;
    }
    else if ( state == 7 )
    {
      if      ( symbol == 'D' ) state = 8;
      else err_state = TRUE;
    }
    else if ( state == 8 )
    {
      if      ( symbol == 'D' ) state = 8;
      else err_state = TRUE;
    }

    i += 1;
  }

  result = !err_state && ( state==4 || state==5 || state==8 );
  return(result);
}

bool bool_from_string(char *s,bool *r_ok)
{
  bool result = FALSE;
  *r_ok = FALSE;

  if ( s != NULL && s[0] != '\0' )
  {
    char c = s[0];
    if ( c >= 'a' && c <= 'z' ) c += 'A' - 'a';
    if ( c == 'Y' || c == 'T' || c == '1' )
    {
      *r_ok = TRUE;
      result = TRUE;
    }
    if ( c == 'N' || c == 'F' || c == '0' )
    {
      *r_ok = TRUE;
      result = FALSE;
    }
  }

  return(result);
}

/*------------------------------------------------*/
int int_from_string(char *s,bool *r_ok)
{
  int result;
  if ( is_a_number(s) )
  {
    result = atoi(s);
    *r_ok = TRUE;
  }
  else
  {
    result = 0;
    *r_ok = FALSE;
  }
  return(result);
}

/*------------------------------------------------*/
double double_from_string(char *s,bool *r_ok)
{
  double result;
  if ( is_a_number(s) )
  {
    result = atof(s);
    *r_ok = TRUE;
  }
  else
  {
    result = 0;
    *r_ok = FALSE;
  }
  return(result);
}

void sensible_limits( double xlo, double xhi, double *res_lo, double *res_hi, double *res_delta )
{
  double scale,rel_hi,rel_delta;

  if ( xlo > xhi )
  {
    double temp = xlo; xlo = xhi; xhi = temp;
  }

  if ( xhi - xlo < 1e-5 )
  {
    double xmid = (xlo + xhi)/2.0;
    xlo = xmid - 1e-5;
    xhi = xmid + 1e-5;
  }
    
  scale = pow(10.0,ceil(-log10(xhi - xlo)));
  rel_hi = scale * (xhi - xlo);

  rel_delta = ( rel_hi < 1.5 ) ? 0.2 :
              ( rel_hi < 2.5 ) ? 0.5 :
              ( rel_hi < 5.0 ) ? 1.0 :
              2.0;

  *res_delta = rel_delta / scale;
  *res_lo = *res_delta * floor(xlo / *res_delta);
  *res_hi = *res_delta * ceil(xhi / *res_delta);
}

#if 0
int next_highest_power_of_two(int n)
{
  int result;

  if (n == 0)
    return (0);
  result = 1;
  n--;
  while (n)
  {
    result <<= 1;
    n >>= 1;
  }
  return result;
}
#endif

#if 1
/* This version is about 3 times as fast.  -jab  */
int next_highest_power_of_two(int n)
{
  int result;

  if (n<=0) return 0;
  if (n<=1) return 1;
  if (n<=2) return 2;
  if (n<=4) return 4;
  if (n<=8) return 8;
  if (n<=16) return 16;
  if (n<=32) return 32;
  if (n<=64) return 64;
  if (n<=128) return 128;
  if (n<=256) return 256;
  if (n<=512) return 512;
  if (n<=1024) return 1024;
  if (n<=2048) return 2048;
  if (n<=4096) return 4096;
  if (n<=8192) return 8192;
  if (n<=16384) return 16384;
  if (n<=32768) return 32768;
  
  result = 65536;
  n--;
  n >>= 16;
  
  while (n)
  {
    result <<= 1;
    n >>= 1;
  }
  return result;
}
#endif

double roundest_number_between(double lo,double hi)
{
  double d = hi - lo;
  double result = -77777.7;
  if ( d < 0.0 )
    my_error("iouadbciusabc");
  else if ( d < 1e-20 )
    result = (lo + hi)/2;
  else
  {
    double dist = 1.0;
    double factors[] = { 0.5 , 0.2 , 0.1 , 0.05 };
    int num_factors = sizeof(factors)/sizeof(double);
    int i;
    bool ok = FALSE;
    lo += 0.001 * d;
    hi -= 0.001 * d;
    while ( dist > d )
      dist *= 0.1;
    while ( dist < d )
      dist *= 10.0;
    /* dist is now a power of ten between 1 and 10 times the size of d */   
    for ( i = 0 ; !ok && i < num_factors ; i++ )
    {
      double scale = dist * factors[i];
      double a1 = (lo + hi)/2;
      double a2 = a1 / scale;
      double a3 = floor(0.5 + a2);
      double a4 = a3 * scale;
      if ( a4 >= lo && a4 <= hi )
      {
        ok = TRUE;
        result = a4;
      }
    }

    if ( !ok ) my_error("Logically MUST be okay");
  }
  return result;
}

/* Returns TRUE if and only if x is NaN or Inf (i.e. returns FALSE
   if and only if x is a completely legal number) */         
bool is_ill_defined(double x)
{
  return FALSE;
}

/* *
   File:        command.h
   Author:      Andrew W. Moore
   Created:     Mon Jun  2 19:58:44 EDT 1997
   Description: Basics for an extensible commandline interface 

   Copyright 1998, SPR
*/

#ifndef command_H
#define command_H

#include "amut.h"
#include "standard.h"
#include "genarray.h"

#define ARG_INT 0
#define ARG_REAL 1
#define ARG_STRING 2
#define ARG_BOOL 3
#define ARG_FILE_OPEN 4
#define ARG_FILE_SAVE 5
#define ARG_INT_LIST 6
#define ARG_REAL_LIST 7
#define ARG_STRING_LIST 8
#define ARG_STRING_SET 9

typedef struct command
{
  string_array *line;
  string_array *args;
  char *s;
  char *a1;
  char *a2;
  char *a3;
  char *a4;
  ivec *attnums;
} command;

command *mk_command_from_string_array(string_array *sa);
command *mk_command_from_string(char *s);
void rebuild_command_from_string_array(command *c,string_array *sa);

char *command_arg_ref(command *c,int arg);

/* Message may be NULL, which means it's ignored.

   If non-null, then if we end up prompting the user,
   we prompt them with the message instead of the key */
char *mk_string_from_command_with_message(command *c,char *key,char *message,
                                          char *default_string);

int command_num_tokens(command *c);
char *mk_string_from_command(command *c,char *key,char *default_string);
int int_from_command(command *c,char *key,int default_val,int lo,int hi);
bool bool_from_command(command *c,char *key,bool default_val);
double double_from_command(command *c,char *key,double default_val,
			   double lo,double hi);

/*HERE is the big daddy function that builds any arguments not specified in
  the command line through subtle use of the gui (if it exists), or text
  prompts (if it does not).  The arguments are as follows:
  1. command *c -- The command.  Any new arguments get inserted into it as 
    they are made.  However, c->attnums does NOT get rebuilt.  Therefore
    c->attnums should only be made AFTER all arguments have been made.
  2. int argnum -- The index into c->line at which the argument should be
    found.  If it has no specific place and is only specified by name, then
    argnum can be set to -1.
  3. char *argname -- The name of the argument.  This is used both to find
    the argument if it's not in its specified place or doesn't have one, and
    also to generate a prompt for it if it isn't there.
  4. int argtype -- Can be any of the constants at the top of this file.
    A LIST differs from a SET in that it can have repeats and is ordered.
  5. double argmin,argmax -- These specify the minimum and maximum values
    of arguments of type ARG_INT and ARG_REAL.  If argmin==argmax==0 then 
    they have no bounds.  
  6. double argdefault -- The default value of the argument.  For arguments of
    type ARG_STRING with argparams, (int)argdefault is an index into argparams.
  7. string_array *argparams -- If argtype==ARG_STRING, then argparams is a list
    of the possible values of the argument.  If it can be anything, set
    argparams = NULL.
  8. bool validation -- Should be TRUE if the argument absolutely has to fit
    into the constraints specified above.

  WARNING: This function will return NULL if the user presses the cancel button.
  For this reason (and because it always returns a pointer, even for a ARG_BOOL),
    it is easiest to free the argument by calling free_arg() below, as this does
    whatever is necessary to free the argument of the given type, and does nothing
    if the argument is NULL.
*/
void* mk_arg_from_command(command *c,int argnum,char *argname,int argtype,
                          double argmin,double argmax,double argdefault,
                          string_array *argparams,bool validation);
void free_arg(void *arg,int argtype);

bool is_ok_number_between(command *c,int arg,double lo,double hi);
bool arg_exists_ok(command *c,int arg);
bool arg_is(command *c,int argnum,char *name);

command *mk_command_from_line(char *prompt);
void free_command(command *c);

void basic_command(void *env,command *c);

/* You MUST put the following function at the end of your main.c
   before you call am_malloc_report
   if you want to avoid a memory-leak complaint */
void free_global_helper();

void command_malloc_report();

void generic_cli(void *env,char *prompt,
                 void (*generic_command)(void *env,command *c),
		 int argc,char *argv[]);

void generic_dispatch(void *env,command *c,char *name,
                      void (*generic_dispatch)(void *env,command *c,bool act));

void multiple_dispatch(void *env,command *c,char *name,
                       void (*generic_dispatch)(void *env,command *c,bool act),
                       void dispatch_command(void *env,command *c));

extern char *Help_string;
extern int Menu_ID;
extern command *Next_Command;
extern bool Stop_Commands;
extern int Num_Iterations;

/* Makes a new null-terminated string using character
   s[start] s[start+1] ... s[end-1]
   strlen(result) = end - start
*/
char *mk_substring(char *s,int start,int end);

/* Useful if you are parsing a command and want to look for a certain key
   and then remove it. Returns TRUE if it found the key. */
bool find_and_remove(string_array *sa,char *key);

/* string_has_suffix("plop.ps",".ps") == TRUE
   string_has_suffix("plop.ps","ps") == TRUE
   string_has_suffix("plops",".ps") == FALSE
   string_has_suffix("plops","ps") == TRUE
   string_has_suffix("plop.ps.hello",".ps") == FALSE */
bool string_has_suffix(char *string,char *suffix);

/* Finds least index such that
     string_has_suffix(sa[index],suffix)
   and returns -1 if no such string exists in sa */
int find_index_of_string_with_suffix(string_array *sa,char *suffix);

/* If there's a string in sa which has the given suffix, returns a
   copy of the leftmost (lowest indexed such string) and removes that
   one string from sa.

   If there's no such string, returns NULL and leaves sa unchanged */
char *mk_find_string_with_suffix_and_remove(string_array *sa,char *suffix);

/* These are to do with "comment sessions". If the comment session is on,
   it means that all graphs created with ag_on will be saved to
   specially named postrscript files. Also the output of the program
   will include some latex directives so that it can be pasted into a 
   latex document. */
bool comment_is_on();
char *comment_session_name();

#endif

#ifndef ADGUI_H
#define ADGUI_H 

#include "standard.h"

/*
   File:        adgui.h
   Author:      Daniel T. Halstead
   Created:     April 18 1999
   Updated:     June 22 1999
   Description: Basic

   Copyright 1996, Schenley Park Research
*/

/*This file contains the functionality needed for a Windows GUI shell to 
  (hopefully) all the SPR programs.  Execution flow enters at WinMain()
  (when the project is originally built as a WinApp), creates windows, and
  also creates a secondary thread.  The primary thread then handles Windows
  messaging, while the secondary begins its execution in the first main() 
  function it can find.

  Projects to use this GUI MUST BE built as Windows Applications instead of 
  Windows Console Applications.  And they ALSO MUST #define PC_MVIS_PLATFORM
  rather than PC_TTY_PLATFORM in utils/standard.h.  And they ALSO MUST 
  enable multi-threading in the project settings under the C/C++ tab.
  */

#ifdef PC_MVIS_PLATFORM
#ifdef __cplusplus
#include <afx.h>
#else
#include <windows.h>
#endif
#include <commctrl.h>
#include <commdlg.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "amstr.h"

/*MINER must be #define'd for all new miner projects*/
#ifdef MINER
#include "../utils/miner_resource.h"
#else
#ifdef UNIX_PLATFORM
#include "miner_resource.h"
#else
#include "utils/planit_resource.h"
#endif
#endif

#define BUFSIZE (1<<18)  /*The number of characters in the text-buffer*/
#define LINESIZE (1<<14)  /*The max number of characters in one line*/
#define COMHISTSIZE 32   /*The number of commands stored in the command history*/

/*These are the different states that the input window can be in.
  It can be waiting for a command, waiting for specifics of a parameter,
  waiting for a key press, or busy processing (not taking any input).*/
#define CLI_COMMAND_INPUT 0           /*Takes a word*/
#define CLI_COMMAND_OR_CLICK_INPUT 1  /*Takes a word or a click*/
#define CLI_SPEC_INPUT 2              /*Takes a possibly empty word*/
#define CLI_KEY_INPUT 3               /*Takes a keypress*/
#define CLI_PROCESSING 4              /*Takes nothing*/

#define WM_SHOWSCROLL WM_USER   /*Message to a window to show/hide its scrollbars*/

#ifdef PC_MVIS_PLATFORM

extern HWND hWnd;
extern HINSTANCE hInstance;
extern void *DlgInfo;

#define gui_dialog(id,proc) DialogBox(hInstance,MAKEINTRESOURCE(id),hWnd,proc)

#endif  /*PC_MVIS_PLATFORM*/

#ifndef UNIX_PLATFORM
/*These override the true printf and fprintf functions, so that anything
  printf'd to standard out may also get sent into the text-window buffer
  or a log file.*/
int fprintf(FILE *fil,const char *str,...);
int printf(const char *str,...);
#endif /* ! UNIX_PLATFORM */

/*These functions give a GUI prompt to the user, then set the
  input-window to wait for a specific user-response using wait_type.
  Set fullsize to TRUE if the prompt should take up the entire input 
  window, keeping the user from being able to type anything.*/
void gui_prompt(char *str,int wait_type,int fullsize);
/*Unprompt must be called after the user-input is recieved, in order
  to get rid of the prompt.*/
void gui_unprompt(void);

/*This functions just like getchar().  Note that '\r's are replaced by '\n's.*/
char gui_getchar(void);

/*This functions just like fgetc(s).  Note that '\r's are replaced by '\n's.*/
char gui_getc(FILE *s);

/*This will MAKE a string from the command-input.  Uses scanf if consolized*/
char* gui_getstring(FILE *s);

/*This flushes the output to the text-window (if fil!=stdout, it does nothing).  
  If this is not called, output does not get sent until it reads a '\n'.*/
void gui_flush(FILE *fil);

void gui_link_command_to_menu(char *name,int id);

/*Call this to gracefully exit the secondary thread, and display a happy
  "Quit." message.*/
void gui_quit(void);

/*This will send all stdout messages to a file as well as to the screen*/
void gui_log(char *filename);

#endif 

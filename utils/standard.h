/*
   File:        standard.h
   Author:      Andrew W. Moore
   Created:     Wed Jun  7 15:29:06 EDT 1995
   Description: Includer of standard libraries

   Copyright (C) 1995, Andrew W. Moore
*/

#ifndef STANDARD_H
#define STANDARD_H


#ifdef WIN32
/* Precisely one of the following constants must be #defined and the
   others must not be #defined. */
/* #define UNIX_TTY_PLATFORM */
 /*#define UNIX_XW_PLATFORM*/
/* #define PC_TTY_PLATFORM */
#define PC_MVIS_PLATFORM
#else
#define UNIX_XW_PLATFORM
#endif

/* 
   UNIX_TTY_PLATFORM -- this code will run under Unix and not produce
                        any graphics.  Any graphics calls will be
                        ignored, or possibly produce a printed warning
   UNIX_XW_PLATFORM  -- this code will run under Unix and is free to
                        produce graphics using damut/amgr.h graphics
                        routines
   PC_TTY_PLATFORM   -- this code will run on PCs (compiled by Visual C++)
                        and must not produce any graphics.  Any
			graphics calls will be ignored, or possibly
			produce a printed warning, or possibly pop up
                        graphics windows, though all control will remain
                        with the standard input.
   PC_MVIS_PLATFORM  -- this code will run on PCs, compiled by Visual
                        C++ as part of a single-document project.  It
                        requires the use of Mary's EYE files.  There
                        is no stdio/stderr output or input.  User
                        communication uses expos, apicts, and aform
			interface functions.  The programmer is free
			to produce graphics using damut/amgr.h graphics
                        routines.
*/

/* The following #defined constants are created as a function of which
   one of the above list is defined.  */
#ifdef UNIX_TTY_PLATFORM
#define UNIX_PLATFORM
#define TTY_PLATFORM
#endif /* UNIX_TTY_PLATFORM */

#ifdef PC_TTY_PLATFORM
#define TTY_PLATFORM
#define PC_PLATFORM
#endif /* PC_TTY_PLATFORM */

#ifdef UNIX_XW_PLATFORM
#define UNIX_PLATFORM
#endif /* UNIX_XW_PLATFORM */

#ifdef PC_MVIS_PLATFORM
#define PC_PLATFORM
#endif /* PC_MVIS_PLATFORM */

#ifdef PC_PLATFORM
   /* !!!!  Comment out the pragma statement on non-Visual C++ platforms.
   I could not put it within a #ifdef PC_MVIS_PLATFORM because xdamut
   and xambl no longer set this in their Build Settings. 

   The pragma statement should disable the level 4 warning: "unreferenced 
   inline function has been removed" */
#pragma warning( disable : 4514)
#endif

/* Directory separation string  e.g. "/" for Unix, "\" for DOS */
#ifdef UNIX_PLATFORM
#define DIRSEP "/"
#else
#define DIRSEP "\\"
#endif /* UNIX_PLATFORM */

#ifndef STDIO_H
#include <stdio.h>
#ifndef STDIO_H
#define STDIO_H
#endif /* STDIO_H inner */
#endif /* STDIO_H outer */

#ifndef STRING_H
#include <string.h>
#ifndef STRING_H
#define STRING_H
#endif /* STRING_H inner */
#endif /* STRING_H outer */

#ifndef STDLIB_H
#include <stdlib.h>
#ifndef STDLIB_H
#define STDLIB_H
#endif /* STDLIB_H inner */
#endif /* STDLIB_H outer */

#ifndef MATH_H
#include <math.h>
#ifndef MATH_H
#define MATH_H
#endif /* MATH_H inner */
#endif /* MATH_H outer */

#ifdef PC_PLATFORM
#ifdef NDEBUG
#define AMFAST
#endif
#endif

#endif /* standrd_h */

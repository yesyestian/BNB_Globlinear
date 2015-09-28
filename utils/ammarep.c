
/*
   File:        ammarep.c
   Author:      Andrew W. Moore
   Created:     Sun Jun 18 10:51:22 EDT 1995
   Description: Explains current malloc situation

   Copyright (C) 1995, Andrew W. Moore
*/

#include "ammarep.h"
#include "amma.h"
#include "amdym.h"
#include "amiv.h"
#include "stats.h"
#include "command.h"

void am_malloc_report()
{
  basic_am_malloc_report();
  dym_malloc_report();
  ivec_malloc_report();
  string_array_malloc_report();
  stats_malloc_report();
  command_malloc_report();
}



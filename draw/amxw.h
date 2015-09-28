

/*
       File:             amxw.h
       Author:           Mary Soon Lee
       Created:          12 Nov 95
       Modified:         23 Mar 96

       Description:     Standard andrew header file for graphics
*/

#ifndef AMXW_H
#define AMXW_H

/*#include <windows.h>*/
#include "./utils/ambs.h"
#include "amgr.h"
#include "./utils/adgui.h"

#define NUM_GRAY_SCALES 32

#define NUM_AG_COLORS 15

#define NUM_COLORS (NUM_GRAY_SCALES + NUM_AG_COLORS)

#define AG_BLACK 0
#define AG_WHITE 1
#define AG_GRAY 2
#define AG_DARKRED 3
#define AG_RED 4
#define AG_OLIVE 5
#define AG_YELLOW 6
#define AG_DARKGREEN 7
#define AG_GREEN 8
#define AG_BLUEGREEN 9
#define AG_CYAN 10
#define AG_DARKBLUE 11
#define AG_BLUE 12
#define AG_PURPLE 13
#define AG_MAGENTA 14
#define AG_GRAYSCALE_BASE (NUM_AG_COLORS)

#ifdef PC_MVIS_PLATFORM
void RedrawPict(HDC hdc);
#endif

int simple_spectrum_color(double frac);
double invert_simple_spectrum_color(int color);
int opposite_simple_spectrum_color(int color);

#endif /* #ifndef AMXW_H */


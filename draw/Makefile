# Makefile for dirlib ~/g/draw
# Created by awm on:  Thu Nov  4 20:51:20 EST 1999

here		= draw

includes	= ongr.h lingraphs.h lingraph.h fgraph.h drac.h apict.h amxw.h amgr.h

sources		= ongr.c lingraphs.c lingraph.c fgraph.c drac.c apict.c amxw.c amgr.c

private_sources = 

siblings	= utils 

t:=debug
include ../Make.common


#Run:  make w=python cmodules
#for Python bindings - sir 9/6/2000
ifeq (python, $w)
	cmodules          = ongr lingraphs lingraph apict amxw amgr #fgraph drac
	interface_headers = $(includes)
	include .g/pywrap/Make.python
endif	

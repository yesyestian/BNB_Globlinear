# Makefile for dirlib ~/g/utils
# Created by awm on:  Sat Feb 28 09:46:50 EST 1998
# Modified by sir on: Thu Jun  8 14:37:59 EDT 2000

here		= utils

includes	= svd.h stats.h standard.h saarray.h randomization.h miner_resource.h imat.h hash.h genarray.h command.h boolvec.h amut.h amstr.h ammarep.h amma.h amiv.h amdyv.h amdym.h amdmex.h ambs.h amar.h adgui.h StdAfx.h 

sources		= svd.c stats.c saarray.c randomization.c imat.c hash.c genarray.c command.c boolvec.c amstr.c ammarep.c amma.c amiv.c amdyv.c amdym.c amdmex.c ambs.c amar.c adgui.c

private_sources = 

siblings	= 

t:=debug
include ../Make.common




#Run:  make w=python cmodules
#for Python bindings - sir 9/6/2000
ifeq (python, $w)
	interface_headers = svd.h stats.h standard.h saarray.h randomization.h \
                            miner_resource.h imat.h hash.h genarray.h command.h \
                            boolvec.h amut.h amstr.h ammarep.h amma.h amiv.h   \
                            amdyv.h amdym.h amdmex.h ambs.h amar.h

	cmodules          = svd stats saarray randomization imat boolvec \
                            ammarep amma amiv amdyv amdym amdmex amar \
			   #adgui hash genarray command amstr ambs

	include ../pywrap/Make.python
endif


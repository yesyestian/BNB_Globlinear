# Makefile for dirlib ~/g/extra
# Created by awm on:  Thu Nov  4 20:51:25 EST 1999

here		= extra

includes	= qopt.h motutils.h hrect.h bayesreg.h 

sources		= qopt.c motutils.c hrect.c bayesreg.c 

private_sources = 

siblings	= draw utils 

t:=debug
include ../Make.common


#Run:  make w=python cmodules
#for Python bindings - sir 9/6/2000
ifeq (python, $w)
  cmodules          = qopt motutils hrect bayesreg
  interface_headers = $(includes)
  include ../pywrap/Make.python
endif	

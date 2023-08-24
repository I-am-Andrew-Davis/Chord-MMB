# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../../../Chombo/lib

makefiles += .._Chord_test_advanced

ebase = test_MMB_interpolation

LibNames = HOMappedConservation \
           HOMappedGeometry \
           MOLAMRTimeDependent \
           AMRTimeDependent \
           HOAMRTools \
           AMRTools \
           BoxTools

XTRALIBS = Chord

include $(CHOMBO_HOME)/mk/Make.test

# additional libraries
LIBFLAGS += $(LAPACKLIBS)

XTRACPPFLAGS += -I../../src -I../../src/support -I../../src/ibc -I../../src/turbmodels




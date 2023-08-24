# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../../../Chombo/lib

makefiles += .._Chord_test_basic

ebase = test_CRDutil_Laplacian \
        test_MMB_fluxAverage \
        test_MMB_meshRefine \
        test_NVector_LDFAB \
        test_PatchMappedFunc_FaceAvgNormDerivFromCellAvg \
        test_PatchMappedFunc_FaceAvgValFromCellAvg \
        test_PatchMappedFunc_CellAvgDerivFromFaceAvg

LibNames := HOMappedConservation \
            HOMappedGeometry \
            MOLAMRTimeDependent \
            AMRTimeDependent \
            HOAMRTools \
            AMRTools \
            BoxTools

XTRALIBS = Chord

include $(CHOMBO_HOME)/mk/Make.test

# additional libraries
LIBFLAGS += $(LAPACKLIBS) $(FFTWLIBS)

XTRACPPFLAGS += -I../../src -I../../src/support -I../../src/turbmodels

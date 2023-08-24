# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += _.._Chord_exec

# the base name(s) of the application(s) in this directory
ebase = chord

# the location of Chombo lib dir
CHOMBO_HOME =../../Chombo/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames := HOMappedConservation \
            HOMappedGeometry \
            MOLAMRTimeDependent \
            AMRElliptic \
            AMRTimeDependent \
            HOAMRTools \
            AMRTools \
            BoxTools

# relative paths to source code directories
base_dir = .

src_dirs = ../src ../src/support ../src/ibc ../src/turbmodels ../src/lbfgsb

# input file for 'run' target
INPUT = chord_test.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# additional libraries
LIBFLAGS += $(LAPACKLIBS) $(FFTWLIBS)
# application-specific variables

# application-specific targets
.PHONY: doc doxygen

# application-specific rules
doc doxygen:
	mkdir -p ../doc/code
	$(DOXYGEN) ../doc/Chord.dox.conf | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)'
	@echo "" ; echo "Point browser at ../doc/code/html/index.html" ; echo ""

ifeq (${CTHR},TRUE)
  XTRACXXFLAGS += -DUSE_PROFILER
endif

#
# user build settings
#

# Clang debug build
#CXX      = clang++
#CXXFLAGS = -O0 -Wall -Wno-unused -g -DCPLUSPLUS

# GCC production build
CXX      = g++
CXXFLAGS = -O3 -funroll-all-loops -Wall -Wno-unused -Wno-write-strings -DCPLUSPLUS


#
# automatic settings and build targets
#

# check that we have GROMACS settings
ifndef GMXBIN
    $(error The "GMXPREFIX" environment variable is not defined. Source GMXRC from the GROMACS installation you want to use)
endif

# GMXPREFIX is not exported by GMXRC, for some reason, reconstruct it here.
GMXPREFIX = $(realpath $(GMXBIN)/..)

INCLUDES = -I$(GMXPREFIX)/include
LDFLAGS  = -L$(GMXLDLIB) -Wl,-rpath,$(GMXLDLIB)
LIBS     = -lgromacs -lm

all: bin/ctcf bin/or-df

bin/or-df: src/or-df.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) -o bin/or-df src/or-df.cpp $(LIBS)

bin/ctcf: src/ctcf.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) -o bin/ctcf src/ctcf.cpp $(LIBS)

clean:
	rm -f bin/ctcf bin/or-df

SRCDIR=./
OUTPUTDIR=../singleRun/
CC=g++
MATLIB=-lm -lmat -lmx
PYLIB=-lpython2.7
LBOOST=-lboost_program_options
OPT=-msse4.2 -march=westmere -fPIC
DEBUG=-g -fmax-errors=2 -fno-omit-frame-pointer
CXXFLAG=-std=c++11 -fPIC $(OPT)

# to do gain curve

## singleRun
sR : nsynBeta.cpp
	$(CC) $(CXXFLAG) $^ $(MATLIB) $(PYLIB) $(LBOOST) -o $(OUTPUTDIR)$@

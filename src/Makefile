SRCDIR=./
#OUTPUTDIR=../singleRun/
OUTPUTDIR=../gainCurve/
CC=g++
RPATH=-Wl,-rpath,/share/apps/NYUAD/matlab/R2015b/bin/glnxa64 -Wl,-rpath,/home/wd554/.conda/envs/py2/lib
MATLAB=-L/share/apps/NYUAD/matlab/R2015b/bin/glnxa64 -lmat -lmx
PYLIB=-L/home/wd554/.conda/envs/py2/lib -lpython2.7
#MATLAB=-lmat -lmx
PYLIB=-lpython2.7
LBOOST=-lboost_program_options
OPT=-msse4.2 -march=westmere -fPIC -fno-omit-frame-pointer -O2
DEBUG=-g -fmax-errors=2 -fno-omit-frame-pointer
CXXFLAG=-std=c++11 $(DEBUG)
OBJ=nNeuroSt.o nNeuroLib.o jumpy_linear.o jumpy_bilinear.o linear_bilinear.o input_args_CA1.o Yale_NEURON_PyAPI.o input_args.o
DEP=matFunc.h typedefs.h heap.h

# to do gain curve
one: $(OBJ) oneTone.o $(DEP)
	$(CC) $(CXXFLAG) $(OBJ) oneTone.o -o $(OUTPUTDIR)$@ $(RPATH) $(PYLIB) $(MATLAB) $(LBOOST)

gainCurve: $(OBJ) variableInput.o $(DEP)
	$(CC) $(CXXFLAG) $(OBJ) variableInput.o -o $(OUTPUTDIR)$@ $(RPATH) $(PYLIB) $(MATLAB) $(LBOOST)

nNeuroSt.o: nNeuroSt.cpp nNeuroSt.h heap.h
	$(CC) $(CXXFLAG) -c $< -o $@

nNeuroLib.o: nNeuroLib.cpp nNeuroLib.h matFunc.h
	$(CC) $(CXXFLAG) -c $< -o $@

jumpy_linear.o: jumpy_linear.cpp jumpy_linear.h jumpy_bilinear.h nNeuroLib.h nNeuroSt.h
	$(CC) $(CXXFLAG) -c $< -o $@

jumpy_bilinear.o: jumpy_bilinear.cpp jumpy_bilinear.h nNeuroLib.h nNeuroSt.h
	$(CC) $(CXXFLAG) -c $< -o $@

linear_bilinear.o: linear_bilinear.cpp linear_bilinear.h nNeuroLib.h nNeuroSt.h
	$(CC) $(CXXFLAG) -c $< -o $@

%.o: %.cpp %.h
	$(CC) $(CXXFLAG) -c $< -o $@

clean:
	rm $(OBJ) variableInput.o oneTone.o

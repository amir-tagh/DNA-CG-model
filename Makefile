#which compiler to use
#CXX=g++
CXX=g++ 
#it is better to add -march flag for specific CPU of the computer like -march=native 
# -O1 level of optimization is the most basic one and -O2 is the best recommended one -Og:enables optimization that do not intefere with debugging 
# -g produces debugging information 
#Wall this flag enables all the warnings about constructions
#wno-deprecated:do not warn about functions,variables and types  
#-I add the directory to the head of the list of directories to be searched for header files.
#-O3 WAS ADDED AT THE END OF CXXFLAGS
CXXFLAGS=-Wall -Wextra -Wno-deprecated -I$(HOME) -DROTORBOX -D_WANGLANDAU -D_SAK_NOBREAK -D_NDEBUG -D_EXT_ADAPTIVE -DINTERCAL  -O3 
#add -pg for profling
PROFLAGS=-Wall  -C  $(CXXFLAGS)

VALGRIND_FLAGS=-Wall -g -DROTORBOX -D_EXT_ADAPTIVE -DINTERCAL

CCFLAGS=$(CXXFLAGS)
CCFLAGS=$(PROFLAGS)
#CCFLAGS=$(VALGRIND_FLAGS)


##original library of the makefile is -lm
LIBS=-llapack -lblas -lgfortran  
#add -pg in linking for profiling
LDFLAGS= -O3
OBJECTS=ran2.o inout.o cell.o hashTable.o intercalate.o switchS.o setup.o FreeEnergyHS.o bond.o  tools.o  dihed.o LJ.o displaceparticle.o sak_potential.o rigidbodyrotation.o topologyrotation.o  LJtopology.o neighbour_debug.o rigidbodydihed.o rigidbond.o dihedtopology.o bondtopology.o Shrink.o Extension.o   

HEADERS=cell.h  FreeEnergyHS.h  hashTable.h  LJ.h  sak_potential.h  tools.h wangLandau.h
 

#########################
#Rules: target: source
#		command
#########################

#.SUFFIXES: .cc .c .C
all: FreeEnergyHS $(HEADERS) 

# $@ and $< are called automatic variables the first is the output variable and the second is called the input variable. -c flag is to generate the .o file (output file) -o means the output.
.C.o: $(HEADERS)
	$(CXX) -c $(CCFLAGS) -o $@ $<


FreeEnergyHS: $(OBJECTS)  
	$(CXX) $(LDFLAGS) $(LIBS) $^ -o FreeEnergyHS.exe

.PHONY: clean
clean: 
	rm -f $(OBJECTS) *~ \#*

	
#########
#gprof options
#########
# -C:to print a tally of functions and the number of times each was called.
# -q:this option is to print the call graph analysis.

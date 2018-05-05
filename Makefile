# Makefile for Hamiltonian Heuristic
# Author: Michel Wan Der Maas
# 05/15/16


CXX = g++ #must support C++11

# sources in alphabetical order
SOURCES = BinaryHeap.cpp
SOURCES += Graph.cpp
SOURCES += Matching.cpp
SOURCES += MST.cpp

CXXFLAGS = -std=gnu++11 -Wall
debugfile = hamiltonian_debug

debug: CXXFLAGS += -O3

other: CXXFLAGS += -O0 -g3 -pg
other: EXEC = hamiltonian_other
other: SOURCES += other_main.cpp
    
linear: CXXFLAGS += -O0 -g3 -pg
linear: EXEC = hamiltonian_linear
linear: SOURCES += linear_main.cpp

OBJECTS = $(SOURCES:.cpp=.o)

other: $(OBJECTS) other_main.o
	${CXX} $(CXXFLAGS) $(OBJECTS) -o $(EXEC) $(LNFLAGS)

linear: $(OBJECTS) linear_main.o
	${CXX} $(CXXFLAGS) $(OBJECTS) -o $(EXEC) $(LNFLAGS)

debug: $(OBJECTS)
	${CXX} $(DBGFLAGS) $(CXXFLAGS) $(OBJECTS) -o $(debugfile) $(LNFLAGS) 

clean:
	$(RM) *.o *~
	#$(RM) $(EXEC)
	$(RM) $(debugfile)

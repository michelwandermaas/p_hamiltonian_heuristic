# Makefile for Hamiltonian Heuristic
# Author: Michel Wan Der Maas
# 05/15/16


CXX = g++-6 #must support C++11

# sources in alphabetical order
SOURCES = BinaryHeap.cpp
SOURCES += Graph.cpp
SOURCES += main.cpp
SOURCES += Matching.cpp
SOURCES += MST.cpp
# SOURCES += MST_main.cpp

OBJECTS = $(SOURCES:.cpp=.o)

CXXFLAGS = -std=gnu++11 -Wall
EXEC = hamiltonian
debugfile = hamiltonian_debug

debug: CXXFLAGS += -O0 -g3 -pg
all: CXXFLAGS += -O3

all: $(OBJECTS)
	${CXX} $(CXXFLAGS) $(OBJECTS) -o $(EXEC) $(LNFLAGS)

debug: $(OBJECTS)
	${CXX} $(DBGFLAGS) $(CXXFLAGS) $(OBJECTS) -o $(debugfile) $(LNFLAGS) 

clean:
	$(RM) *.o *~
	$(RM) $(EXEC)
	$(RM) $(debugfile)

//Reads TSPLIB instances (Only symmetric and asymmetric TSP instances)

#pragma once

#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
using namespace std;

class Graph
{
public:
	Graph();
	Graph(char *filename);
	~Graph();

	void Read(char *filename);

	//Adjacency matrix
	int** M;
	double *X, *Y;
	int n;//number of vertices
};







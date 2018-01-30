#ifndef MST_H
#define MST_H

#include <math.h>
#include <stdio.h>

class MST
{
public:
	MST(int n, int m, int * uE, int * vE, double * c);
	~MST();


	int n;
	int m;
	int * uE;
	int * vE;
	double * c;
	
	int* sol;
	double cost;

	int * parent;
	int * rank;

	int * aux;
	int * sorted;

	void MakeSet(int x);
	void Union(int x, int y);
	int Find(int x);
	void SolveKruskal();
	void MergeSort(int a, int b);
	void Merge(int p, int q, int r);
};

#endif

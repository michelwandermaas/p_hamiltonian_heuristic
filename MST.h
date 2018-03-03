#ifndef MST_H
#define MST_H

#include <math.h>
#include <stdio.h>

typedef struct Edge{
    int v1;
    int v2;
    int cost;
    friend bool operator==(const Edge& e1, const Edge& e2){
	return e1.v1 == e2.v1 && e1.v2 == e2.v2;
    }
    friend bool operator<(const Edge& e1, const Edge& e2){
	return e1.cost < e2.cost;
    }
    friend bool operator>(const Edge& e1, const Edge& e2){
	return e1.cost > e2.cost;
    }
}Edge;

typedef struct{
    Edge* edges;
    double numEdges;
}Edges;

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

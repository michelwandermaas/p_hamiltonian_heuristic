#pragma once

#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BinaryHeap.h"
#include <list>
using namespace std;

#define EVEN 2
#define ODD 1
#define UNLABELED 0

class Matching
{
public:
	Matching();
	Matching(int n);
	void Init();
	~Matching();

	void AddEdge(int u, int v);
	void AddEdge(int u, int v, double c);
	void Grow();
	void Expand(int u);
	void Expand2(int u, int p, int q);
	void OutermostBlocked(int u, int &v);
	void Augment(int u, int v);
	void Reset();
	int GetFreeIndex();
	int Blossom(int u, int v);
	void UpdateDualCosts();
	void SolveMinimumCostPerfectMatching();
	void SolvePerfectMatching();
	void Clear();
	void DestroyBlossom(int t);
	void Open(int u);
	void Open2(int u, int p, int q);
	void SetCost(int u, int v, double c);
	void PrintNeato();
	void PrintNeato2();
	void GrantFeasibility();
	void Heuristic();
	void Heuristic2();
	void PositiveCosts();
	void DeleteEdges();
	int IsInMatching(int u, int v);
	Matching* getClone();

	int *Free;//List of free indices
	int sizeFree;//Size of the list

	int *blossom;//blossom[v] gives the index of the blossom where v is immediatelly contained (default is blossom[v] = v);
	int *outer;//outer[v] gives the index of the blossom that contains v but is not contained in any other blossom (default is outer[v] = v)
	int **deep;//deep[v] is a list of all the original vertices contained inside v
	int *sizeDeep;
	int **shallow;//shallow[v] is a list of the vertices immediately contained inside v, the list has the exact order of the odd circuit of the blossom
	int *sizeShallow;
	int *tip;//tip of the blossom 	
	int *active;

	int *type;//Even, odd, neither (2, 1, 0)
	int *forest;//forest[v] gives the father of v in the alternating forest
	int *root;//root[v] gives the root of the alternating forest 

	int *blocked;//A blossom can be blocked, this means that it behaves as if it were an original vertex and cannot be expanded
	double *dual;//dual multipliers associated to the blossoms, if dual[v] > 0, the blossom is blocked and full
	double *slack;//slack associated to each edge, if slack[e] > 0, the edge cannot be used
	int *mate;//mate[v] gives the mate of v
	int *matching;
	double obj;
	int lastInserted;

	int m, n;
	int *E;//Active edges
	int sizeE;
	int *E1, *E2;//Lists of edges (endpoints)
	double *cost;
	double maxCost;//Maximum possible cost
	double minEdge;

	int feasible;

	int *auxVertexArray1, *auxVertexArray2;

	int perfect;

	int **AdjMat;
	int **Artificial;
	int **I;
	double **C;
	int **AdjList;
	int *sizeAdjList;

	list<int> BFSList;
	int *inList;
	int *visited;

	bool init;

	BinaryHeap* Bheap;
	BHNode **bhnodes;
};


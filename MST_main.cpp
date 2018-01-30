#include "MST.h"
#include <stdio.h>
#include <iostream>
using namespace std;

int main()
{
	//número de vértices
	int n = 5;
	//número de arestas
	int m = 7;

	//arestas: {0,1}, {0,2}, {0,4}, {1,3}, {1,4}, {2,3}, {2,4}
	//vetor de endpoints
	int uE[] = {0, 0, 0, 1, 1, 2, 2};
	//vetor de endpoints
	int vE[] = {1, 2, 4, 3, 4, 3, 4};

	//custos
	double c[] = {10, 5, 2, 4, 1, 3, 1};

	MST mst(n, m, uE, vE, c);
	mst.SolveKruskal();
	printf("------\n");
	for(int i = 0; i < n-1; i++)
		printf("%d %d\n", uE[ mst.sol[i] ], vE[ mst.sol[i] ]);

}

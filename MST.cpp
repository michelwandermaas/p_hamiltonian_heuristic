#include "MST.h"

MST::MST(int n, int m, int * uE, int * vE, double * c)
{
	this->n = n;
	this->m = m;
	this->uE = uE;
	this->vE = vE;
	this->c = c;
	
	sol = new int[n];

	parent = new int[n];
	rank = new int[n];

	aux = new int[m];
	sorted = new int[m];
}

MST::~MST()
{
	delete [] parent;
	delete [] rank;
	delete [] sol;
	delete [] aux;
	delete [] sorted;
}

void MST::MakeSet(int x)
{
	parent[x] = x;
	rank[x] = 0;
};	

void MST::Union(int x, int y)
{
	int xRoot = Find(x);
	int yRoot = Find(y);

	if(xRoot == yRoot) return;

	if(rank[ xRoot ] < rank[ yRoot ]){ parent[ xRoot ] = yRoot; }
	else if(rank[ xRoot ] > rank[ yRoot ]){ parent[ yRoot ] = xRoot; }
	else { parent[ yRoot ] = xRoot; rank[ xRoot ] += 1; }
};

int MST::Find(int x)
{
	if(parent[x] != x)
		parent[x] = Find(parent[x]);
	return parent[x];
};


void MST::SolveKruskal()
{
	for(int i = 0; i < m; i++)
		sorted[i] = i;
	MergeSort(0, m-1);

	for(int i = 0; i < n; i++)
		MakeSet(i);

	int size = 0;
	cost = 0;

	for(int i = 0; i < m; i++)
	{
		int u = uE[ sorted[i] ];
		int v = vE[ sorted[i] ];
		//printf("%d %d\n", u, v);
		if(Find(u) != Find(v))
		{
			Union(u,v);
			sol[size++] = sorted[i];
			cost += c[sorted[i]];
		}
	}
}

void MST::MergeSort(int a, int b)
{
	if(a==b) return;

	int middle = a + (int)( floor( (b-a)/2.0 ) );

	MergeSort(a, middle);
	MergeSort(middle+1, b);
	
	Merge(a, middle, b);
}

void MST::Merge(int p, int q, int r)
{
	for(int i = p; i <= r; i++)
		aux[i] = sorted[i];

	int pos = p;
	int p1 = p;
	int p2 = q+1;	

	while(p1 <= q || p2 <= r)
	{
		if(p2 > r || (p1  <= q && c[aux[p2]] - c[aux[p1]] > 0.00001 ))
		{
			sorted[pos] = aux[p1];
			p1++;
		}
		else
		{
			sorted[pos] = aux[p2];
			p2++;
		}
		pos++;
	}
}


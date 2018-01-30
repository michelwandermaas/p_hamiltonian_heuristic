#include "Matching.h"

Matching::Matching()
{
	Free = NULL;
	blossom = NULL;
	outer = NULL;
	deep = NULL;
	shallow = NULL;
	sizeDeep = NULL;
	sizeShallow = NULL;
	type = NULL;
	forest = NULL;
	root = NULL;
	blocked = NULL;
	dual = NULL;
	slack = NULL;
	mate = NULL;
	E1 = NULL;
	E2 = NULL;
	E = NULL;
	
	init = false;
}

Matching::Matching(int n)
{
	this->n = n;
	Init();
}

void Matching::Init()
{
	init = true;

	Free = new int[2*n];
	blossom = new int[2*n];
	outer = new int[2*n];
	deep = new int*[2*n];
	sizeDeep = new int[2*n];
	shallow = new int*[2*n];
	sizeShallow = new int[2*n];
	type = new int[2*n];
	forest = new int[2*n];
	root = new int[2*n];
	blocked = new int[2*n];
	mate = new int[2*n];
	dual = new double[2*n];
	tip = new int[2*n];
	AdjMat = new int*[n];
	Artificial = new int*[n];
	AdjList = new int*[n];
	sizeAdjList = new int[n];
	I = new int*[n];
	C = new double*[n];
	active = new int[2*n];
	visited = new int[2*n];
	bhnodes = new BHNode *[n];

	for(int i = 0; i < 2*n; i++)
	{
		deep[i] = new int[2*n];	
		shallow[i] = new int[2*n];

		if(i < n)
		{
			bhnodes[i] = new BHNode;
			bhnodes[i]->p = i;

			sizeAdjList[i] = 0;
			AdjMat[i] = new int[n];
			Artificial[i] = new int[n];
			AdjList[i] = new int[n];
			I[i] = new int[n];
			C[i] = new double[n];
			for(int j = 0; j < n; j++)
			{
				I[i][j] = -1;
				C[i][j] = 0;
				AdjMat[i][j] = 0;
				Artificial[i][j] = 0;
			}
		}
	}

	E1 = new int[(n*(n-1))/2];
	E2 = new int[(n*(n-1))/2];
	E = new int[(n*(n-1))/2];
	matching = new int[(n*(n-1))/2];
	slack = new double[(n*(n-1))/2];
	cost = new double[(n*(n-1))/2];
	
	for(int i = 0; i < (n*(n-1))/2; i++)
	{
		cost[i] = 0;
	}
	m = 0;	

	auxVertexArray1 = new int[2*n];
	auxVertexArray2 = new int[2*n];

	Bheap = new BinaryHeap(n+1);
	
	Clear();

	lastInserted = -1;
}

Matching::~Matching()
{
	if(!init) return;

	delete [] Free;
	delete [] blossom;
	delete [] outer;
	delete [] sizeDeep;
	delete [] sizeShallow;
	delete [] type;
	delete [] forest;
	delete [] root;
	delete [] blocked;
	delete [] mate;
	delete [] dual;
	delete [] tip;
	delete [] active;
	delete [] visited;

	for(int i = 0; i < 2*n; i++)
	{
		delete [] deep[i];	
		delete [] shallow[i];
		if(i < n)
		{
			delete bhnodes[i];
			delete [] AdjList[i];
			delete [] AdjMat[i];
			delete [] Artificial[i];
			delete [] C[i];		
			delete [] I[i];
		}
	}
	delete [] deep;
	delete [] shallow;
	delete [] AdjMat;
	delete [] Artificial;
	delete [] AdjList;
	delete [] sizeAdjList;
	delete [] C;
	delete [] I;
	delete [] bhnodes;

	delete [] E1;
	delete [] E2;
	delete [] E;
	delete [] slack;
	delete [] cost;
	delete [] matching;

	delete [] auxVertexArray1;
	delete [] auxVertexArray2;
	delete Bheap;
}

void Matching::AddEdge(int u, int v)
{
	E1[m] = u;
	E2[m] = v;	
	AdjMat[u][v] = AdjMat[v][u] = 1;
	I[u][v] = m; 
	I[v][u] = m++;
	AdjList[u][sizeAdjList[u]++]=v;
	AdjList[v][sizeAdjList[v]++]=u; 
}
void Matching::AddEdge(int u, int v, double c)
{	
	cost[m] = c;
	E1[m] = u;
	E2[m] = v;
	AdjMat[u][v] = AdjMat[v][u] = 1;
	C[u][v] = C[v][u] = c;
	I[u][v] = m; 
	I[v][u] = m++; 
	AdjList[u][sizeAdjList[u]++]=v;
	AdjList[v][sizeAdjList[v]++]=u;
}

int Matching::IsInMatching(int u, int v)
{
	return matching[I[u][v]];
}

void Matching::DeleteEdges()
{
	m = 0;
	for(int i = 0; i < n; i++)
	{
		memset(AdjMat[i], 0, n*sizeof(int));
		memset(Artificial[i], 0, n*sizeof(int));
	}
	memset(sizeAdjList, 0, n*sizeof(int));
}

//Grows an alternating forest
void Matching::Grow()
{
	int hungarian = false;

	Reset();

	while(!BFSList.empty())
	{
		int w = outer[BFSList.front()];
		BFSList.pop_front();

		for(int k = 0; k < sizeDeep[w]; k++)
		{
			int u = deep[w][k];
			int cont = false;
			for(int i = 0; i < sizeAdjList[u]; i++)
			{
				int v = AdjList[u][i];

				//Check if the edge is blocked
				if(slack[I[u][v]] > 0.00001) continue;

				//u is even and v is odd
				if(type[outer[v]] == ODD) continue;	

				//if v is unlabeled
				if(type[outer[v]] != EVEN)
				{
					//We grow the alternating forest
					int vm = mate[outer[v]];

					forest[outer[v]] = u;
					type[outer[v]] = ODD;
					root[outer[v]] = root[outer[u]];
					forest[outer[vm]] = v;
					type[outer[vm]] = EVEN;
					root[outer[vm]] = root[outer[u]];

					if(!visited[outer[vm]])
					{
						BFSList.push_back(vm);
						visited[outer[vm]] = true;
					}
				}
				//If v is even and u and v are on different trees
				//we found an augmenting path
				else if(root[outer[v]] != root[outer[u]])
				{
					Augment(u,v);
					Reset();

					cont = true;
					break;
				}
				//If u and v are even and on the same tree
				//we found a blossom
				else if(outer[u] != outer[v])
				{
					int b = Blossom(u,v);

					BFSList.push_front(deep[b][0]);
					visited[b] = true;

					cont = true;
					break;
				} 
			}
			if(cont) break;
		}
	}
/*
printf("....\n");
	for(int i = 0; i < m; i++)
	{
		if(slack[i] > 0.00001) continue;
		
		printf("%d %d", E1[i], E2[i]);
		if(mate[E1[i]] == E2[i]) printf("**\n");
		else printf("\n");
	}
	printf("----\n");
	for(int i = n; i < 2*n; i++)
	{
		printf("blossom ");
		if(blocked[i] && outer[i] == i)
		{
			for(int j = 0; j < sizeDeep[i]; j++)
			{
				printf("%d ", deep[i][j]);
			}
		}
		printf("\n");
	}
	for(int i = 0; i < n; i++)
	{
		printf("type v%d ov%d tov%d mov%d\n", i, outer[i], type[outer[i]], mate[outer[i]]);
	}
*/
	//We proceed until the forest is hungarian	
/*	while(!hungarian)
	{
		hungarian = true;

		//Look for an even vertex u. connected to an unlabeled or even node v
		for(int k = 0; k < sizeE; k++)
		{
			int i = E[k];

			//Check if the edge is blocked
			if(slack[i] > 0.00001) continue;

			int u = E1[i];
			int v = E2[i];

			//Let u be the one that is possibly even
			if(type[outer[u]] != EVEN){int aux = u; u = v; v = aux;}

			//if none is even
			if(type[outer[u]] != EVEN) continue;
			//u is even and v is odd
			if(type[outer[v]] == ODD) continue;	

			//if v is unlabeled
			if(type[outer[v]] != EVEN)
			{
				//We grow the alternating forest
				int vm = mate[outer[v]];

				forest[outer[v]] = u;
				type[outer[v]] = ODD;
				root[outer[v]] = root[outer[u]];
				forest[outer[vm]] = v;
				type[outer[vm]] = EVEN;
				root[outer[vm]] = root[outer[u]];

				hungarian = false;
				break;
			}
			//If v is even and u and v are on different trees
			//we found an augmenting path
			else if(root[outer[v]] != root[outer[u]])
			{
				Augment(u,v);
				Reset();
				hungarian = false;
				break;
			}
			//If u and v are even and on the same tree
			//we found a blossom
			else if(outer[u] != outer[v])
			{
				Blossom(u,v);

				hungarian = false;
				break;
			} 
		}
	}
*/
	//Check if the matching is perfect
	perfect = true;
	for(int i = 0; i < n; i++)
	{
		if(mate[outer[i]] == -1)
		{
			perfect = false;
		}
	}
}

void Matching::Heuristic()
{
	for(int i = 0; i < m; i++)
	{
		if(slack[i] > 0.00001) continue;

		int u = E1[i];
		int v = E2[i];

		if(outer[u] != outer[v] && mate[outer[u]] == -1 && mate[outer[v]] == -1)
		{
			mate[outer[u]] = v;
			mate[outer[v]] = u;
		}
	}
	Reset();
}

void Matching::Heuristic2()
{
	Bheap->Clear();
	for(int i = 0; i < n; i++)
		bhnodes[i]->key = 0;	

	for(int i = 0; i < m; i++)
	{
		if(slack[i] > 0.00001) continue;

		int u = E1[i];
		int v = E2[i];

		bhnodes[u]->key = bhnodes[u]->key +1;
		bhnodes[v]->key = bhnodes[v]->key +1;
	}
	for(int i = 0; i < n; i++)
	{
		Bheap->Insert(bhnodes[i]);
	}
	BHNode *next =	Bheap->DeleteMin();
	while(next)
	{
		int u = next->p;
		if(mate[outer[u]] == -1)
		{
			int min = -1;
			for(int i = 0; i < sizeAdjList[u]; i++)
			{
				int v = AdjList[u][i];
				if(slack[I[u][v]] > 0.00001) continue;
				if(outer[u] == outer[v]) continue;
				if(mate[outer[v]] != -1) continue;

				if(min == -1) 
				{
					min = v;
				}
				else if(bhnodes[v]->key < bhnodes[min]->key)
				{
					min = v;	
				}
			}
			if(min != -1)
			{
				mate[outer[u]] = min;
				mate[outer[min]] = u;
			}
		}

		next = Bheap->DeleteMin();
	}

	Reset();
}

//Destroys a blossom recursively
void Matching::DestroyBlossom(int t)
{
	if(t < n) return;
	if(blocked[t] && dual[t] > 0.00001) return;

	for(int i = 0; i < sizeShallow[t]; i++)
	{
		int s = shallow[t][i];
		blossom[s] = s;
		outer[s] = s;
		for(int j = 0; j < sizeDeep[s]; j++)
		{
			outer[deep[s][j]] = s;	
		}
	}
	for(int i = 0; i < sizeShallow[t]; i++)
		DestroyBlossom(shallow[t][i]);

	active[t] = false;
	blocked[t] = false;
	Free[sizeFree++] = t;
	mate[t] = -1;
}

void Matching::Expand(int u)
{
	int v = mate[u];

	int index = 10000000;
	int p, q;
	//Find the regular edge {p,q} of minimum index I[p,q] connecting u and its mate
	for(int i = 0; i < sizeDeep[u]; i++)
	{	
		int di = deep[u][i];
		for(int j = 0; j < sizeDeep[v]; j++)
		{
			int dj = deep[v][j];
			if(AdjMat[di][dj] && !(slack[ I[di][dj] ] > 0.00001) && I[di][dj] < index)
			{
				index = I[di][dj];
			
				p = di;
				q = dj;
			}
		}
	}

	mate[u] = q;
	//If u is a regular vertex, we are done
	if(u < n || blocked[u]) return;

	int t;
	//Find the position t of the new tip of the blossom
	for(int i = 0; i < sizeShallow[u]; i++)
	{
		int found = false;
		int si = shallow[u][i];

		for(int j = 0; j < sizeDeep[ si ]; j++)
		{
			if(deep[ si ][j] == p )
			{
				t = i;
				found = true; break;
			}
		}
		if(found) break;
	}

	//Adjust the mate of the tip
	mate[shallow[u][t]] = mate[u];

	//Now we go through the odd circuit adjusting the new mates
	for(int i = t+1; i % sizeShallow[u] != t; i++)
	{
		if((i - t)%2)
		{
			mate[ shallow[u][ i%sizeShallow[u] ] ] = shallow[u][ (i+1)%sizeShallow[u] ];	
		}
		else
		{
			mate[ shallow[u][ i%sizeShallow[u] ] ] = shallow[u][ (i-1)%sizeShallow[u] ];	
		}	
	}

	//We update the sets blossom, shallow, and outer since this blossom is being deactivated
	for(int i = 0; i < sizeShallow[u]; i++)
	{
		int s = shallow[u][i];
		blossom[s] = s;
		outer[s] = s;
		for(int j = 0; j < sizeDeep[s]; j++)
		{
			outer[deep[s][j]] = s;	
		}
	}
	active[u] = false;
	Free[sizeFree++] = u;

	//We expand the vertices in the blossom
	for(int i = 0; i < sizeShallow[u]; i++)
	{
		Expand(shallow[u][i]);
	}
}


void Matching::Expand2(int u, int p, int q)
{
	mate[u] = q;
	//If u is a regular vertex, we are done
	if(u < n || blocked[u]) return;

	int t;
	//Find the position t of the new tip of the blossom
	for(int i = 0; i < sizeShallow[u]; i++)
	{
		int found = false;
		int si = shallow[u][i];

		for(int j = 0; j < sizeDeep[ si ]; j++)
		{
			if(deep[ si ][j] == p )
			{
				t = i;
				found = true; break;
			}
		}
		if(found) break;
	}

	//Adjust the mate of the tip
	mate[shallow[u][t]] = mate[u];

	//Now we go through the odd circuit adjusting the new mates
	for(int i = t+1; i % sizeShallow[u] != t; i++)
	{
		if((i - t)%2)
		{
			mate[ shallow[u][ i%sizeShallow[u] ] ] = shallow[u][ (i+1)%sizeShallow[u] ];	
		}
		else
		{
			mate[ shallow[u][ i%sizeShallow[u] ] ] = shallow[u][ (i-1)%sizeShallow[u] ];	
		}	
	}

	//We update the sets blossom, shallow, and outer since this blossom is being deactivated
	for(int i = 0; i < sizeShallow[u]; i++)
	{
		int s = shallow[u][i];
		blossom[s] = s;
		outer[s] = s;
		for(int j = 0; j < sizeDeep[s]; j++)
		{
			outer[deep[s][j]] = s;	
		}
	}
	active[u] = false;
	Free[sizeFree++] = u;

	//We expand the vertices in the blossom
	for(int i = 0; i < sizeShallow[u]; i++)
	{
		if(i == t) Expand2(shallow[u][i], p, q);
		else Expand(shallow[u][i]);
	}
}

//Augment the path root[u], ..., u, v, ..., root[v]
void Matching::Augment(int u, int v)
{
	memset(auxVertexArray1, 0, 2*n*sizeof(int));
	int sizePath = 0;

	int p, q;

	//We go from u and v to its respective roots, alternating the matching
	p = outer[u];
	q = outer[v];
	mate[p] = q;
	mate[q] = p;
	auxVertexArray1[sizePath++] = p;
	auxVertexArray1[sizePath++] = q;
	while(forest[p] != -1)
	{
		q = outer[forest[p]];
		p = outer[forest[q]];

		mate[p] = q;
		mate[q] = p;
		auxVertexArray1[sizePath++] = p;
		auxVertexArray1[sizePath++] = q;
	}

	p = outer[v];
	while(forest[p] != -1)
	{
		q = outer[forest[p]];
		p = outer[forest[q]];

		mate[p] = q;
		mate[q] = p;
		auxVertexArray1[sizePath++] = p;
		auxVertexArray1[sizePath++] = q;
	}

	//Some vertices may be blossoms, they need to be expanded
	for(int i = 0; i < sizePath; i++)
	{
		Expand(auxVertexArray1[i]);
	}
}

void Matching::Reset()
{
	for(int i = 0; i < 2*n; i++)
	{
		forest[i] = -1;
		root[i] = i;

		if(i >= n && active[i] && outer[i] == i)
			DestroyBlossom(i);
	}

	memset(visited, 0, 2*n*sizeof(int));
	BFSList.clear();
	for(int i = 0; i < n; i++)
	{
		if(mate[outer[i]] == -1)
		{
			type[outer[i]] = 2;
			if(!visited[outer[i]])
				BFSList.push_back(i);
			visited[outer[i]] = true;
		}
		else type[outer[i]] = 0;
	}
}

//Returns a free index for a new blossom
int Matching::GetFreeIndex()
{
	int i = Free[0];
	Free[0] = Free[--sizeFree];
	return i;
}

//Contracts the blossom t, ..., u, v, ..., t, where t is the first vertex that appears in the paths from u and v to their respective roots
int Matching::Blossom(int u, int v)
{
	int t = GetFreeIndex();

	memset(auxVertexArray1, 0, 2*n*sizeof(int));
	memset(auxVertexArray2, 0, 2*n*sizeof(int));

	//Find the tip of the blossom
	int u_ = u, v_ = v;
	while(true)
	{
		auxVertexArray1[outer[u_]] = true;
		auxVertexArray2[outer[v_]] = true;

		if(auxVertexArray1[outer[v_]])
		{
			tip[t] = outer[v_];
			break;
		}
		if(auxVertexArray2[outer[u_]])
		{
			tip[t] = outer[u_];
			break;
		}
		
		if(forest[outer[u_]] != -1) u_ = forest[outer[u_]];
		if(forest[outer[v_]] != -1) v_ = forest[outer[v_]];
	}

	//Find the odd circuit, update shallow, outer, blossom and deep
	//First we construct the set shallow (the odd circuit)
	u_ = outer[u];
	int sizePathu = 0;
	auxVertexArray1[sizePathu++] = u_;
	while(u_ != tip[t])
	{
		u_ = outer[forest[u_]];
		auxVertexArray1[sizePathu++] = u_;
	}

	sizeShallow[t] = 0;
	sizeDeep[t] = 0;
	for(int i = 0; i < sizePathu; i++)
	{
		u_ = auxVertexArray1[sizePathu-i-1];
		shallow[t][i] = u_;
		sizeShallow[t]++;
	}

	v_ = outer[v];
	int sizePathv = 0;
	while(v_ != tip[t])
	{
		shallow[t][sizeShallow[t]] = v_;
		
		v_ = outer[forest[v_]];
		sizeShallow[t]++;
	}
	//Now we construct deep and update outer
	for(int i = 0; i < sizeShallow[t]; i++)
	{
		u_ = shallow[t][i];
		blossom[u_] = t;
		outer[u_] = t;
		for(int j = 0; j < sizeDeep[u_]; j++)
		{
			deep[t][ sizeDeep[t]++ ] = deep[u_ ][j];
			outer[deep[u_][j]] = t;
		}
	}

	forest[t] = forest[tip[t]];
	type[t] = EVEN;
	root[t] = root[tip[t]];
	active[t] = true;
	outer[t] = t;
	mate[t] = mate[tip[t]];

	return t;
}

void Matching::UpdateDualCosts()
{
	double e1, e2, e3;
	int inite1 = false, inite2 = false, inite3 = false;
	for(int i = 0; i < m; i++)
	{
		int u = E1[i], v = E2[i];
		if( (type[outer[u]] == EVEN && type[outer[v]] == UNLABELED) || (type[outer[v]] == EVEN && type[outer[u]] == UNLABELED) )
		{
			if(!inite1 || (e1 - slack[i] > 0.00001))
			{
				e1 = slack[i];
				inite1 = true;
			}
		}
		else if( (outer[u] != outer[v]) && type[outer[u]] == EVEN && type[outer[v]] == EVEN )
		{
			if(!inite2 || (e2 - slack[i] > 0.00001))
			{
				e2 = slack[i];
				inite2 = true;
			}
		}
	}
	for(int i = n; i < 2*n; i++)
	{
		if(active[i] && i == outer[i] && type[outer[i]] == ODD && (!inite3 || (e3 - dual[i] > 0.00001)))
		{
			e3 = dual[i]; 
			inite3 = true;
		}	
	}
	double e;
	if(inite1) e = e1;
	else if(inite2) e = e2;
	else if(inite3) e = e3;

	if(e - e2/2.0 > 0.00001 && inite2)
		e = e2/2.0;
	if(e - e3 > 0.00001 && inite3)
		e = e3;
	 
	for(int i = 0; i < 2*n; i++)
	{
		if(i != outer[i]) continue;

		if(active[i] && type[outer[i]] == EVEN)	
		{
			dual[i] += e; 
		}
		else if(active[i] && type[outer[i]] == ODD)
		{
			dual[i] -= e; 
		}
	}

	sizeE = 0;
	for(int i = 0; i < m; i++)
	{
		int u = E1[i], v = E2[i];			

		if(outer[u] != outer[v])
		{	
			if(type[outer[u]] == EVEN && type[outer[v]] == EVEN)
			{
				slack[i] -= 2.0*e;
			}
			else if(type[outer[u]] == ODD && type[outer[v]] == ODD)
			{
				slack[i] += 2.0*e;
			}	
			else if( (type[outer[v]] == UNLABELED && type[outer[u]] == EVEN) || (type[outer[u]] == UNLABELED && type[outer[v]] == EVEN) )
			{
				slack[i] -= e;
			}	
			else if( (type[outer[v]] == UNLABELED && type[outer[u]] == ODD) || (type[outer[u]] == UNLABELED && type[outer[v]] == ODD) )
			{
				slack[i] += e;
			}
		}
		
		if(!(slack[i] > 0.00001))
		{ E[sizeE++] = i;
			AdjList[u][sizeAdjList[u]++] = v;
			AdjList[v][sizeAdjList[v]++] = u;
		}
	}
	for(int i = n; i < 2*n; i++)
	{
		if(dual[i] > 0.00001)
		{
			if(blocked[i] == false)
			{
				//the blossom is becoming blocked
				mate[i] = mate[tip[i]];	
			}

			blocked[i] = true;
		}
		else if(active[i] && blocked[i])
		{
			//The blossom is becoming unblocked
			if(mate[i] == -1)
			{
				DestroyBlossom(i);
			}
			else
			{
				blocked[i] = false;
				int q = mate[i];
				int p = mate[outer[q]];
				Expand2(i, p, q);
			}
		}
	}	
}

void Matching::SolveMinimumCostPerfectMatching()
{
	Clear();
	PositiveCosts();
	GrantFeasibility();
	feasible = true;

	//Initialize slacks (reduced costs for the edges)
	for(int i = 0; i < m; i++)
		slack[i] = cost[i];

		memset(sizeAdjList, 0, n*sizeof(int));
	while(true)
	{
		Heuristic2();

		//Grow a new hungarian forest
		Grow();
		//If the matching on the compressed graph is perfect, we are done
		if(perfect) break;
		//Otherwise we update the dual costs	
		memset(sizeAdjList, 0, n*sizeof(int));	
		UpdateDualCosts();
		//Set up the algorithm for a new grow step
		Reset();
	}

	//Find the actual matching
	for(int i = 0; i < 2*n; i++)
	{
		if(active[i] && outer[i] == i)
		{
			int q = mate[i];
			int p = mate[outer[q]];
			Open2(i, p, q);
		}
	}

	obj = 0;
	for(int i = 0; i < m; i++)
	{
		int u = E1[i];
		int v = E2[i];
		if(!Artificial[u][v])
		cost[i] += minEdge;
		if(mate[u] == v)
		{
			if(Artificial[u][v])
			{
				feasible = false;
				//printf("->>>>> %d %d\n", u, v); 
				break;
			}

			matching[i] = 1;
			obj += cost[i];
		}
		else
		{
			matching[i] = 0;
		}
	}

	double dualObj = 0;
	for(int i = 0; i < 2*n; i++)
	{
		if(i < n) dualObj += dual[i];
		else if(blocked[i]) dualObj += dual[i];	
	}

	//PrintNeato();

//	printf("primal obj: %.4f  dual obj: %.4f\n", obj, dualObj);
//	if(fabs(obj - dualObj) > 0.0001) printf("diff-------------------------------------------------------------------------\n");
}

void Matching::PositiveCosts()
{
	minEdge = 0;
	for(int i = 0; i < m ;i++)
	{
		if(minEdge - cost[i] > 0.00001) minEdge = cost[i];
	}
//printf("me: %.4f\n", minEdge);
	for(int i = 0; i < m; i++)
	{
		if(!Artificial[E1[i]][E2[i]])
			cost[i] -= minEdge;
//		printf("%.3f ", cost[i]);
	}
//	printf("\n");
}

void Matching::GrantFeasibility()
{
	maxCost = 0;
	for(int i = 0; i < m; i++)
		maxCost += cost[i];
	for(int i = 0; i < n-1; i++)
	{
		if(!AdjMat[i][i+1])
		{
			AddEdge(i, i+1, maxCost);
			Artificial[i][i+1] = true;
			Artificial[i+1][i] = true;
		}
		else if(Artificial[i][i+1])
		{
			SetCost(i, i+1, maxCost);
		}
	}
	if(!AdjMat[0][n-1])
	{
		AddEdge(0, n-1, maxCost+1);
		Artificial[0][n-1] = true;
		Artificial[n-1][0] = true;
	}
	else if(Artificial[0][n-1])
	{
		SetCost(0, n-1, maxCost);
	}
}

void Matching::SolvePerfectMatching()
{
	Clear();
	Grow();
	PrintNeato();
}

//Sets up the algorithm for a new run
void Matching::Clear()
{
	for(int i = 0; i < 2*n; i++)
	{
		sizeFree = 0;
		for(int j = n; j < 2*n; j++)
		{
			Free[sizeFree++] = j;
		}
		blossom[i] = i;
		outer[i] = i;
		sizeDeep[i] = 0;
		if(i<n)
		{
			deep[i][sizeDeep[i]++] = i;
		}
		sizeShallow[i] = 0;
		if(i < n) active[i] = true;
		else active[i] = false;
	
		type[i] = 0;
		forest[i] = -1;
		root[i] = i;

		blocked[i] = false;
		dual[i] = 0;
		mate[i] = -1;
		tip[i] = i;
	}
	for(int i = 0; i < (n*(n-1))/2; i++)
	{
		slack[i] = 0;
		matching[i] = 0;
	}
	
	sizeE = 0;
}

void Matching::Open(int u)
{
	int v = mate[u];
	int index = 10000000;
	int p, q;
	for(int i = 0; i < sizeDeep[u]; i++)
	{	
		int di = deep[u][i];
		for(int j = 0; j < sizeDeep[v]; j++)
		{
			int dj = deep[v][j];
			if(AdjMat[di][dj] && !(slack[ I[di][dj] ] > 0.00001) && I[di][dj] < index)
			{
				index = I[di][dj];
				p = di;
				q = dj;
			}
		}
	}
	mate[u] = q;
	if(u < n) return;

	int t;
	for(int i = 0; i < sizeShallow[u]; i++)
	{
		int found = false;
		int si = shallow[u][i];
		for(int j = 0; j < sizeDeep[ si ]; j++)
		{
			if(deep[ si ][j] == p )
			{
				t = i;
				found = true; break;
			}
		}
		if(found) break;
	}

	mate[shallow[u][t]] = mate[u];

	for(int i = t+1; i % sizeShallow[u] != t; i++)
	{
		if((i - t)%2)
		{
			mate[ shallow[u][ i%sizeShallow[u] ] ] = shallow[u][ (i+1)%sizeShallow[u] ];	
		}
		else
		{
			mate[ shallow[u][ i%sizeShallow[u] ] ] = shallow[u][ (i-1)%sizeShallow[u] ];	
		}	
	}

	for(int i = 0; i < sizeShallow[u]; i++)
	{
		Open(shallow[u][i]);
	}
}

void Matching::Open2(int u, int p, int q)
{
	mate[u] = q;
	if(u < n) return;

	int t;
	for(int i = 0; i < sizeShallow[u]; i++)
	{
		int found = false;
		int si = shallow[u][i];
		for(int j = 0; j < sizeDeep[ si ]; j++)
		{
			if(deep[ si ][j] == p )
			{
				t = i;
				found = true; break;
			}
		}
		if(found) break;
	}

	mate[shallow[u][t]] = mate[u];

	for(int i = t+1; i % sizeShallow[u] != t; i++)
	{
		if((i - t)%2)
		{
			mate[ shallow[u][ i%sizeShallow[u] ] ] = shallow[u][ (i+1)%sizeShallow[u] ];	
		}
		else
		{
			mate[ shallow[u][ i%sizeShallow[u] ] ] = shallow[u][ (i-1)%sizeShallow[u] ];	
		}	
	}

	for(int i = 0; i < sizeShallow[u]; i++)
	{
		if(i == t) Open2(shallow[u][i], p, q);
		else Open(shallow[u][i]);
	}
}

void Matching::PrintNeato()
{
//	printf("graph G{\n");
//	for(int i = 0; i < n; i++)
//	{
//		printf("n%d [height=0.1, width=0.1, label=\"\"];\n", i);
//	}
	for(int i = 0; i < m; i++)
	{
		if(mate[E1[i]] != E2[i] && !(slack[i] > 0.00001))
			printf("n%d -- n%d [color = lightgrey]; \n", E1[i], E2[i]);
		else if(mate[E1[i]] == E2[i] && mate[E2[i]] == E1[i])
			printf("n%d -- n%d [color = red]; \n", E1[i], E2[i]);

//		if(mate[E1[i]] == E2[i])
//			printf("n%d -- n%d [color = red]; \n", E1[i], E2[i]);
//		else 
//			printf("n%d -- n%d [color = lightgrey]; \n", E1[i], E2[i]);
	}
	printf("\n}\n");
}

void Matching::PrintNeato2()
{
	printf("graph G{\n");
	for(int i = 0; i < n; i++)
	{
		printf("n%d;\n", i);
	}
	for(int i = 0; i < m; i++)
	{
		if(mate[E1[i]] != E2[i])
			printf("n%d -- n%d; \n", E1[i], E2[i]);
		else
			printf("n%d -- n%d [color = red]; \n", E1[i], E2[i]);
	}
	printf("\n}\n");
}

void Matching::SetCost(int u, int v, double c)
{
	C[u][v] = C[v][u] = c;
	cost[I[u][v]] = c;			
}

Matching* Matching::getClone(){
	return NULL;
}


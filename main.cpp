#include "Matching.h"
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <vector>
#include <set>

#include "MST.h"

#define MAX_COST 100
#define MIN_COST 1


Edge* generateRandomCompleteGraph(int numVertices){ //non-direction graph
	assert(numVertices % 2 == 0);
    Edge* edges = (Edge*) malloc (sizeof(Edge) * numVertices * (numVertices-1));
    int count = 0;
    for(int i=0;i<numVertices; ++i){
        for(int j=i+1;j<numVertices;++j){
            edges[count].v1 = i;
            edges[count].v2 = j;
            edges[count].cost = (rand() % (MAX_COST-MIN_COST))+(MIN_COST);
            count++;
        }
    }
    //printf("added %d edges\n", count);
    return edges;
}

void printEdgesInMatching(Edge* edges, int numEdges, Matching* M){
    for(int i=0;i<numEdges;++i) 
	    if(M->IsInMatching(edges[i].v1, edges[i].v2)) printf("{%d, %d} ", edges[i].v1,edges[i].v2);
	printf("\n");
}

void printEdges(Edges edges){
	for(int i=0;i<edges.numEdges;++i)
		printf("{%d, %d, %d} ", edges.edges[i].v1,edges.edges[i].v2, edges.edges[i].cost);
	printf("\n");
}

void printEdgesSimple(Edges edges){
	for(int i=0;i<edges.numEdges;++i)
		printf("{%d, %d} ", edges.edges[i].v1,edges.edges[i].v2);
	printf("\n");
}

void addEdges(Edges edges, Matching* M){
    for(int i=0;i<edges.numEdges;++i){
	    M->AddEdge(edges.edges[i].v1, edges.edges[i].v2, edges.edges[i].cost);
    }
}

int compareEdges(const void * a, const void * b)
{
  if ( (*(Edge*)a).cost <  (*(Edge*)b).cost ) return -1;
  if ( (*(Edge*)a).cost == (*(Edge*)b).cost ) return 0;
  if ( (*(Edge*)a).cost >  (*(Edge*)b).cost ) return 1;
}

void ordenarArestas(Edge* edges, int numEdges){
	qsort (edges, numEdges, sizeof(Edge), compareEdges);
}

Edges getMST(Edges edges, int cost){
	std::set<int> Vn; //Vn = vertices incident to a nonpositive weight edge
	std::vector<Edge> En; //En = edges with both endpoints in Vn
	for(int i=0;i<edges.numEdges;i++){
		if (edges.edges[i].cost <= cost){ //this is equivalent to checking if the new cost is less than or equal to zero
			Vn.insert(edges.edges[i].v1);
			Vn.insert(edges.edges[i].v2);
			En.push_back(edges.edges[i]);
		}
	}
	//prepare structures to call MST
	int uE[En.size()], vE[En.size()];
	double c[En.size()];
	int numVerticesMST = Vn.size();
	int numEdgesMST = En.size();
	/*
	 * Insert the elements into uE and vE after the translation.
	 * This translation is done so we our edges are all between 0 and numVerticesMST.
	 */
	for(int i=0;i<En.size();i++){
		uE[i] = std::distance(Vn.begin(),Vn.find(En[i].v1));
		vE[i] = std::distance(Vn.begin(),Vn.find(En[i].v2));
		c[i] = En[i].cost;
	}
	MST mst(numVerticesMST, numEdgesMST, uE, vE, c);
	mst.SolveKruskal();

	Edges result;
	result.numEdges = numVerticesMST-1;
	Edge* edgesResult = new Edge[result.numEdges];
	result.edges = edgesResult;
	for(int i = 0; i < result.numEdges; i++){
		edgesResult[i] = En[mst.sol[i]];
	}
	return result;
}

Edges getMinCover(Edges edges, int cost){
	/* TODO: I need my paper so I can lookup how I will do this.
	 * I know I will have to provide it with different costs and also add some zero-cost edges but I do not remember all the details.
	 */
	/*
	The steps are:
	Find the edges and vertices for the MinCover in edges. (Vp and Ep)
	For every vertex v find bv, such that bv is the cost of the edge incident to v of smallest cost.
	For every edge e(incident to u and v), calculate the new cost Ze, such that Ze = Ce - bv - bu.
	Create a MST object adding the edges whose cost is greater than 0, and complete the graph with zero-cost-edges
	*/
	//printf("Start min cover algorithm");
	std::set<int> Vp; //Vp = vertices incident to a positive weight edge
	std::vector<Edge> Ep; //Ep = edges with at least one endpoint in Vp
	for(int i=0;i<edges.numEdges;i++){
		if (edges.edges[i].cost > cost){ //this is equivalent to checking if the new cost is greater than zero
			Vp.insert(edges.edges[i].v1);
			Vp.insert(edges.edges[i].v2);
		}
	}
	std::set<int> V_used(Vp);
	for(int i=0;i<edges.numEdges;i++){
		if (Vp.count(edges.edges[i].v1) == 1 || Vp.count(edges.edges[i].v2) == 1){ //if one of the endpoints is in Vp
			Ep.push_back(edges.edges[i]);
			V_used.insert(edges.edges[i].v1); //add edge if to V_used if not in Vp
			V_used.insert(edges.edges[i].v2);
		}
	}
	//printf("Get best costs");
	int b_cost[V_used.size()]; //best costs associated to each vertex
	int aux_cost, currVertex;
	std::vector<int> Vp_vector( V_used.begin(), V_used.end() );
	for(int i=0;i<V_used.size();i++){
		aux_cost = MAX_COST;
		currVertex = Vp_vector[i];
		for(int j=0;j<Ep.size();j++){
			if (Ep[j].v1 == currVertex || Ep[j].v2 == currVertex) //if one of the endpoints is the current vertex
				if (Ep[j].cost < aux_cost) //if the cost is less than the best one so far
					aux_cost = Ep[j].cost;
		}
		//printf("vertex %d b_cost %d\n", Vp_vector[i], aux_cost);
		b_cost[std::distance(V_used.begin(),V_used.find(currVertex))] = aux_cost; //change according to new index position
	}
	int new_cost[V_used.size() + (V_used.size() % 2)][V_used.size() + (V_used.size() % 2)];
	//initialize all new_costs with zero, therefore completing the graph also
	for(int i=0;i<V_used.size() + (V_used.size() % 2);i++){
		for(int j=i+1;j<V_used.size() + (V_used.size() % 2);j++){ //matching expects non-directional edges
			new_cost[i][j] = 0;
		}
	}
	//change the costs of the edges if the new cost is less than zero
	int v1_aux, v2_aux;
	for(int i=0;i<Ep.size();++i){
		//if (Vp.find(Ep[i].v1) == Vp.end() || Vp.find(Ep[i].v2) == Vp.end()) continue; //if one of the vertices are not in Vp, do not change cost
		v1_aux = std::distance(V_used.begin(),V_used.find(Ep[i].v1)); //index change according to the new array
		v2_aux = std::distance(V_used.begin(),V_used.find(Ep[i].v2));
		if (v1_aux > v2_aux) std::swap(v1_aux, v2_aux); //swap the numbers if necessary, since it is nondirectional
		aux_cost = Ep[i].cost - b_cost[v1_aux] - b_cost[v2_aux];
		if (aux_cost < 0)
			new_cost[v1_aux][v2_aux] = aux_cost;
	}
	//printf("Vp:\n");
//	for(int i=0;i<V_used.size();i++){
//		if (Vp.count(Vp_vector[i]))
//			printf("%d ", Vp_vector[i]);
//	}
	//printf("\nVp size: %d", Vp.size());
	//assert(Vp.size() % 2 == 0); //make sure the number of vertices is even. it is treated by adding another vertex at the end
	Matching *M = new Matching((V_used.size()) + ((V_used.size()) % 2)); //possibly add another vertex if not even and include the number of edges not in Vp
	//add all edges and the new costs
	for(int i=0;i<V_used.size() + (V_used.size() % 2);i++){
		for(int j=i+1;j<V_used.size() + (V_used.size() % 2);j++){  //matching expects non-directional edges
			//printf("Added {%d, %d, %d}\n ", i,j, new_cost[i][j]);
			M->AddEdge(i, j, new_cost[i][j]); //add the edges that have both endpoints in Vp
		}
    }

	M->SolveMinimumCostPerfectMatching();
	Edges minCover;
	std::set<int> Vp_left(Vp);
	minCover.edges = (Edge*) malloc(sizeof(Edge)*Vp.size()-1);
	int edgesAdded = 0;
	//printf("\nEdges in matching for min cover\n");
    for(int i=0;i<Ep.size();++i){
    	if(Vp.count(Ep[i].v1) && Vp.count(Ep[i].v2)){ //if both vertices are in Vp. otherwise, I do not care about this edge
			v1_aux = std::distance(V_used.begin(),V_used.find(Ep[i].v1)); //index change according to the new array
			v2_aux = std::distance(V_used.begin(),V_used.find(Ep[i].v2));
			if (v1_aux > v2_aux) std::swap(v1_aux, v2_aux); //swap the numbers if necessary, since it is nondirectional
			if (new_cost[v1_aux][v2_aux] < 0){
				if(M->IsInMatching(v1_aux, v2_aux)){
					Vp_left.erase(Ep[i].v1);
					Vp_left.erase(Ep[i].v2);
					minCover.edges[edgesAdded++] = Ep[i];
					//printf("{%d, %d, %d, %d}\n ", Ep[i].v1,Ep[i].v2, Ep[i].cost, new_cost[v1_aux][v2_aux]);
				}
			}
    	}
    }
	//printf("\n");
//	for(int i=0;i<edgesAdded;++i){
//		printf("{%d, %d, %d} ", minCover.edges[i].v1,minCover.edges[i].v2, minCover.edges[i].cost);
//	}
	//printf("\n");
	//add minimum edges to complete the graph
	//printf("Vertices left: ");
	std::vector<int> Vp_left_vector( Vp_left.begin(), Vp_left.end() );
//	for(int i=0;i<Vp_left_vector.size();i++){
//		printf("%d ", Vp_left_vector[i]);
//	}
	//printf("\n");
	while(Vp_left.size() > 0 ){
		for(int i=0;i<Ep.size();++i){
			//printf("v1 %d v2 %d begin %d cost %d b_cost %d\n", Ep[i].v1, Ep[i].v2, *Vp_left.begin(), Ep[i].cost, b_cost[std::distance(Vp_left.begin(),Vp_left.find(*Vp_left.begin()))]);
			if (Ep[i].v1 == *Vp_left.begin() || Ep[i].v2 == *Vp_left.begin()){ //if edge incident to the first vertex
				if (b_cost[std::distance(V_used.begin(),V_used.find(*Vp_left.begin()))] == Ep[i].cost){ //if it is the best edge
					//printf("found\n");
					Vp_left.erase(*Vp_left.begin());
					minCover.edges[edgesAdded++] = Ep[i];
					break;
				}
			}
		}
		//break;
	}
	for(int i=0;i<edgesAdded;++i){
		printf("{%d, %d, %d} ", minCover.edges[i].v1,minCover.edges[i].v2, minCover.edges[i].cost);
	}
	printf("\n");
}

void printEdges(Edge* edges, int numEdges){
	for(int i=0;i<numEdges;++i)
		printf("%d %d %d\n", edges[i].v1, edges[i].v2, edges[i].cost);
}

int main(int argc, char* argv[])
{
    srand (time(NULL));

	//Número de vértices, TEM que ser par.
	int n = 20;

	//Instancia-se o objeto passando-se o número de vértices desejado.
	Matching *M = new Matching(n);

    Edge* edges = generateRandomCompleteGraph(n);

    int numEdges = (n)*(n-1)/2; //somatorio

    Edges allEdges;
    allEdges.edges = edges;
    allEdges.numEdges = numEdges;

    addEdges(allEdges, M);

    M->SolveMinimumCostPerfectMatching();

    /*
     * P-forest algorithm:
     *
     * Order edges.
     *
     * for v in edges:
     * 	find k forest cover on reduced costs: (I dont really need to reduce the costs, just calculate the edges and vertices accordingly)
     * 		Vn = vertices incident to a nonpositive weight edge
     * 		En = edges with both endpoints in Vn
     * 		G1 = MST in Vn over En
     * 		Vp = V - Vn
     * 		Ep = edges with at least one endpoint in Vp
     * 		G2 = Minimum cover over Vp in Ep
     * 		G = G1 + G2
     * 		if |V| - |G| == k:
     * 			return solution
     */

//    ordenarArestas(edges, numEdges);
//    printf("All Edges:\n");
//    printEdges(allEdges);
//    printf("Matching:\n");
//    printEdgesInMatching(edges, numEdges, M);
//    printf("MST:\n");
//    Edges mst = getMST(allEdges, 100);
//    printEdgesSimple(mst);

    getMinCover(allEdges, 95);

    //printEdges(edges, numEdges);

    //variables used for the MST class
    /*
	The structure used in the MST class and the Matching class are completly different.
	MST code assumes I will have vertices numbered from 0 to m, such that m is the number of vertices.
	That means that I might have to do some translation back and forth to make it into this notation.
	I will have to figure out a way to extract the solution out of it.
    */

    /*
    int* uE, vE, c;
    int numVerticesMST;
    int numEdgesMST;

    for(int i=0;i<numEdges;++i){
    	MST* mst = getMST(allEdges, allEdges.edges[i].cost);
    	Matching* minimumCover = getMinCover(allEdges, allEdges.edges[i].cost);
    }
    */




	delete M;

	return 0;
}

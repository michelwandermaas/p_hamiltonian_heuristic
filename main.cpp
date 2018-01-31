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
    printf("added %d edges\n", count);
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
	std::set<int> Vp; //Vp = vertices incident to a positive weight edge
	std::vector<Edge> Ep; //En = edges with at least one endpoint in Vp
	for(int i=0;i<edges.numEdges;i++){
		if (edges.edges[i].cost > cost){ //this is equivalent to checking if the new cost is greater than zero
			Vp.insert(edges.edges[i].v1);
			Vp.insert(edges.edges[i].v2);
		}
	}

	for(int i=0;i<edges.numEdges;i++){
		if (Vp.count(edges.edges[i].v1) == 1 || Vp.count(edges.edges[i].v2) == 1) //if one of the endpoints is in Vp
			Ep.push_back(edges.edges[i]);
	}
	int b_cost[Vp.size()]; //best costs associated to each vertex
	int aux_cost, currVertex;
	std::vector<int> Vp_vector( Vp.begin(), Vp.end() );
	for(int i=0;i<Vp.size();i++){
		aux_cost = MAX_COST;
		currVertex = Vp_vector[i];
		for(int j=0;j<Ep.size();j++){
			if (Vp.count(Ep[j].v1) == currVertex || Vp.count(Ep[j].v2) == currVertex) //if one of the endpoints is the current vertex
				if (Ep[j].cost < aux_cost) //if the cost is less than the best one so far
					aux_cost = Ep[j].cost;
		}
		b_cost[i] = aux_cost;
	}
	int new_cost[Ep.size()];
	for(int i=0;i<Ep.size();i++){
		new_cost[i] = Ep[i].cost - b_cost[Ep[i].v1] - b_cost[Ep[i].v2];	
	}
	assert(Vp.size() % 2 == 0); //make sure the number of vertices is even
	Matching *M = new Matching(Vp.size());
	//add positive edges
	for(int i=0;i<Ep.size();++i){
		if (new_cost[i] > 0)
			M->AddEdge(Ep[i].v1, Ep[i].v2, new_cost[i]);
    	}
	//add edges to complete the graph
}

void printEdges(Edge* edges, int numEdges){
	for(int i=0;i<numEdges;++i)
		printf("%d %d %d\n", edges[i].v1, edges[i].v2, edges[i].cost);
}

int main(int argc, char* argv[])
{
    srand (time(NULL));

	//Número de vértices, TEM que ser par.
	int n = 10;

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

    ordenarArestas(edges, numEdges);
    printf("All Edges:\n");
    printEdges(allEdges);
    printf("Matching:\n");
    printEdgesInMatching(edges, numEdges, M);
    printf("MST:\n");
    Edges mst = getMST(allEdges, 100);
    printEdgesSimple(mst);

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

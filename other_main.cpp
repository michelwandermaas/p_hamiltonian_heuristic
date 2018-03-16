#include "Matching.h"
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <vector>
#include <set>
#include <cstdio>
#include <limits.h>

#include "MST.h"
#include "Graph.h"

//#define MAX_COST 100
double MAX_COST = LONG_MAX;
double MIN_COST = 0;
//#define MIN_COST 1


int numVertices;

Edge* generateRandomCompleteGraph(int numVertices){ //non-direction graph
    //assert(numVertices % 2 == 0);
    Edge* edges = (Edge*) malloc (sizeof(Edge) * numVertices * (numVertices-1));
    int count = 0;
    for(int i=0;i<numVertices; ++i){
        for(int j=i+1;j<numVertices;++j){
            edges[count].v1 = i;
            edges[count].v2 = j;
            edges[count].cost = (rand() % int(MAX_COST-MIN_COST))+(MIN_COST);
            count++;
        }
    }
    //printf("added %d edges\n", count);
    return edges;
}

Edge* generateCompletePCost(int numVertices, int p_cost){ //non-direction graph
    //assert(numVertices % 2 == 0);
    Edge* edges = (Edge*) malloc (sizeof(Edge) * numVertices * (numVertices-1));
    int count = 0;
    for(int i=0;i<numVertices; ++i){
        for(int j=i+1;j<numVertices;++j){
            edges[count].v1 = i;
            edges[count].v2 = j;
            edges[count].cost = p_cost;
            count++;
        }
    }
    //printf("added %d edges\n", count);
    return edges;
}


void printEdgesInMatching(Edge* edges, int numEdges, Matching* M){
    	int costSum = 0;
    	for(int i=0;i<numEdges;++i){
		if(M->IsInMatching(edges[i].v1, edges[i].v2)){
			printf("{%d, %d} ", edges[i].v1,edges[i].v2);
			costSum += edges[i].cost;
		}
	}
	printf("\n");
	printf("total cost: %d\n", costSum);
}

void printEdges(Edges edges){
    	int costSum = 0;
	for(int i=0;i<edges.numEdges;++i){
		printf("{%d, %d, %d} ", edges.edges[i].v1,edges.edges[i].v2, edges.edges[i].cost);
		costSum += edges.edges[i].cost;
	}
	printf("\n");
	printf("total cost: %d\n", costSum);
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

Edges getMST(Edges edges, double cost){
	/*
	 *      Vn = vertices incident to a nonpositive weight edge
     * 		En = edges with both endpoints in Vn
     * 		G1 = MST in Vn over En
	 */
	Edges result;
	std::set<int> Vn; //Vn = vertices incident to a nonpositive weight edge
	std::vector<Edge> En; //En = edges with both endpoints in Vn
	for(int i=0;i<edges.numEdges;i++){
		if (edges.edges[i].cost <= cost){ //this is equivalent to checking if the new cost is less than or equal to zero
			Vn.insert(edges.edges[i].v1);
			Vn.insert(edges.edges[i].v2);
		}
	}

	for(int i=0;i<edges.numEdges;i++){
		if (Vn.count(edges.edges[i].v1) == 1 && Vn.count(edges.edges[i].v2) == 1){ //both endpoints in Vn
			En.push_back(edges.edges[i]);
		}
	}

	if (Vn.size() == 0 || En.size() == 0){ //it might never happen
		printf("Nothing for MST.\n");
		result.edges = NULL;
		result.numEdges = 0;
		return result;
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
		c[i] = En[i].cost - cost;
	}
	MST mst(numVerticesMST, numEdgesMST, uE, vE, c);
	mst.SolveKruskal();

	Edge* edgesResult = new Edge[numVerticesMST-1];
	result.edges = edgesResult;
	int numEdges = 0;

	std::vector<Edge> zeroCostEdges; 

	std::set<int> verticesLeftToCover(Vn);

	for(int i = 0; i < numVerticesMST-1; i++){
		if (c[mst.sol[i]] <= 0){ //only add nonpositive edges
			edgesResult[numEdges++] = En[mst.sol[i]];
			verticesLeftToCover.erase(En[mst.sol[i]].v1);
			verticesLeftToCover.erase(En[mst.sol[i]].v2);
		}else if (c[mst.sol[i]] == cost){
			zeroCostEdges.push_back(En[mst.sol[i]]);
		}
	}

	//if (zeroCostEdges.size() > 0){
		//printf("opa, tem custo zero aqui. %d vertices to cover", verticesLeftToCover.size());		
	//}	

	result.numEdges = numEdges;
	return result;
}

Edges getMinCoverCalculation(Edges edges, double cost, std::set<int>& Vp, std::set<int>& V_used, std::vector<Edge>& Ep){

	Edges minCover; //edge for the the return
	//printf("Get best costs");
	double b_cost[V_used.size()]; //best costs associated to each vertex
	double aux_cost;
	int currVertex;
	std::vector<int> Vp_vector( V_used.begin(), V_used.end() );
	//Edge aux;
	for(int i=0;i<V_used.size();i++){
		aux_cost = MAX_COST;
		currVertex = Vp_vector[i];
		for(int j=0;j<Ep.size();j++){
			if (Ep[j].v1 == currVertex || Ep[j].v2 == currVertex) //if one of the endpoints is the current vertex
				if (Ep[j].cost < aux_cost){ //if the cost is less than the best one so far
					aux_cost = Ep[j].cost;
					//aux = Ep[j];
				}
		}
		//printf("vertex %d b_cost %d edge:%d,%d\n", Vp_vector[i], aux_cost-cost, aux.v1, aux.v2);
		b_cost[std::distance(V_used.begin(),V_used.find(currVertex))] = aux_cost - cost; //change according to new index position, subtracting the cost
	}
	double new_cost[V_used.size() + (V_used.size() % 2)][V_used.size() + (V_used.size() % 2)];
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
		aux_cost = (Ep[i].cost - cost) - b_cost[v1_aux] - b_cost[v2_aux];
		if (aux_cost < 0)
			new_cost[v1_aux][v2_aux] = aux_cost;
	}

//	printf("Vp:\n");
//	for(int i=0;i<V_used.size();i++){
//		if (Vp.count(Vp_vector[i]))
//			printf("%d ", Vp_vector[i]);
//	}
//	printf("\n");

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

	std::set<int> Vp_left(Vp);
	minCover.edges = (Edge*) malloc(sizeof(Edge)*V_used.size()-1);

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

	printf("Edges by matching: \n");
	for(int i=0;i<edgesAdded;++i){
		printf("{%d, %d, %d} ", minCover.edges[i].v1,minCover.edges[i].v2, minCover.edges[i].cost);
	}
	printf("\n");
	//add minimum edges to complete the graph
	//printf("Vertices left: ");
	std::vector<int> Vp_left_vector( Vp_left.begin(), Vp_left.end() );
//	for(int i=0;i<Vp_left_vector.size();i++){
//		printf("%d ", Vp_left_vector[i]);
//	}
//	printf("\n");
	while(Vp_left.size() > 0 ){
		//printf("vertex left: %d\n", *Vp_left.begin()); 	
		for(int i=0;i<Ep.size();++i){
			//printf("v1 %d v2 %d begin %d cost %d b_cost %d\n", Ep[i].v1, Ep[i].v2, *Vp_left.begin(), Ep[i].cost, b_cost[std::distance(Vp_left.begin(),Vp_left.find(*Vp_left.begin()))]);
			if (Ep[i].v1 == *Vp_left.begin() || Ep[i].v2 == *Vp_left.begin()){ //if edge incident to the first vertex
				//printf("v1 %d v2 %d begin %d cost %d b_cost %d\n", Ep[i].v1, Ep[i].v2, *Vp_left.begin(), Ep[i].cost-cost, b_cost[std::distance(V_used.begin(),V_used.find(*Vp_left.begin()))]);
				if (b_cost[std::distance(V_used.begin(),V_used.find(*Vp_left.begin()))] == Ep[i].cost-cost){ //if it is the best edge
					//printf("found\n");
					Vp_left.erase(Ep[i].v1);
					Vp_left.erase(Ep[i].v2);
					minCover.edges[edgesAdded++] = Ep[i];
					break;
				}
			}
		}
		//break;
	}
	//printf("%d %d %d\n", edgesAdded, minCover.numEdges, Vp.size());
	//assert(edgesAdded == minCover.numEdges);
	minCover.numEdges = edgesAdded;

//	for(int i=0;i<edgesAdded;++i){
//		printf("{%d, %d, %d} ", minCover.edges[i].v1,minCover.edges[i].v2, minCover.edges[i].cost);
//	}
//	printf("\n");

	return minCover;

}


Edges getMinCover(Edges edges, double cost){
	/*
	The steps are:
	Find the edges and vertices for the MinCover in edges. (Vp and Ep)
	For every vertex v find bv, such that bv is the cost of the edge incident to v of smallest cost.
	For every edge e(incident to u and v), calculate the new cost Ze, such that Ze = Ce - bv - bu.
	Create a MST object adding the edges whose cost is greater than 0, and complete the graph with zero-cost-edges
	*/
	//printf("Start min cover algorithm");

	
	std::set<int> Vn; //Vn = vertices incident to a nonpositive weight edge
	for(int i=0;i<edges.numEdges;i++){
		if (edges.edges[i].cost <= cost){ //this is equivalent to checking if the new cost is less than or equal to zero
			Vn.insert(edges.edges[i].v1);
			Vn.insert(edges.edges[i].v2);
		}
	}
	std::set<int> Vp; //Vp = V - Vn
	for(int i=0;i<edges.numEdges;i++){
		if (Vn.count(edges.edges[i].v1) == 0) //edge not in Vn
			Vp.insert(edges.edges[i].v1);
		if (Vn.count(edges.edges[i].v2) == 0) //edge not in Vn
			Vp.insert(edges.edges[i].v2);
	}

	std::vector<Edge> Ep; //Ep = edges with at least one endpoint in Vp
	std::set<int> V_used(Vp);
	for(int i=0;i<edges.numEdges;i++){
		if (Vp.count(edges.edges[i].v1) == 1 || Vp.count(edges.edges[i].v2) == 1){ //if one of the endpoints is in Vp
			Ep.push_back(edges.edges[i]);
			V_used.insert(edges.edges[i].v1); //add edge to V_used if not in Vp
			V_used.insert(edges.edges[i].v2);
		}
	}

	if (Vp.size() == 0 || Ep.size() == 0){
		//printf("Nothing for Min. Cover.\n");
		Edges minCover; //edge for the the return
		minCover.edges = NULL;
		minCover.numEdges = 0;
		return minCover;
	}

	return getMinCoverCalculation(edges, cost, Vp, V_used, Ep);

}

Edges doubleEdges(Edges edges){ //return all the edges doubled, as if it were directional
	Edges ret;
	ret.edges = (Edge*) malloc(sizeof(Edge)*edges.numEdges*2);
	ret.numEdges = edges.numEdges*2;
	
	for(int i=0;i<edges.numEdges;++i){
		ret.edges[2*i] = edges.edges[i];	
		ret.edges[(2*i)+1].v1 = edges.edges[i].v2;	
		ret.edges[(2*i)+1].v2 = edges.edges[i].v1;	
		ret.edges[(2*i)+1].cost = edges.edges[i].cost;	
	}
	return ret;
}

std::vector<int>* removeEqual(std::vector<int>& vector){
	std::vector<int>* returnVector = new std::vector<int>();
	bool found;
	for(int i=0;i<vector.size();++i){
		found = false;
		for(int j=0;j<returnVector->size();++j){
			if (*(returnVector->begin()+j) == vector[i]){
				found = true;
				break;
			}
		}
		if (!found){
			returnVector->push_back(vector[i]);
		}
	}
	return returnVector;
}

std::vector<Edge>* removeEqualEdge(std::vector<Edge>& vector){
	std::vector<Edge>* returnVector = new std::vector<Edge>();
	bool found;
	for(int i=0;i<vector.size();++i){
		found = false;
		for(int j=0;j<returnVector->size();++j){
			if (*(returnVector->begin()+j) == vector[i]){
				found = true;
				break;
			}
		}
		if (!found){
			returnVector->push_back(vector[i]);
		}
	}
	return returnVector;
}



Edges eulerianCircuit(Edges edges, Edges allEdges){
	std::set<int> unvisitedVertices;	
	std::vector<Edge> unvisitedEdges;
	for(int i=0;i<edges.numEdges;i++){
		unvisitedEdges.push_back(edges.edges[i]);	
		unvisitedVertices.insert(edges.edges[i].v1);
		unvisitedVertices.insert(edges.edges[i].v2);
	}
	int initialVertex = *unvisitedVertices.begin();
	int currVertex = initialVertex;
	unvisitedVertices.erase(currVertex);
	std::vector<int>* path;
	path = new std::vector<int>;
	bool found;
	std::vector<Edge>* edgesUsed = new std::vector<Edge>();
	std::vector<std::vector<int>*> allPaths;
	allPaths.push_back(path);
	while (unvisitedVertices.size() > 0){	
		found = false;
		path->push_back(currVertex);
		for(int i=0;i<unvisitedEdges.size();i++){
			if (unvisitedEdges[i].v1 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v2) == 1){ //if begin at currVertex and ends in a unvisitedVertex
				currVertex = unvisitedEdges[i].v2;
				edgesUsed->push_back(unvisitedEdges[i]);
				unvisitedEdges.erase(unvisitedEdges.begin()+i);
				unvisitedVertices.erase(currVertex);
				found = true;
				break;
			}
		}
		if (!found){ //if no unvisited edge is available, go to some already visited
			for(int i=0;i<unvisitedEdges.size();i++){
				if (unvisitedEdges[i].v1 == currVertex){// && unvisitedVertices.count(edges.edges[i].v2) == 1){ //if begin at currVertex and ends in a unvisitedVertex
					currVertex = unvisitedEdges[i].v2;
					edgesUsed->push_back(unvisitedEdges[i]);
					unvisitedEdges.erase(unvisitedEdges.begin()+i);
					found = true;
					break;
				}
			}
		}
		if (!found){ //that means that that cycle is done, move to next one
			printf("new path\n");
			path->push_back(currVertex);
			std::vector<int>* aux = removeEqual(*path);
			aux->push_back(initialVertex);
			allPaths.push_back(aux);
			path = new std::vector<int>;
			initialVertex = *unvisitedVertices.begin();
			currVertex = initialVertex;
			unvisitedVertices.erase(currVertex);
		}
	}
	path->push_back(currVertex);
	std::vector<int>* aux = removeEqual(*path);
	aux->push_back(initialVertex);
	allPaths.push_back(aux);

	//shortcut
	std::vector<Edge>* edgesInPath = new std::vector<Edge>();
	std::vector<int>* shortcutPath;	
	for(int k=0;k<allPaths.size();++k){
		shortcutPath = allPaths[k];	
		for(int i=0;i<shortcutPath->size()-1;++i){
			for(int j=0;j<allEdges.numEdges;++j){
				if (allEdges.edges[j].v1 == *(shortcutPath->begin()+i) && allEdges.edges[j].v2 == *(shortcutPath->begin()+i+1)){
					edgesInPath->push_back(allEdges.edges[j]);
					break;
				}
			}
		}
	}

	edgesInPath = removeEqualEdge(*edgesInPath);
	Edges circuit;
	circuit.edges = &edgesInPath->front();
	circuit.numEdges = edgesInPath->size();
	return circuit;
}


void printEdges(Edge* edges, int numEdges){
	for(int i=0;i<numEdges;++i)
		printf("%d %d %d\n", edges[i].v1, edges[i].v2, edges[i].cost);
}

Edges readEdges(char* filename){
	Graph graph(filename);
	Edge* edges = (Edge*) malloc (sizeof(Edge) * graph.n * (graph.n-1));
	Edges edgesReturn;
	edgesReturn.edges = edges;
	edgesReturn.numEdges = 	(graph.n/2) * (graph.n-1);
	int numAdded = 0;
	for(int i=0;i<graph.n;++i){
		for(int j=i+1;j<graph.n;++j){
			Edge edge;
			edge.v1 = i;
			edge.v2 = j;
			edge.cost = graph.M[i][j];
			edgesReturn.edges[numAdded++] = edge;
		}
	}
	//printf("%d %d\n", numAdded, edgesReturn.numEdges);
	assert(numAdded == edgesReturn.numEdges);
	numVertices = graph.n;
	return edgesReturn;
}

Edges two_factor(Edges edges){

}

int countCycles(Edges edges){

}

int countCyclesMinVertices(Edges edges, int minVertices){
//for every cycle with more than (minVertices*2), result += (floor(numCycles / minVertices) - 1)

}

int main(int argc, char* argv[]){

	int minArg = 3, maxArg = 3;
	int numTrees;

	if(argc < 2 || argc > maxArg){
		printf("Use --help (or -h) to see possible arguments.\n");
		return 0;
	}

	std::string filename;

	for(int i=1;i<argc;i++){
		std::string x,y;
		x.assign(argv[i]);
		if (!x.compare("--help") || !x.compare("-h")){
			printf("Usage: ./program FILE_NAME NUM_TREES\n");
			return 0;
		}else{
			if (argc < minArg){
				printf("Use --help (or -h) to see possible arguments.\n");
				return 0;
			}
				
			filename.assign(argv[1]);	
			numTrees = strtol(argv[2], NULL, 10);	
			break;
		}
	}

	srand (time(NULL));

	Edges allEdges = readEdges((char*)filename.c_str());

	int costsSum = 0;

	for(int i=0;i<allEdges.numEdges;++i){
		costsSum += allEdges.edges[i].cost;	
	}

    /*
    The algorithm consists of forming an initial 2-factor, which is a set of cycles
    such that each vertex is in exactly one cycle. It can be found in polynomial time using
    matching techniques.
    After finding the 2-factor (call it F) we will then have q cycles, and we will have 3 cases:
    Case 1: q == p. That will be our solution. We guarantee this is the best solution for the problem.
    Case 2: q > p. It means that we need to reduce the number of cycles by adding edges.
        We find a T with is the MST of G. We need to find (q-p) edges in T which we will double and then
        add to F, such that each of these edges connects 2 cycles in F. We guarantee a 3-approximation this way.
    Case 3: q < p. It means that we need to increase the number of cycles by removing edges.
        As long as there are at least (p-q) cycles with 6 or more vertices, we can do this by removing two edges (a,b)
        and (c,d) and then adding edges (a,c) and (b,d), for example. Therefore there's a limit to how many trees we
        can guarantee in this case. We guarantee a 2-approximation this way.

    2-factor algorithm:
        1.For every edge (a,b) in graph, add two vertices ab1 and ab2, remove edge (a,b) and add edges (a,ab1),
        edge (ab1,ab2) and edge (ab2,b), such that c(a,ab1)=c(ab2,b)=c(a,b)/2 and c(ab1,ab2)=0.
        2.Solve Minimum Perfect Matching for the new Graph
        3.If (a, ab1) or (ab2,b) in matching, then edge (a,b) is in the 2-factor.
    */

	if (numTrees < 0 || numTrees > numVertices/2){
		std::cout << "Número inválido de árvores. Número deve estar entre 1-"<< numVertices/2 << std::endl;
		return 0;
	}

    //start the algorithm here

    Edges factor2 = two_factor(allEdges); 
    int numCycles2Factor = countCycles(factor2);

    if (numCycles2Factor == numTrees){
        std::cout << "Solução ótima encontrada." << std::endl;
        printEdges(factor2);
        return 0;
    }else if (numCycles2Factor < numTrees){
    /*
    Case 2: q > p. It means that we need to reduce the number of cycles by adding edges.
        We find a T with is the MST of G. We need to find (q-p) edges in T which we will double and then
        add to F, such that each of these edges connects 2 cycles in F. We guarantee a 3-approximation this way.
    */
        int numCyclesToAdd = countCyclesMinVertices(factor2, 3);
        if (numCyclesToAdd < (numTrees - numCycles2Factor)){
            std::cout << "Não é possível encontrar solução para esse número de árvores." << std::endl;
            return 0;
        }
    }else{
    /*
    Case 2: q > p. It means that we need to reduce the number of cycles by adding edges.
        We find a T with is the MST of G. We need to find (q-p) edges in T which we will double and then
        add to F, such that each of these edges connects 2 cycles in F. We guarantee a 3-approximation this way.
    */
    }
	return 0;
}

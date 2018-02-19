#include "Matching.h"
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <vector>
#include <set>
#include <cstdio>

#include "MST.h"
#include "Graph.h"

#define MAX_COST 100
#define MIN_COST 1


Edge* generateRandomCompleteGraph(int numVertices){ //non-direction graph
    //assert(numVertices % 2 == 0);
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
	for(int i = 0; i < numVerticesMST-1; i++){
		if (c[mst.sol[i]] <= cost) //only add nonpositive edges
			edgesResult[numEdges++] = En[mst.sol[i]];
	}
	result.numEdges = numEdges;
	return result;
}

Edges getMinCover(Edges edges, int cost){
	/*
	The steps are:
	Find the edges and vertices for the MinCover in edges. (Vp and Ep)
	For every vertex v find bv, such that bv is the cost of the edge incident to v of smallest cost.
	For every edge e(incident to u and v), calculate the new cost Ze, such that Ze = Ce - bv - bu.
	Create a MST object adding the edges whose cost is greater than 0, and complete the graph with zero-cost-edges
	*/
	//printf("Start min cover algorithm");

	Edges minCover; //edge for the the return

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
		minCover.edges = NULL;
		minCover.numEdges = 0;
		return minCover;
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
		b_cost[std::distance(V_used.begin(),V_used.find(currVertex))] = aux_cost - cost; //change according to new index position, subtracting the cost
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
		for(int i=0;i<Ep.size();++i){
			//printf("v1 %d v2 %d begin %d cost %d b_cost %d\n", Ep[i].v1, Ep[i].v2, *Vp_left.begin(), Ep[i].cost, b_cost[std::distance(Vp_left.begin(),Vp_left.find(*Vp_left.begin()))]);
			if (Ep[i].v1 == *Vp_left.begin() || Ep[i].v2 == *Vp_left.begin()){ //if edge incident to the first vertex
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
		//printf("%d\n", currVertex);
		//printf("unvisited: ");
		//std::set<int>::iterator x;
		//for(x=unvisitedVertices.begin();x!=unvisitedVertices.end();x++)
		//	printf("%d ", *x);
		//printf("\n");
		found = false;
		path->push_back(currVertex);
		for(int i=0;i<unvisitedEdges.size();i++){
			//printf("%d %d\n", unvisitedEdges[i].v1, unvisitedEdges[i].v2);  
			if (unvisitedEdges[i].v1 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v2) == 1){ //if begin at currVertex and ends in a unvisitedVertex
				currVertex = unvisitedEdges[i].v2;
				//printf("found unvisited %d\n", currVertex);
				edgesUsed->push_back(unvisitedEdges[i]);
				unvisitedEdges.erase(unvisitedEdges.begin()+i);
				unvisitedVertices.erase(currVertex);
				found = true;
				break;
			}
		}
		if (!found){ //if no unvisited edge is available, go to some already visited
			for(int i=0;i<unvisitedEdges.size();i++){
				//printf("%d %d\n", unvisitedEdges[i].v1, unvisitedEdges[i].v2);  
				if (unvisitedEdges[i].v1 == currVertex){// && unvisitedVertices.count(edges.edges[i].v2) == 1){ //if begin at currVertex and ends in a unvisitedVertex
					currVertex = unvisitedEdges[i].v2;
					//printf("found visited %d\n", currVertex);
					edgesUsed->push_back(unvisitedEdges[i]);
					unvisitedEdges.erase(unvisitedEdges.begin()+i);
					//unvisitedVertices.erase(currVertex);
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
		//break;
	}
	path->push_back(currVertex);
	std::vector<int>* aux = removeEqual(*path);
	aux->push_back(initialVertex);
	allPaths.push_back(aux);

	//shortcut
	std::vector<Edge>* edgesInPath = new std::vector<Edge>();
	std::vector<int>* shortcutPath;	
	//shortcutPath->push_back(initialVertex);
	for(int k=0;k<allPaths.size();++k){
		shortcutPath = allPaths[k];	
		for(int i=0;i<shortcutPath->size()-1;++i){
			//printf("%d %d\n", *(shortcutPath->begin()+i), *(shortcutPath->begin()+i+1) ); 
			for(int j=0;j<allEdges.numEdges;++j){
				if (allEdges.edges[j].v1 == *(shortcutPath->begin()+i) && allEdges.edges[j].v2 == *(shortcutPath->begin()+i+1)){
					edgesInPath->push_back(allEdges.edges[j]);
					break;
				}
			}
		}
	}

//	for(int k=0;k<edgesUsed->size();++k)
//		printf("{%d, %d, %d} ", (*(edgesUsed->begin()+k)).v1, (*(edgesUsed->begin()+k)).v2, (*(edgesUsed->begin()+k)).cost);
//	printf("\n");


//	for(int i=0;i<path.size();i++)
//		printf("%d ", path[i]);
//	printf("\n");
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
	edgesReturn.numEdges = 	graph.n * (graph.n-1);
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
	assert(numAdded == edgesReturn.numEdges);
	return edgesReturn;
}

int main(int argc, char* argv[])
{
    srand (time(NULL));

	//Número de vértices, TEM que ser par.
	int numVertices = 10;

	//Instancia-se o objeto passando-se o número de vértices desejado.
	Matching *M = new Matching(numVertices);

    Edge* edges = generateRandomCompleteGraph(numVertices);
    //Edge* edges = generateCompletePCost(numVertices, 20);

    int numEdges = (numVertices)*(numVertices-1)/2; //somatorio

    int costsSum = 0;

    for(int i=0;i<numEdges;++i){
	costsSum += edges[i].cost;	
    }

    Edges allEdges;
    allEdges.edges = edges;
    allEdges.numEdges = numEdges;

/* testing a certain graph
    for(int i=0;i<numEdges;++i){
	int x1 = edges[i].v1;
	int x2 = edges[i].v2;
	if(x1 == 0 && x2 == 8)
		edges[i].cost = 1;
	if(x1 == 1 && x2 == 2)
		edges[i].cost = 2;
	if(x1 == 2 && x2 == 9)
		edges[i].cost = 3;
	if(x1 == 2 && x2 == 3)
		edges[i].cost = 3;
	if(x1 == 4 && x2 == 5)
		edges[i].cost = 3;
	if(x1 == 6 && x2 == 10)
		edges[i].cost = 2;
	if(x1 == 7 && x2 == 8)
		edges[i].cost = 2;
	if(x1 == 3 && x2 == 11)
		edges[i].cost = 2;
    }
*/

    addEdges(allEdges, M);

    //M->SolveMinimumCostPerfectMatching();

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
//    printf("All Edges:\n");
//    printEdges(allEdges);
//    printf("Matching:\n");
//    printEdgesInMatching(edges, numEdges, M);
//    printf("MST:\n");
//    Edges mst = getMST(allEdges, 100);
//    printEdgesSimple(mst);

    //getMinCover(allEdges, 0);

    //printEdges(edges, numEdges);

    //variables used for the MST class
    /*
	The structure used in the MST class and the Matching class are completly different.
	MST code assumes I will have vertices numbered from 0 to m, such that m is the number of vertices.
	That means that I might have to do some translation back and forth to make it into this notation.
	I will have to figure out a way to extract the solution out of it.
    */


    
	int numTrees;
	std::cout << "Digite o número de árvores procuradas: ";
	std::cin >> numTrees;

	int numEdgesNeeded = numVertices - numTrees;


	std::vector<Edge> solution;
	int i = -1;

	int minLimit = -2*costsSum;
	int maxLimit = 0;

	for(int i=0;i<numEdges;++i){
		if (edges[i].cost > maxLimit)
			maxLimit = edges[i].cost;
	}

	printf("min %d max %d", minLimit, maxLimit);


	//do a first pass with the costs as they are
	solution.clear();
//	Edges mst = getMST(allEdges, 0);
//	for(int k=0;k<mst.numEdges;++k)
//		solution.push_back(mst.edges[k]);
//	Edges minimumCover = getMinCover(allEdges, 0);
//	for(int k=0;k<minimumCover.numEdges;++k)
//		solution.push_back(minimumCover.edges[k]);
//	printf("num trees: %d\n", numVertices - solution.size()); 
	int currCost = (minLimit+maxLimit)/2;
	int currTrees, zeroCostEdges;
	while(true){
	//if (!(numVertices - solution.size() <= numTrees)){ //if solution not found, continue
		//for(i=0;i<numEdges;++i){
		   solution.clear();
		   printf("min %d max %d\n", minLimit, maxLimit);
		   printf("lambda: %d\n", currCost);
		   //currCost = allEdges.edges[i].cost;
		   Edges mst = getMST(allEdges, currCost);
		   for(int k=0;k<mst.numEdges;++k)
		   	solution.push_back(mst.edges[k]);
		   printf("MST:\n");
		   printEdges(mst);
		   Edges minimumCover = getMinCover(allEdges, currCost);
		   for(int k=0;k<minimumCover.numEdges;++k)
		   	solution.push_back(minimumCover.edges[k]);
		   printf("Min Cover:\n");
		   printEdges(minimumCover);
		
		   Edges auxedges;
		   auxedges.numEdges = solution.size();
		   auxedges.edges =  &solution.front();
	           printEdges(auxedges);

		   //printf("MST size %d Min Cover size %d solution size %d\n", mst.numEdges, minimumCover.numEdges, solution.size());
		   //print solution
		   currTrees = numVertices - solution.size();
		   zeroCostEdges = 0;
		   int numEdgesIncident[numVertices];
		   for(i=0;i<numVertices;++i){
			numEdgesIncident[i] = 0;
		   }

	       	   for(i=0;i<solution.size();++i){
			numEdgesIncident[solution[i].v1]++;
			numEdgesIncident[solution[i].v2]++;
		   }
	       	   for(i=0;i<solution.size();++i){
				if (solution[i].cost == currCost && numEdgesIncident[solution[i].v1] > 1 && numEdgesIncident[solution[i].v2] > 1){
					numEdgesIncident[solution[i].v1]--;
					numEdgesIncident[solution[i].v2]--;
					zeroCostEdges++;
				}	
		   }

		   printf("num trees: %d-%d\n", currTrees, currTrees+zeroCostEdges); 
		   if (numTrees >= currTrees && numTrees <= currTrees+zeroCostEdges){
			printf("found\n");
			break; //found solution
		   }else if (currTrees+zeroCostEdges > numTrees){ //I have less trees than expected
			printf("more trees\n");
			minLimit = currCost;
			currCost += (maxLimit - currCost)/2;
		   }else{ //I have more trees than expected
			printf("less trees\n");
			maxLimit = currCost;
			currCost -= (currCost - minLimit)/2;
		   }
		   //if (numVertices - solution.size() <= numTrees){ //tree number is less than or equal to the amount expected
		   //     printf("stop\n");
		   //	break;
		   //}
		//}
	}

	if (currTrees < numTrees){ //I have less trees than expected
		/*
		This means that I should remove some edges. If I remove n edges without having a vertex uncovered,
		such that n is the difference between the amount of trees I have and the expected amount,
		I will have the amount of trees expected. These edges will have cost equal to the last cost calculated
		(zero-cost edges in the paper). They should also be the edges of smallest cost that I should be able
		to remove.rne Stroustrup: Why I Created C++
		*/
		int removed = 0;
		int expected = (numTrees - currTrees);
		int numEdgesIncident[numVertices];
	       	for(i=0;i<numVertices;++i){
			numEdgesIncident[i] = 0;
		}
	       	for(i=0;i<solution.size();++i){
			numEdgesIncident[solution[i].v1]++;
			numEdgesIncident[solution[i].v2]++;
		}
		while (removed < expected){
			for(i=0;i<solution.size();++i){
				if (solution[i].cost == currCost && numEdgesIncident[solution[i].v1] > 1 && numEdgesIncident[solution[i].v2] > 1){
					numEdgesIncident[solution[i].v1]--;
					numEdgesIncident[solution[i].v2]--;
					removed++;
					solution.erase(solution.begin()+i);
					break;
				}	
			}
		}	
	}

	printf("num trees: %d\n", numVertices - solution.size()); 
	for(int k=0;k<solution.size();++k)
		printf("{%d, %d, %d} ", solution[k].v1, solution[k].v2, solution[k].cost);
	
	printf("\n");

	Edges toDouble;
	toDouble.edges = (Edge*)&solution.front();
	toDouble.numEdges = solution.size();

	printEdges(eulerianCircuit(doubleEdges(toDouble), doubleEdges(allEdges)));//	TODO: I need to do this for every tree in the graph
	
	delete M;

	return 0;
}

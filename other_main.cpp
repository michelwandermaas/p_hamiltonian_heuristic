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
		printf("%d %d %lf\n", edges[i].v1, edges[i].v2, edges[i].cost);
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
    std::set<int> vertices;
	for(int i=0;i<edges.numEdges;i++){
		vertices.insert(edges.edges[i].v1);
		vertices.insert(edges.edges[i].v2);
	}

    std::set<int> usedVertices(vertices);
    std::vector<Edge> newVertexToEdge; //newVertex to the equivalent edge
    std::set<int> addedVertices;
    int numVertex = vertices.size(); //assume 0 based, with all vertices
    std::vector<Edge> zeroCostEdges;

    std::vector<Edge> allEdges;
    for(int i=0;i<edges.numEdges;i++){
        Edge aux;
        aux.v1 = edges.edges[i].v1;
        aux.v2 = numVertex;
        aux.cost = edges.edges[i].cost / 2.0;

        Edge aux1;
        aux1.v1 = edges.edges[i].v2;
        aux1.v2 = numVertex+1;
        aux1.cost = edges.edges[i].cost / 2.0;

	    newVertexToEdge.push_back(edges.edges[i]);	
	    newVertexToEdge.push_back(edges.edges[i]);	

        Edge aux2;
        aux2.v1 = numVertex;
        aux2.v2 = numVertex + 1;
        aux2.cost = 0.0;

        zeroCostEdges.push_back(aux2);

        //add vertices to set
        addedVertices.insert(numVertex); 
        usedVertices.insert(numVertex++); 
        addedVertices.insert(numVertex); 
        usedVertices.insert(numVertex++); 

        //add edges to allEdges
        allEdges.push_back(aux);
        allEdges.push_back(aux1);
        allEdges.push_back(aux2);
	}

    if (((usedVertices.size()) % 2) == 1){ //add dummy edges if need be
	    for(int i=0;i<usedVertices.size();i++){
           Edge aux;
           aux.v1 = i;
           aux.v2 = usedVertices.size();
           aux.cost = 0;
           allEdges.push_back(aux); 
        }
    }
    
    Matching *M = new Matching((usedVertices.size()) + ((usedVertices.size()) % 2));
	for(int i=0;i<allEdges.size();i++){
		    //printf("{%d, %d, %d}\n", allEdges[i].v1,allEdges[i].v2, allEdges[i].cost);
			M->AddEdge(allEdges[i].v1, allEdges[i].v2, allEdges[i].cost); //add the edges that have both endpoints in Vp
    }

	M->SolveMinimumCostPerfectMatching();

    std::vector<Edge> twoFactor;

    std::vector<int> verticesVector(vertices.begin(), vertices.end());
    std::vector<int> addedVector(addedVertices.begin(), addedVertices.end());

	for(int i=0;i<verticesVector.size();i++){
	    for(int j=0;j<addedVector.size();j++){
            //printf("%d %d\n",  verticesVector[i], addedVector[j]-edges.numEdges);
		    if(M->IsInMatching(verticesVector[i], addedVector[j])){
                Edge x = newVertexToEdge[addedVector[j]-vertices.size()];
                //printf("%d %d\n",  addedVector[j], addedVector[j]-vertices.size());
                //printf("%d %d\n", verticesVector[i], addedVector[j]);
		        printf("{%d, %d, %lf}\n", x.v1,x.v2, x.cost);
                twoFactor.push_back(newVertexToEdge[addedVector[j]-vertices.size()]);
            }
        }
    }

    Edges circuit;
	circuit.edges = &twoFactor.front();
	circuit.numEdges = twoFactor.size();

    printEdgesInMatching(circuit.edges, circuit.numEdges, M);

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

    if (allEdges.numEdges == 0){
        printf("Error.\n");
        return 0;
    }

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

    return 0;
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

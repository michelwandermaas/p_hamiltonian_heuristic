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

int compareEdges(const void * a, const void * b)
{
  if ( (*(Edge*)a).cost <  (*(Edge*)b).cost ) return -1;
  if ( (*(Edge*)a).cost == (*(Edge*)b).cost ) return 0;
  if ( (*(Edge*)a).cost >  (*(Edge*)b).cost ) return 1;
}


void orderEdges(Edge* edges, int numEdges){
	qsort (edges, numEdges, sizeof(Edge), compareEdges);
}

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
    	double costSum = 0;
	for(int i=0;i<edges.numEdges;++i){
		printf("{%d, %d, %.0lf} ", edges.edges[i].v1,edges.edges[i].v2, edges.edges[i].cost);
		costSum += edges.edges[i].cost;
	}
	printf("\n");
	printf("total cost: %.0lf\n", costSum);
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


Edges getMST(Edges edges){
	/*
	 *      Vn = vertices incident to a nonpositive weight edge
     * 		En = edges with both endpoints in Vn
     * 		G1 = MST in Vn over En
	 */
	Edges result;
	std::set<int> Vn; //Vn = vertices incident to a nonpositive weight edge
	std::vector<Edge> En; //En = edges with both endpoints in Vn
	for(int i=0;i<edges.numEdges;i++){
		Vn.insert(edges.edges[i].v1);
		Vn.insert(edges.edges[i].v2);
	}

	for(int i=0;i<edges.numEdges;i++){
		En.push_back(edges.edges[i]);
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

	Edge* edgesResult = new Edge[numVerticesMST-1];
	result.edges = edgesResult;
	int numEdges = 0;

	for(int i = 0; i < numVerticesMST-1; i++){
		edgesResult[numEdges++] = En[mst.sol[i]];
	}

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

Edges removeEqualEdge(Edges vector){
	std::vector<Edge>* returnVector = new std::vector<Edge>();
	bool found;
	for(int i=0;i<vector.numEdges;++i){
		found = false;
		for(int j=0;j<returnVector->size();++j){
			if (*(returnVector->begin()+j) == vector.edges[i]){
				found = true;
				break;
			}
		}
		if (!found){
			returnVector->push_back(vector.edges[i]);
		}
	}
	Edges ret;
	ret.edges=&returnVector->front();
	ret.numEdges=returnVector->size();
	return ret;
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

    std::vector<int> verticesVector(vertices.begin(), vertices.end()); //vector with original vertices

    int numVerticesOriginal = vertices.size();

	for(int i=0;i<numVerticesOriginal;i++){ //add dummy original vertices
        vertices.insert(i+numVerticesOriginal);
    }

    std::set<int> usedVertices(vertices);
    std::vector<Edge> newVertexToEdge; //newVertex to the equivalent edge
    std::set<int> addedVertices;
    int numVertex = vertices.size(); //assume 0 based, with all vertices
    std::vector<Edge> zeroCostEdges;

    std::vector<Edge> allEdges;
    for(int i=0;i<edges.numEdges;i++){

        //add first edge

        Edge aux;
        aux.v1 = edges.edges[i].v1;
        aux.v2 = numVertex;
        aux.cost = edges.edges[i].cost;

        allEdges.push_back(aux);

        //add edge to dummy
        Edge aux3;
        aux3.v1 = edges.edges[i].v1+numVerticesOriginal;
        aux3.v2 = numVertex;
        aux3.cost = edges.edges[i].cost;

        allEdges.push_back(aux3);


        //add second edge

        Edge aux1;
        aux1.v1 = edges.edges[i].v2;
        aux1.v2 = numVertex+1;
        aux1.cost = edges.edges[i].cost;

        allEdges.push_back(aux1);

        //add edge to dummy
        Edge aux4;
        aux4.v1 = edges.edges[i].v2+numVerticesOriginal;
        aux4.v2 = numVertex+1;
        aux4.cost = edges.edges[i].cost;

        allEdges.push_back(aux4);



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
    //printf("size V %d\n", (usedVertices.size()) + ((usedVertices.size()) % 2));
    //printf("size %d\n", allEdges.size());
	for(int i=0;i<allEdges.size();i++){
            //printf("%d\n", i);
		    //printf("{%d, %d, %d}\n", allEdges[i].v1,allEdges[i].v2, allEdges[i].cost);
            //std::cout << std::flush;
			M->AddEdge(allEdges[i].v1, allEdges[i].v2, allEdges[i].cost); //add the edges that have both endpoints in Vp
    }
    //assert(0==1); 
    //std::cout << "solving" << std::endl << std::flush;
	M->SolveMinimumCostPerfectMatching();
    //std::cout << "solved" << std::endl << std::flush;

    std::vector<Edge>* twoFactor = new std::vector<Edge>();

    std::vector<int> addedVector(addedVertices.begin(), addedVertices.end());

	for(int i=0;i<verticesVector.size();i++){
	    for(int j=0;j<addedVector.size();j++){
            //printf("%d %d\n",  verticesVector[i], addedVector[j]-edges.numEdges);
		    if(M->IsInMatching(verticesVector[i], addedVector[j]) || M->IsInMatching(verticesVector[i]+numVerticesOriginal, addedVector[j])){
                //Edge x = newVertexToEdge[addedVector[j]-vertices.size()];
                //printf("%d %d\n",  addedVector[j], addedVector[j]-vertices.size());
                //printf("%d %d\n", verticesVector[i], addedVector[j]);
		        //printf("{%d, %d, %lf}\n", x.v1,x.v2, x.cost);
                twoFactor->push_back(newVertexToEdge[addedVector[j]-vertices.size()]);
            }
        }
    }

    Edges solution;
	solution.edges = &twoFactor->front();
	solution.numEdges = twoFactor->size();

    //printEdges(solution);

    return solution;
}

int countCycles(Edges edges){
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
			if ((unvisitedEdges[i].v1 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v2) == 1) || (unvisitedEdges[i].v2 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v1) == 1) ){ //if begin at currVertex and ends in a unvisitedVertex
				if (unvisitedEdges[i].v1 == currVertex)
					currVertex = unvisitedEdges[i].v2;
				else
					currVertex = unvisitedEdges[i].v1;
				edgesUsed->push_back(unvisitedEdges[i]);
				unvisitedEdges.erase(unvisitedEdges.begin()+i);
				unvisitedVertices.erase(currVertex);
				found = true;
				break;
			}
		}
		if (!found){ //if no unvisited edge is available, go to some already visited
			for(int i=0;i<unvisitedEdges.size();i++){
				if ((unvisitedEdges[i].v1 == currVertex || unvisitedEdges[i].v2 == currVertex)){// && unvisitedVertices.count(edges.edges[i].v2) == 1){ //if begin at currVertex and ends in a unvisitedVertex
					if (unvisitedEdges[i].v1 == currVertex)
						currVertex = unvisitedEdges[i].v2;
					else
						currVertex = unvisitedEdges[i].v1;
					edgesUsed->push_back(unvisitedEdges[i]);
					unvisitedEdges.erase(unvisitedEdges.begin()+i);
					found = true;
					break;
				}
			}
		}
		if (!found){ //that means that that cycle is done, move to next one
			//path->push_back(firstVertexPath);
			std::vector<int>* aux = removeEqual(*path);
			path = aux;
			path->push_back(initialVertex);
			allPaths.push_back(path);
			path = new std::vector<int>;
			initialVertex = *unvisitedVertices.begin();
			currVertex = initialVertex;
			unvisitedVertices.erase(currVertex);
		}
	}

    return allPaths.size();
}

int countCyclesMinVertices(Edges edges, int minVertices){ //this count might be wrong
    //for every cycle with more than (minVertices*2), result += (floor(numCycles / minVertices) - 1)
    int count = 0;

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
			if ((unvisitedEdges[i].v1 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v2) == 1) || (unvisitedEdges[i].v2 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v1) == 1) ){ //if begin at currVertex and ends in a unvisitedVertex
				if (unvisitedEdges[i].v1 == currVertex)
					currVertex = unvisitedEdges[i].v2;
				else
					currVertex = unvisitedEdges[i].v1;
				edgesUsed->push_back(unvisitedEdges[i]);
				unvisitedEdges.erase(unvisitedEdges.begin()+i);
				unvisitedVertices.erase(currVertex);
				found = true;
				break;
			}
		}
		if (!found){ //if no unvisited edge is available, go to some already visited
			for(int i=0;i<unvisitedEdges.size();i++){
				if ((unvisitedEdges[i].v1 == currVertex || unvisitedEdges[i].v2 == currVertex)){// && unvisitedVertices.count(edges.edges[i].v2) == 1){ //if begin at currVertex and ends in a unvisitedVertex
					if (unvisitedEdges[i].v1 == currVertex)
						currVertex = unvisitedEdges[i].v2;
					else
						currVertex = unvisitedEdges[i].v1;
					edgesUsed->push_back(unvisitedEdges[i]);
					unvisitedEdges.erase(unvisitedEdges.begin()+i);
					found = true;
					break;
				}
			}
		}

		if (!found){ //that means that that cycle is done, move to next one
			path->push_back(currVertex);
            		//printf("cycle size %d\n", path->size());
            		if (path->size() >= minVertices * 2)
                		count += (path->size() / minVertices) - 1; 
			std::vector<int>* aux = removeEqual(*path);
			aux->push_back(initialVertex);
			allPaths.push_back(aux);
			path = new std::vector<int>;
			initialVertex = *unvisitedVertices.begin();
			currVertex = initialVertex;
			unvisitedVertices.erase(currVertex);
		}
	}
    if (path->size() >= minVertices * 2)
        count += (path->size() / minVertices) - 1; 

    return count;
}

std::vector<std::vector<int>*>* breakPath(std::vector<int>* path, int numBreak, int minVertices){
	std::vector<std::vector<int>*>* allPaths = new std::vector<std::vector<int>*>();
	int breakPos = 0;
	std::vector<int>* newPath;
	for (int i=0;i<numBreak;i++){
		newPath = new std::vector<int>();
		for(int j=breakPos; j<breakPos+minVertices; j++){
			newPath->push_back((*path)[j]);	
		}
		newPath->push_back((*path)[breakPos]);
		allPaths->push_back(newPath);
		breakPos += minVertices;	
	}
	if (breakPos < path->size()){
		newPath->pop_back();
		for(int j=breakPos; j<path->size(); j++){
			newPath->push_back((*path)[j]);	
		}
		newPath->push_back((*path)[breakPos-minVertices]);
	}
	return allPaths;
}

Edges breakCycles(Edges edges, Edges allEdges, int minVertices, int numToBreak){
    //similar to the one above and the euclidean circuit, but break as it goes
	std::set<int> unvisitedVertices;	
	std::vector<Edge> unvisitedEdges;
	for(int i=0;i<edges.numEdges;i++){
		unvisitedEdges.push_back(edges.edges[i]);	
		unvisitedVertices.insert(edges.edges[i].v1);
		unvisitedVertices.insert(edges.edges[i].v2);
	}

	std::vector<Edge>* auxVector = removeEqualEdge(unvisitedEdges);
	unvisitedEdges = *auxVector;
	int initialVertex = *unvisitedVertices.begin();
	int currVertex = initialVertex;
	unvisitedVertices.erase(currVertex);
	std::vector<int>* path;
	path = new std::vector<int>;
	bool found;
	int count;
	std::vector<Edge>* edgesUsed = new std::vector<Edge>();
	std::vector<std::vector<int>*> allPaths;
	allPaths.push_back(path);
	int firstVertexPath = currVertex;
	while (unvisitedVertices.size() > 0){	
		found = false;
		path->push_back(currVertex);
		for(int i=0;i<unvisitedEdges.size();i++){
			//printf("%d %d -> %d\n", unvisitedEdges[i].v1, unvisitedEdges[i].v2, currVertex); 
			if ((unvisitedEdges[i].v1 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v2) == 1) || (unvisitedEdges[i].v2 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v1) == 1) ){ //if begin at currVertex and ends in a unvisitedVertex
				if (unvisitedEdges[i].v1 == currVertex)
					currVertex = unvisitedEdges[i].v2;
				else
					currVertex = unvisitedEdges[i].v1;
				edgesUsed->push_back(unvisitedEdges[i]);
				unvisitedEdges.erase(unvisitedEdges.begin()+i);
				unvisitedVertices.erase(currVertex);
				//printf("found\n");
				found = true;
				break;
			}
		}
		if (!found){ //if no unvisited edge is available, go to some already visited
			for(int i=0;i<unvisitedEdges.size();i++){
				if ((unvisitedEdges[i].v1 == currVertex || unvisitedEdges[i].v2 == currVertex)){// && unvisitedVertices.count(edges.edges[i].v2) == 1){ //if begin at currVertex and ends in a unvisitedVertex
					if (unvisitedEdges[i].v1 == currVertex)
						currVertex = unvisitedEdges[i].v2;
					else
						currVertex = unvisitedEdges[i].v1;
					edgesUsed->push_back(unvisitedEdges[i]);
					unvisitedEdges.erase(unvisitedEdges.begin()+i);
					found = true;
					break;
				}
			}
		}
		if (!found){ //that means that that cycle is done, move to next one
			//printf("new path\n");
			//path->push_back(firstVertexPath);
			std::vector<int>* aux = removeEqual(*path);
			path = aux;
			if (numToBreak > 0 && path->size() >= minVertices * 2){
                		count = fmin((path->size() / minVertices) - 1, numToBreak); 
				numToBreak -= count;
				//printf("break %d cycles\n", count);
				//when I use count + 1, it means that I want to add count cycles, but the result should be count + 1 cycles, since it is currently just one
				std::vector<std::vector<int>*>* brokenPaths = breakPath(path, count+1, 3);	
				assert(brokenPaths->size() == count+1);
				for (int k=0;k<count+1; k++){
					//for(int p=0;p<(*brokenPaths)[k]->size();p++)
						//printf("%d ", (*brokenPaths)[k]->at(p));	
					//printf("\n");
					allPaths.push_back((*brokenPaths)[k]);
				}
			}else{
				path->push_back(initialVertex);
				allPaths.push_back(path);
			}
			path = new std::vector<int>;
			initialVertex = *unvisitedVertices.begin();
			currVertex = initialVertex;
			firstVertexPath = currVertex;
			unvisitedVertices.erase(currVertex);
		}
	}
	path->push_back(currVertex);
	std::vector<int>* aux = path;
	aux->push_back(initialVertex);
	allPaths.push_back(aux);

	//shortcut
	std::vector<Edge>* edgesInPath = new std::vector<Edge>();
	std::vector<int>* shortcutPath;	
	int j;
	//printf("\n\n");
	allPaths.erase(allPaths.begin()); //delete first one, do not ask why
	for(int k=0;k<allPaths.size();++k){
		shortcutPath = allPaths[k];	
		//for(int i=0;i<shortcutPath->size();++i)
		//	printf("%d ", *(shortcutPath->begin()+i) );
		//printf("\n");
		for(int i=0;i<shortcutPath->size()-1;++i){
			for(j=0;j<allEdges.numEdges;++j){
				if (allEdges.edges[j].v1 == *(shortcutPath->begin()+i) && allEdges.edges[j].v2 == *(shortcutPath->begin()+i+1)){
					edgesInPath->push_back(allEdges.edges[j]);
					break;
				}
			}
			
			if (j == allEdges.numEdges){
				for(int k=0;k<shortcutPath->size();++k)
					printf("%d ", *(shortcutPath->begin()+k) );
				printf("\ncould not find edge %d %d\n", *(shortcutPath->begin()+i), *(shortcutPath->begin()+i+1) );
				assert(0==1); //it did not find the correct edge
			}
			
		}
	}

	//edgesInPath = removeEqualEdge(*edgesInPath);
	Edges circuit;
	circuit.edges = &edgesInPath->front();
	circuit.numEdges = edgesInPath->size();
	return circuit;
}

Edges addMSTEdges(Edges two_factor_sol, Edges mst, int cyclesToReduce){
	/*beginning is just like the countCycles.
	then with the paths set up, I will go over the mst edges in cost order
	and I will add edges (doubled) that are between two cycles, and merge the paths until I have
	merged the right number of cycles.
	*/
	std::set<int> unvisitedVertices;	
	std::vector<Edge> unvisitedEdges;
	for(int i=0;i<two_factor_sol.numEdges;i++){
		unvisitedEdges.push_back(two_factor_sol.edges[i]);	
		unvisitedVertices.insert(two_factor_sol.edges[i].v1);
		unvisitedVertices.insert(two_factor_sol.edges[i].v2);
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
			if ((unvisitedEdges[i].v1 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v2) == 1) || (unvisitedEdges[i].v2 == currVertex && unvisitedVertices.count(unvisitedEdges[i].v1) == 1) ){ //if begin at currVertex and ends in a unvisitedVertex
				if (unvisitedEdges[i].v1 == currVertex)
					currVertex = unvisitedEdges[i].v2;
				else
					currVertex = unvisitedEdges[i].v1;
				edgesUsed->push_back(unvisitedEdges[i]);
				unvisitedEdges.erase(unvisitedEdges.begin()+i);
				unvisitedVertices.erase(currVertex);
				found = true;
				break;
			}
		}
		if (!found){ //if no unvisited edge is available, go to some already visited
			for(int i=0;i<unvisitedEdges.size();i++){
				if ((unvisitedEdges[i].v1 == currVertex || unvisitedEdges[i].v2 == currVertex)){// && unvisitedVertices.count(edges.edges[i].v2) == 1){ //if begin at currVertex and ends in a unvisitedVertex
					if (unvisitedEdges[i].v1 == currVertex)
						currVertex = unvisitedEdges[i].v2;
					else
						currVertex = unvisitedEdges[i].v1;
					edgesUsed->push_back(unvisitedEdges[i]);
					unvisitedEdges.erase(unvisitedEdges.begin()+i);
					found = true;
					break;
				}
			}
		}
		if (!found){ //that means that that cycle is done, move to next one
			//path->push_back(firstVertexPath);
			std::vector<int>* aux = removeEqual(*path);
			path = aux;
			path->push_back(initialVertex);
			allPaths.push_back(path);
			path = new std::vector<int>;
			initialVertex = *unvisitedVertices.begin();
			currVertex = initialVertex;
			unvisitedVertices.erase(currVertex);
		}
	}

	Edges result;
	result.numEdges = two_factor_sol.numEdges;// + (cycleToReduce*2)
	Edge* resultEdges = new Edge[result.numEdges + (cyclesToReduce*2)];
	result.edges = resultEdges;

	orderEdges(two_factor_sol.edges, two_factor_sol.numEdges); //order edges to add

	for(int i=0;i<two_factor_sol.numEdges;++i)
		resultEdges[i] = two_factor_sol.edges[i];


	//TODO: this can be done efficently, not this way
	//it should be done in a manner similar to the MST idea to avoid cycles
	int numAdded = 0;
	for(int i=0;i<mst.numEdges;++i){
		int v1 = mst.edges[i].v1;
		int v2 = mst.edges[i].v2;
		std::vector<int>* path1;
		std::vector<int>* path2;
		int path2_pos = 0;
		for(int j=0;j<allPaths.size();++j){
			for(int k=0;k<allPaths[j]->size();++k){
				if ((*allPaths[j])[k] == v1)
					path1 = allPaths[j];		
				if ((*allPaths[j])[k] == v2){
					path2 = allPaths[j];	
					path2_pos = k;
				}	
			}
		}
		if (path1 != path2){ //they are in different cycles
			numAdded++;
			printf("%d\n", result.numEdges);
			result.edges[result.numEdges++] = mst.edges[i];
			Edge invertEdge;
			invertEdge.v1 = v2;
			invertEdge.v2 = v1;
			invertEdge.cost = mst.edges[i].cost;
			result.edges[result.numEdges++] = invertEdge;
			//merge paths
			for(int k=0;k<path2->size();++k){
				path1->push_back((*path2)[k]);
			}
			allPaths.erase(allPaths.begin()+path2_pos);
		}	
		if (numAdded == cyclesToReduce)
			break;
	}
	
	return result;
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

    Edges factor2 = removeEqualEdge(two_factor(allEdges)); 

    int numCycles2Factor = countCycles(factor2);

    printf("2-factor %d cycles\n", numCycles2Factor);

    printEdges(factor2);

    if (numCycles2Factor == numTrees){
        std::cout << "Solução ótima encontrada." << std::endl;
        printEdges(factor2);
        return 0;
    }else if (numCycles2Factor < numTrees){
    /*
    Case 3: q < p. It means that we need to increase the number of cycles by removing edges.
        As long as there are at least (p-q) cycles with 6 or more vertices, we can do this by removing two edges (a,b)
        and (c,d) and then adding edges (a,c) and (b,d), for example. Therefore there's a limit to how many trees we
        can guarantee in this case. We guarantee a 2-approximation this way.
    */
        //std::cout << "p < q" << std::endl;
        int numCyclesToAdd = countCyclesMinVertices(factor2, 3);
        //printf("can add up to %d cycles\n", numCyclesToAdd);
        if (numCyclesToAdd < (numTrees - numCycles2Factor)){
            std::cout << "Não é possível encontrar solução para esse número de árvores." << std::endl;
            return 0;
        }
	printf("Solution with %d cycles\n", numTrees);
	factor2 = breakCycles(factor2, doubleEdges(allEdges), 3, numCyclesToAdd);
	assert(countCycles(factor2) == numTrees);
        printEdges(factor2);

    }else{
    /*
    Case 2: q > p. It means that we need to reduce the number of cycles by adding edges.
        We find a T with is the MST of G. We need to find (q-p) edges in T which we will double and then
        add to F, such that each of these edges connects 2 cycles in F. We guarantee a 3-approximation this way.
    */
        std::cout << "p > q" << std::endl;
	Edges mst = getMST(allEdges);
	Edges sol = addMSTEdges(factor2, mst, numCycles2Factor - numTrees);
        printEdges(sol);
    }
    return 0;
}

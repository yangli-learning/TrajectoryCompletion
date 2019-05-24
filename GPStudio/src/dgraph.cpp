/* Directed Graphs
 * ----------------------------------------------------------------------------
 * Author:  Shane Saunders
 * modified by yang li
 */
#include <cstdio>
#include "dgraph.h"
#include <iostream>
/*--- DGraph ----------------------------------------------------------------*/

/* --- Constructor ---
 * Creates a DGraph object containing n vertices.
 */
DGraph::DGraph(int n)
{
    nVertices = n;
    
    vertices = new DGraphVertex[n];
    initVertices();
}

/* --- Destructor ---
 */
DGraph::~DGraph()
{
    clear();
    delete [] vertices;
}

/* --- clear() ---
 * Clears all edges from the graph.
 */
void DGraph::clear()
{
    DGraphEdge *edge, *nextEdge;
    for(int i = 0; i < nVertices; i++) {
        edge = vertices[i].outHead;

        while(edge) {
            nextEdge = edge->nextOut;
            delete edge;
            edge = nextEdge;
        }
    }
    initVertices();
}

void DGraph::initVertices()
{
    for(int i = 0; i < nVertices; i++) {
        vertices[i].outHead = vertices[i].outTail = 0;
        vertices[i].inHead = vertices[i].inTail = 0;
        vertices[i].outSize = vertices[i].inSize = 0;
    }
}
    
/* --- addNewEdge() ---
 * Adds a new edge from vertex 'source' to vertex 'target' with
 * with a corresponding distance of dist.
 */
void DGraph::addNewEdge(int source, int target, long dist)
{
    DGraphEdge *newEdge = new DGraphEdge;
    newEdge->source = source;
    newEdge->target = target;
    newEdge->dist = dist;
    newEdge->nextOut = NULL;
    newEdge->nextIn = NULL;
    
    DGraphVertex *vertex = &vertices[source];
    if(vertex->outTail) {
        vertex->outTail->nextOut = newEdge;
    }
    else {
        vertex->outHead = newEdge;
    }
    vertex->outTail = newEdge;
    vertex->outSize++;

    vertex = &vertices[target];
    if(vertex->inTail) {
        vertex->inTail->nextIn = newEdge;
    }
    else {
        vertex->inHead = newEdge;
    }
    vertex->inTail = newEdge;
    vertex->inSize++;
};
long int DGraph::getEdgeWeight (int src, int dest){
    const DGraphEdge *edge = vertices[src].outHead;
    while(edge) {
        if(edge->target == dest) {
            return edge->dist;
        }
        edge = edge->nextOut;
    }
    std::cout << "such reverse.."<< std::endl;
    const DGraphEdge *edgeR = vertices[dest].outHead;
    while(edgeR) {
        if(edgeR->target == src) {
            return edgeR->dist;
        }
        edgeR = edgeR->nextOut;
    }
    std::cout << "no such edge "<< src<< ", " << dest << std::endl;
    return 0;
}

bool DGraph::edgeExists(int v, int w) const
{
    /* Scan all existing edges from v to determine whether an edge to w exists.
     */
    const DGraphEdge *edge = vertices[v].outHead;
    while(edge) {
        if(edge->target == w) return true;
        edge = edge->nextOut;
    }
    return false;
}

void DGraph::removeEdge(int v,int w) {
	DGraphEdge *edge = vertices[v].outHead;
	if (edge==NULL) return;
	if (edge->target == w){
		// edge to be removed is head
		vertices[v].outHead = edge->nextOut;
	}else{
		// find the edge before the edge to be removed
		while(edge->nextOut){
			if ( (edge->nextOut)->target == w)
				break;
			edge = edge->nextOut;
		}
		if (vertices[v].outTail == edge->nextOut)
			vertices[v].outTail = edge; // edge to be removed is tail
		edge->nextOut = (edge->nextOut)->nextOut;
	}
	vertices[v].outSize --;
	DGraphEdge *toberemoved;
	edge = vertices[w].inHead;
	if (edge->source == v){
		toberemoved = edge;
		vertices[w].inHead = edge->nextIn;
	}else{
		while(edge->nextIn){
			if ((edge->nextIn)->source ==v)
				break;
			edge = edge->nextIn;
		}
		toberemoved = edge->nextIn;
		if (vertices[w].inTail == edge->nextIn)
			vertices[w].inTail = edge;
		edge->nextIn = (edge->nextIn)->nextIn;
	}
	vertices[w].inSize --;
	delete toberemoved;

}
/* --- reachable() ---
 * Test whether all vertices are reachable from the source vertex s.
 */
bool DGraph::reachable(int s) const
{
    int *stack = new int[nVertices];
    int tos = 0;
    
    int *visited = new int[nVertices];
    for(int i = 0; i < nVertices; i++) visited[i] = 0;

    int vertexCount = 0;
    visited[s] = 1;
    stack[tos++] = s;
    DGraphEdge *edge;
    int v, w;
    while(tos) {
        v = stack[--tos];
        vertexCount++;
        edge = vertices[v].outHead;
        while(edge) {
            w = edge->target;
            if(!visited[w]) {
                visited[w] = 1;
                stack[tos++] = w;
            }
            edge = edge->nextOut;
        }
    }

    delete [] stack;
    delete [] visited;
    std::cout << " number of reachable vertices from vertex "<<s<< ": "
         << vertexCount<< std::endl;
    return vertexCount == nVertices;
}


/* --- print() ---
 * Prints a text representation of the graph to the standard output.
 */
void DGraph::print() const
{
    const DGraphEdge *edge;

    printf("Graph (vertex: edge list) = \n");

    for(int i = 0; i < nVertices; i++) {
        printf("%d: ", i);
        edge = vertices[i].outHead;
        while(edge) {
            printf(" %d", edge->target);
            edge = edge->nextOut;
        }
        putchar('\n');
    }

    printf("Graph (vertex: edge{dist} list) = \n");

    for(int i = 0; i < nVertices; i++) {
        printf("%d: ", i);
        edge = vertices[i].outHead;
        while(edge) {
            long int tt = edge->target;
            long int dd = edge->dist;

            printf("%ld {%ld} ", tt, dd);
//            printf("%l %l", edge->target, edge->dist);
            edge = edge->nextOut;
        }
        putchar('\n');
    }
}


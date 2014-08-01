/*
 * Main developer: Nico Van Cleemput
 * 
 * Copyright (C) 2014 Nico Van Cleemput.
 * Licensed under the GNU GPL, read the file LICENSE.txt for details.
 */

/* This program reads simple graphs from standard in and
 * checks whether they are pancyclic. This program is designed for large and dense
 * graphs which have a large automorphism group. For small and sparse graphs or
 * for graphs which have near-to-trivial symmetry, it might be faster to just
 * generate all cycles.
 * This program does not support graphs which contain vertices of degree 2.
 * 
 * 
 * Compile with:
 * 
 *     cc -o pancyclic -O4 pancyclic.c
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>

#include "nauty/nauty.h"

#include "shared/multicode_base.h"
#include "shared/multicode_input.h"
#include "shared/multicode_output.h"

typedef int VERTEXPAIR[2];

/* Nauty worksize */
#define WORKSIZE 50 * MAXM

/** Nauty variables */
int lab[MAXN], ptn[MAXN], orbits[MAXN];
static DEFAULTOPTIONS_GRAPH(options);
statsblk stats;
setword workspace[WORKSIZE];

graph ng[MAXN*MAXM]; /* nauty graph datastructure */
graph ng_canon[MAXN*MAXM]; /* nauty graph datastructure */

set verticesInCycle[MAXM];

int generators[MAXN+1][MAXN/2][MAXN];
int generatorCount[MAXN+1];
boolean generatorsDetermined[MAXN+1];

int vertexCount = 0;
int addedVerticesCount = 0;
int m;

// debugging methods

void printGenerators(int depth, int vertexCount){
    int i, j;
    fprintf(stderr, "Generators:\n");
    for(i = 0; i < generatorCount[depth]; i++){
        for(j = 0; j < vertexCount; j++){
            fprintf(stderr, "%d ", generators[depth][i][j]);
        }
        fprintf(stderr, "\n");
    }
}

// end debugging methods

/**
 * Method which is called each time nauty finds a generator.
 */
void storeGenerators(int count, permutation perm[], nvector orbits[], int numorbits, int stabvertex, int n) {
    memcpy(generators[addedVerticesCount] + generatorCount[addedVerticesCount], perm, sizeof (permutation) * n);

    generatorCount[addedVerticesCount]++;
}

void initNautyRelatedVariables(){
    
    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    options.userautomproc = storeGenerators;

}

inline void prepareNautyCall(){
    generatorCount[addedVerticesCount] = 0;
}

/* This method translates the internal data structure to nauty's dense graph
 * data structure, so the graph can be passed to nauty.
 */
inline void translateGraphToNautyDenseGraph(GRAPH graph, ADJACENCY adj){
    int n, i, j;
    
    n = graph[0][0];
    //the largest graphs we can get contain 2*n vertices (i.e. if the graph is hamiltonian)
    
    if(n > MAXN/2){
        fprintf(stderr, "We only support graphs with up to %d vertices - exiting!\n", MAXN/2);
        exit(EXIT_FAILURE);
    }
    
    m = SETWORDSNEEDED(2*n);
    
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    
    EMPTYGRAPH(ng,m,2*n);
    
    for(i = 1; i <= graph[0][0]; i++){
        for(j = 0; j < adj[i]; j++){
            if(i < graph[i][j]){
                ADDONEEDGE(ng, i - 1, graph[i][j] - 1, m);
            }
        }
    }
}

inline void callNauty(){
    
    //call nauty
    prepareNautyCall();
    nauty((graph*) &ng, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, vertexCount + addedVerticesCount, (graph*) &ng_canon);
    
    generatorsDetermined[addedVerticesCount] = TRUE;
}

/**
 * Adds the edge uv to the cycle. This method assumes that uv is an existing edge
 * in the graph.
 */
inline void addEdgeToCycle(int u, int v){
    setword *gu,*gv;
    
    //determine label of new vertex
    int newVertex = vertexCount + addedVerticesCount;
    addedVerticesCount++;
    
    //remove old edge
    gu = GRAPHROW(ng, u, m);
    gv = GRAPHROW(ng, v, m);
    DELELEMENT(gu,v);
    DELELEMENT(gv,u);
    
    //add two new edges
    ADDONEEDGE(ng, u, newVertex, m);
    ADDONEEDGE(ng, v, newVertex, m);
}

/**
 * Removes the edge uv from the cycle. This assumes that uv was the last edge
 * to be added to the cycle.
 */
inline void removeLastEdgeFromCycle(int u, int v){
    setword *gu,*gv, *gn;
    
    //determine label of subdividing vertex
    addedVerticesCount--;
    int newVertex = vertexCount + addedVerticesCount;
    
    //remove subdivided edges
    gu = GRAPHROW(ng, u, m);
    gv = GRAPHROW(ng, v, m);
    gn = GRAPHROW(ng, newVertex, m);
    DELELEMENT(gu,newVertex);
    DELELEMENT(gv,newVertex);
    DELELEMENT(gn,u);
    DELELEMENT(gn,v);
    
    //add original edge
    ADDONEEDGE(ng, u, v, m);
}

int findRootOfElement(int forest[], int element) {
    //find with path-compression
    if(element!=forest[element]){
        forest[element]=findRootOfElement(forest, forest[element]);
    }
    return forest[element];
}

void unionElements(int forest[], int treeSizes[], int *numberOfComponents, int element1, int element2){
    int root1 = findRootOfElement(forest, element1);
    int root2 = findRootOfElement(forest, element2);

    if(root1==root2) return;

    if(treeSizes[root1]<treeSizes[root2]){
        forest[root1]=root2;
        treeSizes[root2]+=treeSizes[root1];
    } else {
        forest[root2]=root1;
        treeSizes[root1]+=treeSizes[root2];
    }
    (*numberOfComponents)--;
}

void determineExtensionEdgeOrbits(VERTEXPAIR edges[], int edgesCount, int edgeOrbits[], int *orbitsCount) {
    int i, j, k;
    int orbitSize[edgesCount];

    //initialization of the variables
    for(i=0; i<edgesCount; i++){
        edgeOrbits[i]=i;
        orbitSize[i]=1;
    }
    *orbitsCount=edgesCount;

    if(generatorCount[addedVerticesCount]==0){
        //if the automorphism group is trivial
        return;
    }

    int *permutation;
    VERTEXPAIR edge;
    
    for(i = 0; i < generatorCount[addedVerticesCount]; i++) {
        permutation = generators[addedVerticesCount][i];

        for(j = 0; j<edgesCount; j++){
            //apply permutation to current edge
            edge[0] = permutation[edges[j][0]];
            edge[1] = permutation[edges[j][1]];

            //assert: edge[0] == endpoint1 or endpoint2

            //search the pair in the list
            for(k = 0; k<edgesCount; k++){
                if(edge[0] == edges[k][0] && edge[1] == edges[k][1]){
                    unionElements(edgeOrbits, orbitSize, orbitsCount, j, k);
                    break; //the list of edges doesn't contain any duplicates so we can stop
                }
            }
        }
    }

    //make sure that each element is connected to its root
    for(i = 0; i < edgesCount; i++){
        findRootOfElement(edgeOrbits, i);
    }
}

void determineGeneralEdgeOrbits(VERTEXPAIR edges[], int edgesCount, int edgeOrbits[], int *orbitsCount) {
    int i, j, k, temp;
    int orbitSize[edgesCount];

    //initialization of the variables
    for(i=0; i<edgesCount; i++){
        edgeOrbits[i]=i;
        orbitSize[i]=1;
    }
    *orbitsCount=edgesCount;

    if(generatorCount[addedVerticesCount]==0){
        //if the automorphism group is trivial
        return;
    }

    int *permutation;
    VERTEXPAIR edge;
    
    for(i = 0; i < generatorCount[addedVerticesCount]; i++) {
        permutation = generators[addedVerticesCount][i];

        for(j = 0; j<edgesCount; j++){
            //apply permutation to current edge
            edge[0] = permutation[edges[j][0]];
            edge[1] = permutation[edges[j][1]];

            //canonical form of edge
            if(edge[0]>edge[1]){
                temp = edge[1];
                edge[1] = edge[0];
                edge[0] = temp;
            }

            //search the pair in the list
            for(k = 0; k<edgesCount; k++){
                if(edge[0] == edges[k][0] && edge[1] == edges[k][1]){
                    unionElements(edgeOrbits, orbitSize, orbitsCount, j, k);
                    break; //the list of edges doesn't contain any duplicates so we can stop
                }
            }
        }
    }

    //make sure that each element is connected to its root
    for(i = 0; i < edgesCount; i++){
        findRootOfElement(edgeOrbits, i);
    }
}

boolean isLastAddedEdgeCanonicalCycle(){
    //at the moment we don't use colours: just call nauty and check that the
    //last edge is in the orbit of the smallest edge
    callNauty();
    
    int reverseLabelling[MAXN];
    int i;
    for (i = 0; i < vertexCount + addedVerticesCount; i++) {
        reverseLabelling[lab[i]]=i;
    }
    
    int lastEdgeVertex = vertexCount + addedVerticesCount - 1;
    
    int smallestLabelInOrbitLastEdge = reverseLabelling[lastEdgeVertex];
    int smallestEdgeLabel = MAXN;
    for(i = vertexCount; i < vertexCount + addedVerticesCount; i++){
        if(orbits[i] == orbits[lastEdgeVertex]){
            if(reverseLabelling[i] < smallestLabelInOrbitLastEdge){
                smallestLabelInOrbitLastEdge = reverseLabelling[i];
            } else if(reverseLabelling[i] < smallestEdgeLabel){
                smallestEdgeLabel = reverseLabelling[i];
            }
        }
    }
    
    if(smallestEdgeLabel < smallestLabelInOrbitLastEdge){
        return FALSE;
    } else {
        return TRUE;
    }
}

boolean isLastAddedEdgeCanonicalPath(int otherEndPoint){
    //at the moment we don't use colours: just call nauty and check that the
    //last edge is in the orbit of the smallest edge
    callNauty();
    
    int reverseLabelling[MAXN];
    int i;
    for (i = 0; i < vertexCount + addedVerticesCount; i++) {
        reverseLabelling[lab[i]]=i;
    }
    
    int lastEdgeVertex = vertexCount + addedVerticesCount - 1;
    
    if(orbits[lastEdgeVertex] == orbits[otherEndPoint]){
        //the two end edges are in the same orbit
        return TRUE;
    } else {
        //if they are not in the same orbit, then they each form an orbit by themselves
        return reverseLabelling[lastEdgeVertex] < reverseLabelling[otherEndPoint];
    }
}

void handleFinishedCycle(){

}

void closeCycle(int endpoint1, int endpoint2){
    addEdgeToCycle(endpoint1, endpoint2);
    
    //check that last edge was canonical
    if(isLastAddedEdgeCanonicalCycle()){
        handleFinishedCycle();
    }
    
    removeLastEdgeFromCycle(endpoint1, endpoint2);
}

void extendCycle(int endpoint1, int edgeVertex1, int endpoint2, int edgeVertex2){
    int i, edgeCount = 0, edgeOrbitCount;
    VERTEXPAIR edges[2*MAXN];
    int edgeOrbits[2*MAXN];
    setword *gep1,*gep2;
    
    //load the neighbourhoods of the two endpoints
    gep1 = GRAPHROW(ng, endpoint1, m);
    gep2 = GRAPHROW(ng, endpoint2, m);
    
    //make the possible extensions
    if(!generatorsDetermined[addedVerticesCount]){
        //first call Nauty
        callNauty();
    }
    
    //find all possible extensions
    for (i = -1; (i = nextelement(gep1,m,i)) >= 0;){
        if((!ISELEMENT(verticesInCycle, i)) && i != endpoint2 && i < vertexCount){
            edges[edgeCount][0] = endpoint1;
            edges[edgeCount][1] = i;
            edgeCount++;
        }
    }
    for (i = -1; (i = nextelement(gep2,m,i)) >= 0;){
        if((!ISELEMENT(verticesInCycle, i)) && i != endpoint1 && i < vertexCount){
            edges[edgeCount][0] = endpoint2;
            edges[edgeCount][1] = i;
            edgeCount++;
        }
    }
    
    //partition the extensions into orbits
    determineExtensionEdgeOrbits(edges, edgeCount, edgeOrbits, &edgeOrbitCount);
    
    //make the extensions
    for (i = 0; i < edgeCount; i++) {
        if(orbits[i]==i){
            addEdgeToCycle(edges[i][0], edges[i][1]);
            ADDELEMENT(verticesInCycle, edges[i][1]);
            
            //check that last edge was canonical
            if(edges[i][0]==endpoint1){
                if(isLastAddedEdgeCanonicalPath(edgeVertex2)){
                    //if yes: continue extending
                    extendCycle(edges[i][1], vertexCount + addedVerticesCount - 1,
                            endpoint2, edgeVertex2);
                }
            } else {
                if(isLastAddedEdgeCanonicalPath(edgeVertex1)){
                    //if yes: continue extending
                    extendCycle(endpoint1, edgeVertex1,
                            edges[i][1], vertexCount + addedVerticesCount - 1);
                }
            }
            
            DELELEMENT(verticesInCycle, edges[i][1]);
            removeLastEdgeFromCycle(edges[i][0], edges[i][1]);
        }
    }
    
    //check whether we can close this cycle
    if(ISELEMENT(gep1, endpoint2)){
        closeCycle(endpoint1, endpoint2);
    }
}

void startBuildingCycles(){
    int i, v, edgeCount = 0, edgeOrbitCount;
    VERTEXPAIR edges[MAXN * (MAXN - 1)/2];
    int edgeOrbits[MAXN * (MAXN - 1)/2];
    setword *gv;
    
    //build a list of all edges
    for(v = 0; v < vertexCount; v++){
        gv = GRAPHROW(ng, v, m);
        for (i = -1; (i = nextelement(gv,m,i)) >= 0;){
            if(v < i){
                edges[edgeCount][0] = v;
                edges[edgeCount][1] = i;
                edgeCount++;
            }
        }
    }
    
    //call Nauty so we can build the orbits of edges
    callNauty();
    
    //partition the edges into orbits
    determineGeneralEdgeOrbits(edges, edgeCount, edgeOrbits, &edgeOrbitCount);
    
    //clear set of vertices in the cycle
    EMPTYSET(verticesInCycle, m);
    
    //make the extensions
    for (i = 0; i < edgeCount; i++) {
        if(orbits[i]==i){
            addEdgeToCycle(edges[i][0], edges[i][1]);
            ADDELEMENT(verticesInCycle, edges[i][0]);
            ADDELEMENT(verticesInCycle, edges[i][1]);
            
            //a single edge is always canonical, so no canonicity check is needed
            extendCycle(edges[i][0], vertexCount, edges[i][1], vertexCount);
            
            DELELEMENT(verticesInCycle, edges[i][0]);
            DELELEMENT(verticesInCycle, edges[i][1]);
            removeLastEdgeFromCycle(edges[i][0], edges[i][1]);
        }
    }
}

boolean isGraphPancyclic(GRAPH graph, ADJACENCY adj){
    vertexCount = graph[0][0];
    
    translateGraphToNautyDenseGraph(graph, adj);
    
    startBuildingCycles();
    
    return FALSE;
}

//====================== USAGE =======================

void help(char *name) {
    fprintf(stderr, "The program %s checks whether simple graphs are pancyclic.\n\n", name);
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "\nThis program can handle graphs up to %d vertices and degrees up to %d.\n\n", MAXN, MAXVAL);
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
         {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case '?':
                usage(name);
                return EXIT_FAILURE;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }

    GRAPH graph;
    ADJACENCY adj;
    int i;
    
    unsigned short code[MAXCODELENGTH];
    int length;
    while (readMultiCode(code, &length, stdin)) {
        decodeMultiCode(code, length, graph, adj);
        
        //check that there are no vertices of degree 2
        for(i = 1; i <= graph[0][0]; i++){
            if(adj[i]==2){
                fprintf(stderr, "Vertices of degree 2 are not supported -- exiting!\n");
                return EXIT_FAILURE;
            }
        }
        
        if(isGraphPancyclic(graph, adj)){
            
        }
    }
    
    return EXIT_SUCCESS;
}

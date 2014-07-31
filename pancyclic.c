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

permutation generators[MAXN+1][MAXN/2][MAXN];
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
    gn = GRAPHROW(newVertex, v, m);
    DELELEMENT(gu,newVertex);
    DELELEMENT(gv,newVertex);
    DELELEMENT(gn,u);
    DELELEMENT(gn,v);
    
    //add original edge
    ADDONEEDGE(ng, u, v, m);
}

void closeCycle(int endpoint1, int endpoint2){
    addEdgeToCycle(endpoint1, endpoint2);
    
    //check that last edge was canonical
    //handle graph
    
    removeLastEdgeFromCycle(endpoint1, endpoint2);
}

void extendCycle(int endpoint1, int endpoint2){
    setword *gep1,*gep2;
    gep1 = GRAPHROW(ng, endpoint1, m);
    gep2 = GRAPHROW(ng, endpoint2, m);
    
    //check whether we can close this cycle
    if(ISELEMENT(gep1, endpoint2)){
        closeCycle(endpoint1, endpoint2);
    }
    
    //make the other possible extensions
    if(!generatorsDetermined[addedVerticesCount]){
        //first call Nauty
    }
    
    //find all possible extensions
    
    //partition the extensions into orbits
    
    //make the extensions
}

boolean isGraphPancyclic(GRAPH graph, ADJACENCY adj){
    vertexCount = graph[0][0];
    
    translateGraphToNautyDenseGraph(graph, adj);
    
    EMPTYSET(verticesInCycle, m);
    
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

#ifndef _DJS_DEGENERACY_ALGORITHM_CLIQUES_H_
#define _DJS_DEGENERACY_ALGORITHM_CLIQUES_H_

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"degeneracy_helper.h"

void init_MC(int n, LinkedList** adjListLinked); 

void printsizetree(int k, int n);

void listAllCliquesDegeneracyRecursive_A(  treeNode *root, int depth, 
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR, int max_k, 
                                               int SCT_k, 
                                               int rsize, int drop);

void listAllCliquesDegeneracy_A( treeNode *root, NeighborListArray** orderingArray,
                                      int size, int max_k, int SCT_k); 


void hybridList(NeighborListArray** orderingArray, int size, int max_k, int curr_k, int *corenumber,    
            double nCr[1001][401], long long unsigned int *rho, int *hold_set, int *hold_len, int *pivot_set, int *pivot_len ); 

void visitSmallTreesRecursive( int lastAddedVertex, int depth, int* vertexSets, int* vertexLookup,
                                            int** neighborsInP, int* numNeighbors,
                                            int beginX, int beginP, int beginR, 
                                            int max_k, int curr_k, int rsize, int drop, 
                                            long long unsigned int *rho, int *hold_set, int *hold_len, int *pivot_set, int *pivot_len 
                                            ); 

void hybrid_compute_yi(NeighborListArray** orderingArray, int size, int max_k, int curr_k, int *corenumber,    
        int *y, int *level, int *hold_set, int *hold_len, int *pivot_set, int *pivot_len ); 


void computeYiRecursive( int lastAddedVertex, int depth, int* vertexSets, int* vertexLookup,
                                            int** neighborsInP, int* numNeighbors,
                                            int beginX, int beginP, int beginR, 
                                            int max_k, int curr_k, int rsize, int drop, 
                                            int *y, int *level, 
                                            int *hold_set, int *hold_len, int *pivot_set, int *pivot_len 
                                            );

#endif

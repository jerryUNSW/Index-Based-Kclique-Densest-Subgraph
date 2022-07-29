#ifndef _DJS_DEGENERACY_ALGORITHM_CLIQUES_E_H_
#define _DJS_DEGENERACY_ALGORITHM_CLIQUES_E_H_


#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"degeneracy_helper.h"


void listAllCliquesDegeneracyRecursive_E(double* cliqueCounts,
                                               int *ordering,
                                               int *CSCindex,
                                               int *CSCedges,
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR, int keep, int drop, int *keepV, int *dropV, int max_k);

void listAllCliquesDegeneracy_E(double* cliqueCounts, 
                                      NeighborListArray**,
                                      int *ordering,
                                      int *CSCindex,
                                      int *CSCedges,
                                      int size, int max_k );

#endif

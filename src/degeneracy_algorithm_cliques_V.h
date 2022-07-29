#ifndef _DJS_DEGENERACY_ALGORITHM_CLIQUES_V_H_
#define _DJS_DEGENERACY_ALGORITHM_CLIQUES_V_H_

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"degeneracy_helper.h"




void listAllCliquesDegeneracyRecursive_V(double *,
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR, int keep, int drop, int *keepV, int *dropV, int max_k);

void listAllCliquesDegeneracy_V(double *, NeighborListArray**,
                                      int size, int max_k );

#endif

#ifndef _DJS_MISC_H_
#define _DJS_MISC_H_

#include"LinkedList.h"
#include"degeneracy_helper.h"
#include"degeneracy_algorithm_cliques_A.h"

#define max(x,y) (x > y? x:y)
#define min(x,y) (x < y? x:y)
#define MAX_CSIZE 400

static int ck_compare(const void *a, const void *b) ;

bool compareTreeNode(treeNode* x, treeNode* y);

void populate_nCr();

int nodeComparator(void* node1, void* node2);

void printArray(int* array, int size);

void printArrayOfLinkedLists(LinkedList** listOfLists, int size);

void printInt(void* integer);

void load_balance_1(); 
void load_balance_laying_bricks(); 

void destroyCliqueResults(LinkedList* cliques);

LinkedList** readInGraphAdjList(int* n, int* m);

LinkedList** readInGraphAdjListToDoubleEdges(int* n, int* m, char *fpath);


void runAndPrintStatsCliques(LinkedList** adjListLinked,
                               int n, const char * gname, 
                               char T, int max_k, int flag_d);


int findNbrCSC(int u, int v, int *CSCindex, int *CSCedges);

void moveFromRToXDegeneracyCliques( int vertex, 
                                    int* vertexSets, int* vertexLookup, 
                                    int* pBeginX, int* pBeginP, int* pBeginR);

void moveToRDegeneracyCliques( int vertex, 
                               int* vertexSets, int* vertexLookup, 
                               int** neighborsInP, int* numNeighbors,
                               int* pBeginX, int *pBeginP, int *pBeginR, 
                               int* pNewBeginX, int* pNewBeginP, int *pNewBeginR);

void fillInPandXForRecursiveCallDegeneracyCliques( int vertex, int orderNumber,
                                                   int* vertexSets, int* vertexLookup, 
                                                   NeighborListArray** orderingArray,
                                                   int** neighborsInP, int* numNeighbors,
                                                   int* pBeginX, int *pBeginP, int *pBeginR, 
                                                   int* pNewBeginX, int* pNewBeginP, int *pNewBeginR);

int findBestPivotNonNeighborsDegeneracyCliques( int** pivotNonNeighbors, int* numNonNeighbors,
                                                int* vertexSets, int* vertexLookup,
                                                int** neighborsInP, int* numNeighbors,
                                                int beginX, int beginP, int beginR);

#endif


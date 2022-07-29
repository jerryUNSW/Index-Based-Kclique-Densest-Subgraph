#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"

int main()
{

    int n; // number of vertices
    int m; // 2x number of edges

    LinkedList** adjacencyList = readInGraphAdjList(&n,&m);

    int d = computeDegeneracy(adjacencyList, n);

    int i = 0;

    while(i<n)
    {
        destroyLinkedList(adjacencyList[i]);
        i++;
    }

    Free(adjacencyList); 

    fprintf(stderr, "Degeneracy is %d\n", d);
    return 0;
}

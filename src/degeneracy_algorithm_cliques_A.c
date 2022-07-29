#include<limits.h>
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<vector>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"degeneracy_helper.h"
#include"degeneracy_algorithm_cliques_A.h"

using namespace std; 

std::vector<treeNode*> Invalid_Branches; 

std::vector<int> tree_possible_vertices; 

bool CN_PRUNE = true ;      

bool OUT_DEG_PRUNE = true ; 

bool SHUFFLE = false;

bool COUNT = true; 

bool KCC_COUNT = true ; 

unsigned long long int hybrid_updates = 0 ;

unsigned long long int size_of_tree = 1 ; 

unsigned long long int size_subtree = 0 ;

extern double nCr[1001][401];

extern vector<int> hset; 
extern vector<int> pset; 
extern vector<int> corenumber;

vector<bool> MC_lookup; 
vector<int> curr_mc; 
vector<int> max_k_per_ver ; 
vector<double> cliq_cnt; 
vector<int> stable_set ; 

vector<vector<double>> per_ver_cnt ; 

// in this function, we initialize the bool look up table
void init_MC(int n, LinkedList** adjListLinked)
{
    // printf("size of G' = %d \n", curr_mc.size());

    // fill(max_cliq_per_ver.begin(),max_cliq_per_ver.end(), 0); 
    MC_lookup.resize(n);
    fill(MC_lookup.begin(),MC_lookup.end(), false); 
    // for(auto v : curr_mc){ MC_lookup[v] = true ; }
}

// max_k here is just deg+1. SCT_k is the minK in the index
// instead we can declare corenumber as a global variable (vector)
void listAllCliquesDegeneracy_A(treeNode *root, NeighborListArray** orderingArray, 
                                      int size, int max_k, int SCT_k)
{
    
    // revisiting the memery consumption here!!! 
    if(COUNT){
        cliq_cnt.resize(max_k+1); 
        printf("max_k = %d \n", max_k);
        fill(cliq_cnt.begin(), cliq_cnt.end(), 0 ); 

        // KCC related variables. we do not need to store them. just compute them when needed. 
        if(KCC_COUNT){
            max_k_per_ver.resize(size);
            fill(max_k_per_ver.begin(), max_k_per_ver.end(),0);

            // when do we need to count this? only need it for KCC analysis. 
            per_ver_cnt.resize(size);
            for(int i=0;i<size;i++)
            {
                // per_ver_cnt[i].resize (orderingArray[i]->laterDegree+2) ;  // this is causing problems 
                per_ver_cnt[i].resize (max_k+1 ) ; // maybe this can be shrinked.
                assert(orderingArray[i]->laterDegree < max_k);
                fill(per_ver_cnt[i].begin(),per_ver_cnt[i].end(), 0); 
            } 
        }
    }

    // vertex sets are stored in an array like this:
    // |--X--|--P--|
    int* vertexSets = (int *)Calloc(size, sizeof(int));

    // vertex i is stored in vertexSets[vertexLookup[i]]
    int* vertexLookup = (int *)Calloc(size, sizeof(int));

    int** neighborsInP = (int **)Calloc(size, sizeof(int*));
    int* numNeighbors = (int *)Calloc(size, sizeof(int));
   
    int i = 0;

    while(i<size)
    {
        vertexLookup[i] = i;
        vertexSets[i] = i;
        neighborsInP[i] = (int *)Calloc(1, sizeof(int));
        numNeighbors[i] = 1;
        i++;
    }

    int beginX = 0;
    int beginP = 0;
    int beginR = size;

    // for each vertex
    int NN = 100; 
    bool visited[NN];
    for(int i=0;i<NN;i++) 
        visited[i]=false;

    for(i=0;i<size;i++)
    {
        // if(orderingArray[i]==NULL) continue ; 
        size_subtree = 0 ;

        int vertex = (int)orderingArray[i]->vertex;

        if(CN_PRUNE && corenumber[vertex]+1 < SCT_k)
            continue;
        if(OUT_DEG_PRUNE && orderingArray[i]->laterDegree +1 < SCT_k)
            continue; 
         
        treeNode *currtree = (treeNode*) Calloc(1, sizeof(treeNode));
        currtree->id = vertex ; 
        currtree->label = 0 ; // initialize as a hold node.
        hset.push_back(vertex);

        size_of_tree++;
        size_subtree++;
        // printf("\nprocessing vertex %d \n", vertex );        
        int newBeginX, newBeginP, newBeginR;

        fillInPandXForRecursiveCallDegeneracyCliques( i, vertex, 
                                               vertexSets, vertexLookup, 
                                               orderingArray,
                                               neighborsInP, numNeighbors,
                                               &beginX, &beginP, &beginR, 
                                               &newBeginX, &newBeginP, &newBeginR);


        // recursively compute maximal cliques containing vertex, some of its later neighbors, and avoiding earlier neighbors
        int drop = 0;
        int rsize = 1;
        int depth = 1 ; 

        currtree->max_depth = 1 ; 

        // per_ver = 0 ;

        listAllCliquesDegeneracyRecursive_A(currtree, depth, 
                                                  vertexSets, vertexLookup,
                                                  neighborsInP, numNeighbors,
                                                  newBeginX, newBeginP, newBeginR, 
                                                  max_k, SCT_k, rsize, drop); 

        root->max_depth = root->max_depth > currtree->max_depth ?  root->max_depth : currtree->max_depth; 

        root->child_node.push_back(currtree);

        hset.pop_back(); // or hset.clear()
       
        beginR = beginR + 1;

    }

    for(auto v : curr_mc) MC_lookup[v] = true ; 

    // in the end show the size of the tree.
    printsizetree(SCT_k, size);

    // printf("Now need to free some data structured not to be used. \n"); 
    Free(vertexSets);
    Free(vertexLookup);

    for(i = 0; i<size; i++)
    {
        Free(neighborsInP[i]);
    }

    Free(neighborsInP);
    Free(numNeighbors);
    // printf("number of empty slots = %d \n", tree_possible_vertices.size());

    // printf("num of maximals = %llu \n", maximal_cliq_cnt);
    return;
}

void printsizetree(int k, int n)
{

    printf("currmc = %d \n",  curr_mc.size());

    // for(auto v : curr_mc) MC_lookup[v] = true ; 

    printf("current tree size = %ld, k = %d, max k = %d , n = %d \n", size_of_tree, k, cliq_cnt.size(), n ); 

    for(int i=3;i<cliq_cnt.size();i++){

        if(cliq_cnt[i]>0) {
            // printf("k=%d, clique count = %0.f\n", i , cliq_cnt[i]);

            for(int ver = 0; ver<n; ver++){
                
                assert(per_ver_cnt[ver][i]>=0);
            }
        }

    }
    
}

void listAllCliquesDegeneracyRecursive_A( treeNode *root, int depth, 
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR, int max_k, 
                                               int SCT_k, 
                                               int rsize, int drop)
{   
    // out-degree based pruning
    // only do this when not hybrid 
    if(OUT_DEG_PRUNE && depth +beginR - beginP < SCT_k ){
        if(COUNT){
            for(int i=drop; (i>=0)&&(rsize-i<=max_k);i--) cliq_cnt[rsize-i]+=nCr[drop][i];
        }        
        return;
    }
        // return ; 
    // we should also add up the kclique numbers in these cases. 
    // for small k-cliques. 

    // at leave node, need to recurse back 
    if ((beginP >= beginR) || (rsize-drop > max_k)) 
    {
        if(COUNT){

            // compute the total number of kcliques. 
            for(int i=drop; (i>=0)&&(rsize-i<=max_k);i--) {
                cliq_cnt[rsize-i]+=nCr[drop][i];
            }

            if(KCC_COUNT){
                // compute the number of kcliques each vertex is contained. 
                for(int i = 0; i<= pset.size(); i++){
                    for(auto v : hset){
                        per_ver_cnt[v][hset.size()+i] += nCr[pset.size()][i]; 
                        max_k_per_ver[v] = max_k_per_ver[v] > hset.size()+i ? max_k_per_ver[v] : hset.size()+i ; 
                    }
                }
                for(int i = 0; i<= pset.size()-1; i++){
                    for(auto v : pset){
                        per_ver_cnt[v][hset.size()+i+1] += nCr[pset.size()-1][i]; 
                        max_k_per_ver[v] = max_k_per_ver[v] > hset.size()+i+1 ? max_k_per_ver[v] : hset.size()+i+1 ; 
                    }
                }
            }
            // report the current maximum clique size: 
            if(curr_mc.size() < rsize){
                // update curr_mc
                curr_mc.clear();
                curr_mc = hset;
                curr_mc.insert(curr_mc.end(), pset.begin(), pset.end());
            }
        }
        return;
    }
    
    depth++; 

    int* myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough = 0;

    // get the candidates to add to R to make a maximal clique
    int pivot = findBestPivotNonNeighborsDegeneracyCliques( &myCandidatesToIterateThrough,
                                         &numCandidatesToIterateThrough,
                                         vertexSets, vertexLookup,
                                         neighborsInP, numNeighbors,
                                         beginX, beginP, beginR);
    // printf("pivot chosen = %d  \n", pivot);
    // exit(0);
    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    if(numCandidatesToIterateThrough != 0)
    {
        int iterator = 0;
        
        // expandChildren(root, numCandidatesToIterateThrough) ;
        // printf("root = %d, numCandidatesToIterateThrough:%d\n", root->id, numCandidatesToIterateThrough); 

        while(iterator < numCandidatesToIterateThrough)
        {
            // vertex to be added to the partial clique
            int vertex = myCandidatesToIterateThrough[iterator];

            // add a new node to root. 
            treeNode *currtree = (treeNode*) Calloc(1, sizeof(treeNode));
            currtree->id = vertex ; 
            currtree->max_depth = depth ;
            size_of_tree++;
            size_subtree++;

            int newBeginX, newBeginP, newBeginR;

            // add vertex into partialClique, representing R.
            // printf("add vertex into partialClique, representing R.\n");

            // swap vertex into R and update all data structures 
            // printf("swap vertex into R and update all data structures \n");
            moveToRDegeneracyCliques( vertex, 
                               vertexSets, vertexLookup, 
                               neighborsInP, numNeighbors,
                               &beginX, &beginP, &beginR, 
                               &newBeginX, &newBeginP, &newBeginR);

            // recursively compute maximal cliques with new sets R, P and X
            if (vertex == pivot){
                currtree->label = 1 ; 
                pset.push_back(vertex);
                listAllCliquesDegeneracyRecursive_A(currtree, depth, 
                                                      vertexSets, vertexLookup,
                                                      neighborsInP, numNeighbors,
                                                      newBeginX, newBeginP, newBeginR, max_k, SCT_k, rsize+1, drop+1);
                pset.pop_back();
            }
            else{
                currtree->label = 0 ; 
                hset.push_back(vertex);
                listAllCliquesDegeneracyRecursive_A(currtree, depth, 
                                                      vertexSets, vertexLookup,
                                                      neighborsInP, numNeighbors,
                                                      newBeginX, newBeginP, newBeginR, max_k, SCT_k, rsize+1, drop);
                hset.pop_back();
            }

            root->max_depth = root->max_depth > currtree->max_depth ?  root->max_depth : currtree->max_depth; 
            root->child_node.push_back(currtree);

            moveFromRToXDegeneracyCliques( vertex, 
                                    vertexSets, vertexLookup,
                                    &beginX, &beginP, &beginR );

            iterator++;
        }

        // swap vertices that were moved to X back into P, for higher recursive calls.
        iterator = 0;
        while(iterator < numCandidatesToIterateThrough)
        {
            int vertex = myCandidatesToIterateThrough[iterator];
            int vertexLocation = vertexLookup[vertex];
            beginP--;
            vertexSets[vertexLocation] = vertexSets[beginP];
            vertexSets[beginP] = vertex;
            vertexLookup[vertex] = beginP;
            vertexLookup[vertexSets[vertexLocation]] = vertexLocation;
            iterator++;
        }
    }
    Free(myCandidatesToIterateThrough);
    return;
}
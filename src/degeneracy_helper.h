#ifndef _DJS_DEGENERACY_HELPER_H_
#define _DJS_DEGENERACY_HELPER_H_

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include <unistd.h>

#include <chrono>
#include <random>
#include<algorithm>
#include<vector>
#include <queue>
#include <math.h>
#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"

// added from maxflow code
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
using namespace std;
using namespace boost;


class Network {
private:
  typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
public:
  typedef Traits::vertex_descriptor Vertex;
  Network();
  Vertex AddVertex();
  void AddEdge(Vertex &v1, Vertex &v2, const long capacity);
  unsigned long long MaxFlow(Vertex &s, Vertex &t);
private:
  typedef adjacency_list < vecS, vecS, directedS,
    property < vertex_name_t, std::string,
    property < vertex_index_t, long,
    property < vertex_color_t, boost::default_color_type,
    property < vertex_distance_t, long,
    property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,

    property < edge_capacity_t, long long,
    property < edge_residual_capacity_t, long long,
    property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;
  Graph g;
  property_map < Graph, edge_reverse_t >::type rev;
};

struct treeNode
{
    int id ;    // -1 for root, range is between -1 and max vertex id. 
    int label  ; // 1 for pivot type, 0 for hold type, always initialize to 0. when pivot, update to 1 .
    
    vector<treeNode*> child_node; 
    treeNode* parent;
    int pos; // curr->parent->child_nodes[pos] == curr
    
    int max_depth;         // max depth is just 51 on Orkut
    long int num_cliques ; // the number of k-cliques under this subtree
};

vector<int> get_curr_mc(); 

void pt(int xx);

bool maxflow_check_optimal(treeNode *root) ; 

bool goldberg_check_optimal( double density, int num_iter_run); 

void assign_weight_MC(long long unsigned int *rho, int k); 
void listYi_MC(int r, int *level, int *y);
void MC_list_yi(int *y, long long unsigned int *rho); 
void MC_summary();
void combineYi_MC(int *data, int start, int end, int index, int r, int *level, int *y); 

void change_mc_lookup(vector<bool> input); 



// exact algorithm 
void collect_cliques(int kq, treeNode *root, int n );


int find_func(int i);

void union_func();

void UFKclique(treeNode *root, int depth, int max_k);

void init_work_by_thread();

void showemptybranches();

void shuffle(int *array, size_t n);

int comparitor (const void* a, const void* b); 

void update_length( treeNode *root); 


// sampling based 
int get_sampled_cliques(int remainder, int r); 
int get_sampled_cliques_helper(int remainder, int *visited, int *data, int start, int end, int index, int r); 


// maxflow related code
void list_enumeration_cliques_in_approximate_solution(int r);
void combine_enumeration_cliqs_N_approx(int *data, int start, int end, int index, int r);
void enumerate_cliques_in_approximate_solution(treeNode *root, int num_nodes, int max_k, int *level);

// compute approximate densest subgraph.
void listForYi(int r,int *level, unsigned long long int *y);
void combineForYi(int *data, int start, int end, int index, int r, int *level, unsigned long long int *y, int *argmaxhold, int *maxhold);

void printNumUpdates(double *timer__, double seqtime, int max_k);

void show_arrays(int tmp, long long unsigned int *rho, int *hold_set, int *hold_len, int *pivot_set, int *pivot_len);

void even_distribute_weight(double updates, int count1, long long unsigned int *rho, int *argmin_h_idxs);

// kclist++
void updateRho(int curr_itr, long long unsigned int *rho, int r); 

void updateRhoHelper(int curr_itr, int *visited, long long unsigned int *rho, int *data, int start, int end, int index, int r, int *argminhold, int *minhold);


int combinationUtil_exact(int *data, int start, int end, int index, int r) ; 

int printCombination_exact(int r);

// kclique lisitng
int printCombination(int r); 
int combinationUtil(int *data, int start, int end, int index, int r);

void kcliquePlusPlus(int n, long long unsigned int *rho, treeNode *root, int depth, int max_k);


// void batch_rho_update(int depth, int remainder, double nCr[1001][401], long long unsigned int *rho, int *hold_set, int *hold_len, int *pivot_set, int *pivot_len, int max_k);
void batch_rho_update(int depth, int remainder, double nCr[1001][401], long long unsigned int *rho, int max_k); 
        
// void compute_yi(treeNode *root, int depth, int *hold_set, int *hold_len, int *pivot_set, int *pivot_len, int max_k, int *level, int *y);
void compute_yi(treeNode *root, int depth, int max_k, int *level, unsigned long long int *y);

void destroyTreeFromRoot(treeNode *root); 

long int ListRootToLeavePath(treeNode *root, int depth, int max_k, int base_k);

void showTree(treeNode *root, int depth, int max_k, long long unsigned int *rho, bool *lookup);

/*! \struct NeighborList

    \brief For a given ordering, this stores later neighbors and earlier neighbors
           for a given vertex in linked lists.
*/

struct NeighborList
{
    int vertex; //!< the vertex that owns this neighbor list
    LinkedList* earlier; //!< a linked list of neighbors that come before this vertex in the ordering
    LinkedList* later; //!< a linked list of neighbors that come after this vertex in the ordering
    int orderNumber; //!< the position of this verex in the ordering
};

typedef struct NeighborList NeighborList;

/*! \struct NeighborListArray

    \brief For a given ordering, this stores later neighbors and earlier neighbors
           for a given vertex in arrays.

    This version of the NeighborList structure is more cache efficient.
*/

struct NeighborListArray
{
    int vertex; //!< the vertex that owns this neighbor list
    int* earlier; //!< an array of neighbors that come before this vertex in an ordering
    int earlierDegree; //!< the number of neighbors in earlier
    int* later; //!< an array of neighbors that come after this vertex in an ordering
    int laterDegree; //!< an array of neighbors that come after this vertex in an ordering
    int orderNumber; //!< the position of this verex in the ordering
    // it is useful to include the core number here
};
typedef struct NeighborListArray NeighborListArray;

int computeDegeneracy(LinkedList** list, int size);

NeighborList** computeDegeneracyOrderList(LinkedList** list, int size);

NeighborListArray** computeDegeneracyOrderArray(LinkedList** list, int size, int max_k);

int neighborListComparator(int* nl1, int* nl2);

#endif

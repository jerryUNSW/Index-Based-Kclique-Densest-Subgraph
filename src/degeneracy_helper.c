#include<limits.h>
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"degeneracy_helper.h"


// added from maxflow code
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/graph_utility.hpp>
// #include "MaxFlow.hpp"

using namespace std;
using namespace boost;

extern double nCr[1001][401];
extern int iter_num; // this is the current iteration number starting from 0. 

extern vector<int> corenumber;

extern vector<bool> MC_lookup; 
extern vector<int> curr_mc; 
extern vector<double> cliq_cnt; 
extern vector<vector<double>> per_ver_cnt ; 
extern vector<double> cliq_cnt; 
extern vector<int> stable_set ; 

// max flow based variables
extern Network *curr_network ; 
extern long int m__ ;
extern int num_nodes ;
extern int *level ; 
extern Network::Vertex *source_ptr; 
extern vector<Network::Vertex> R;
extern double last_checked_density; 

//
extern vector<double> Ck ; 
extern vector<double> Ck_next ; 

// union-find-related data structures
vector<int> group; 
vector<int> uf_rank;

// KCC related:
vector<bool> KCC_lookup; 
 
vector<vector<int>> all_mc; 

extern bool KCC;

// when sampling, these switches need to be turned off.
bool GRAPH_REDUCE = false ; 
bool BATCH_UPDATE = true ;
bool SAMPLE = false ; 
// the switch for sampling approximation algorithm 
extern bool EXACT  ; 
bool RAMDOM = true ; 

extern int sample_size; 
unsigned long long int cliq_idx =0 ; 
extern vector<vector<int>> sampled_kcliques; 
extern vector<bool> ver_exists_in_approx; 


bool CUT_SCT = true ; 
bool MAX_DEPTH = true ;
int RECURSION_LIMIT = 50;

// essentail containers:
extern double sub_den; 
vector<int> hset; 
vector<int> pset; 

// for sampling
int num_to_enumerate = -1 ;
double sampling_ratio = 0.01; 

// for graph reduction
vector<int> pset__; 

int num_of_empty_branches = 0 ;
int size_of_pruned_tree=0; 

long int total_num_updates = 0 ;
int tmp = 0 ;

double yyy ; // number of normal updates

double time1=0, time2=0, time3=0; 
double sort_time = 0 ; 

long int work_by_thread[32]; 

// exact algorithm related: 
vector<double> alpha ; 
vector<int> ck_exact;
vector<double> rho_exact ; 

// sampling algorithm related:
vector<int> ck_sampling;

extern int *reordered; 
extern long long unsigned int *rho ; 
extern int cliquesize_ ; 
extern int graph_num_vertices;

// this function has issues. 
// On email when k =5, it returns a solution of 265 vertices. 
bool maxflow_check_optimal(treeNode *root){

    printf("Check optimal by maxflow...\n");
    Network network;
    curr_network = &network; 
    
    // unsigned *id_in_network = (unsigned *)malloc(n* sizeof(unsigned));
    m__ = 0;
    R.clear();

    Network::Vertex s = network.AddVertex(), t = network.AddVertex();

    source_ptr = &s;  

    printf("num_nodes = %d , cliquesize_ = %d\n", num_nodes, cliquesize_);

    for (unsigned i_ = 0; i_ < num_nodes; ++i_) 
        R.push_back(network.AddVertex());

    // for each clique contained in the approximate solution.
    enumerate_cliques_in_approximate_solution(root, num_nodes, cliquesize_, level);

    for (unsigned i = 0; i < num_nodes; ++i)
        network.AddEdge(R[i], t, m__);

    return network.MaxFlow(s, t) >= m__ * num_nodes; 
}

// this function is not costly.
bool goldberg_check_optimal( double density, int num_iter_run){
    
    printf("num_iter_run = %d \n", num_iter_run);

    bool is_optimal = true;
    bool skip = true;
    double sum_rho = 0, sum_ck = 0 ;
    double jck = 0; // j choose k
    for (unsigned j = 1; j < num_nodes; ++j) 
    {
        sum_rho += (double) rho[reordered[j - 1]]*1.0/num_iter_run;
        if (skip) {
            if (j == cliquesize_)
                jck = 1;
            else if (j > cliquesize_)
                jck = (jck * j) / (j - cliquesize_);

            if (jck / j > density )
                skip = false;
                // fprintf(stderr, "Jump to j = %u\n", j);
            else
                continue; // in these cases, LHS <= 0. 
        }
        double ub ; 
        ub = sum_rho/j;
        if ( ub - density >= 1.0 / num_nodes / j && ub - density >= (ceil(density * j) - density * j ) / j) {
            // printf("ub = %lf, j = %u\n", ub, j, sum_rho); 
            is_optimal = false;
            break;
        }
    }
    return is_optimal;
}


static int exact_cdf_NodeCmp(const void *a, const void *b) 
{
  double x = rho_exact[*(const unsigned *)a];
  double y = rho_exact[*(const unsigned *)b];
  if (x > y) return -1;
  if (x < y) return 1;
  return 0;
}

vector<int> get_curr_mc(){
    // for(auto xxx:curr_mc){
    //     printf("mc : %d \n", xxx); 
    // }
    return curr_mc;
}

void change_mc_lookup(vector<bool> input){
    // MC_lookup = input;
    // int xxxxx = 0; 
    // for(int i=0;i<MC_lookup.size();i++){
    //     if(MC_lookup[i]==true) xxxxx++;
    // }
    
    // printf("xxxx = %d \n", xxxxx);
}

void collect_cliques(int kq, treeNode *root, int n){

    unsigned long long int cnt_clique = cliq_cnt[kq];

    printf("Collect all %d-cliques, count = %llu\n", kq, cnt_clique);

    alpha.resize(kq*cnt_clique);
    fill(alpha.begin(), alpha.end(), 0); 

    rho_exact.resize(n);
    fill(rho_exact.begin(), rho_exact.end(), 0 );

    ck_exact.resize(kq*cnt_clique);

    unsigned long int cliq_tmp = ListRootToLeavePath(root, 0, kq, 0); 

    assert(cliq_tmp == cnt_clique);

    printf("Shuffle %u k-cliques\n", cnt_clique);
    srand(time(NULL));
    for (unsigned i = 1; i < cnt_clique; ++i) {
        int rand_index = rand() % (i + 1);
        for(int idx = 0; idx < kq; idx++){
            int temp = ck_exact[i*kq+idx];
            ck_exact[i*kq+idx] = ck_exact[rand_index*kq+idx];
            ck_exact[rand_index*kq+idx] = temp;
        }
    }

  for (unsigned num_iter = 1; ; num_iter <<= 1) {

    printf("Step 1: run the Frank-Wolfe based algorithm for %u rounds\n", num_iter);
    for (unsigned t = num_iter / 2 + 1; t <= num_iter; ++t) {
        if (t % 10 == 0) {
            fprintf(stderr, "Run round %u...\n", t);
        }
        // go through each clique. 
        for(unsigned ic = 0 ; ic<cnt_clique; ic++){

            unsigned node_index = 0;
            for (unsigned i = 1; i < kq; ++i) {
                if (rho_exact[ck_exact[ic*kq+node_index]] > rho_exact[ck_exact[ic*kq+i]]){
                    node_index = i;
                }
            }
            rho_exact[ ck_exact[ic*kq+node_index] ] += 1.0;
            alpha[ic*kq+node_index] += 1.0;
        }
    }

    printf("Step 2: give a tentative decomposition\n");
    int *reordered_exact = (int*) Calloc(n, sizeof(int));
    
    for (int i = 0; i < n; ++i) {
        reordered_exact[i]  = i; 
    }

    vector<int> level__, y__; 
    level__.resize(n); 
    y__.resize(n);

    // before sorting make a copy of rho exact
    vector<double> rho_cp = rho_exact ; 
    qsort(reordered_exact, n, sizeof(int), exact_cdf_NodeCmp);

    for (int i = 0; i < n; ++i) {
        level__[reordered_exact[i]] = i;
        y__[i]=0;
    }

    // go through each clique again to compute y. 
    for(unsigned ic = 0 ; ic<cnt_clique; ic++){
        unsigned node_index = 0;
        for (unsigned i = 1; i < kq; ++i) {
            if (level__[ck_exact[ic*kq+node_index]] < level__[ck_exact[ic*kq+i]]){
                node_index = i;
            }
        }
        y__[ ck_exact[ic*kq+node_index] ] += 1.0;
    }
    // compute density
    unsigned long long int sum = 0;
    int num_nodes = 1; 
    double density = -1; 
    for (int i = 0; i < n; ++i) 
    {
        sum += y__[reordered_exact[i]];
        if ((double)sum / (i + 1) >= density) { // should we use >= or >
            density = (double)sum / (i + 1);
            num_nodes = i + 1; 
        }
    }
    fprintf(stderr, "Approximate densest subgraph: %u nodes, density = %lf.\n", num_nodes, density);

    // after sorting. recover the old order.
    rho_exact = rho_cp; 

    // Step 3: Check stability and optimality
    /*
    if (CDF_CheckStability(g, partition->nag[0], k)) {
      fprintf(stderr, "The potential densest set is stable!\n");
      if (CDF_CheckDensestGoldberg(g, partition->nag[0], k)) {
        fprintf(stderr, "The first %u nodes forms a densest subgraph by criteria A!\n", partition->nag[0]);
        fprintf(ofp, "[Number of Iterations]\t[Stopping Condition]\t[Number of Nodes]\t[Number of Edges]\t[k-Clique Density]\t[Number of Max-Flow Calls]\t[Time (seconds)]\n");
        fprintf(ofp, "%u\tGoldberg\t%u\t%u\t%.12f\t%u\t", num_iter, partition->nag[0], CountEdges(g, partition->nag[0], reordered), partition->val[0], cnt_max_flow);
        break;
      }
      else if (CDF_CheckDensestMaxFlow(g, partition->nag[0], k, &cnt_max_flow)) {
        fprintf(stderr, "The first %u nodes forms a densest subgraph by criteria B!\n", partition->nag[0]);
        fprintf(ofp, "[Number of Iterations]\t[Stopping Condition]\t[Number of Nodes]\t[Number of Edges]\t[k-Clique Density]\t[Number of Max-Flow Calls]\t[Time (seconds)]\n");
        fprintf(ofp, "%u\tMax Flow\t%u\t%u\t%.12f\t%u\t", num_iter, partition->nag[0], CountEdges(g, partition->nag[0], reordered), partition->val[0], cnt_max_flow);
        break;
      }
      else {
        fprintf(stderr, "Cannot guarantee it is densest by either criteria A or criteria B.\n");
        fprintf(ofp, "%u\t%u\t%.12f\t%ld\tSTABLE BUT NOT DENSEST\n", num_iter, partition->nag[0], partition->val[0], time(NULL) - t0);
      }
    } else {
      fprintf(stderr, "The potential densest subset is not stable!\n");
      fprintf(ofp, "%u\t%u\t%.12f\t%ld\tNOT STABLE\n", num_iter, partition->nag[0], partition->val[0], time(NULL) - t0);
    }
    printf("\n");
    */

  }


}

void MC_list_yi(int *y, long long unsigned int *rho){
    // I think the following is problematic
    for(auto v : curr_mc){
        y[v] = rho[v]; 
    }
}

void assign_weight_MC(long long unsigned int *rho, int k)
{
    // if(!MC_REDUCED) return ; 

    int M = curr_mc.size();
    for(auto vertex : curr_mc){
        rho[vertex] = nCr[M][k] / M ; 
        printf("ver = %d, rho = %lf \n", vertex, rho[vertex]);
    }
}

inline int LargeRand() {
  if (RAND_MAX == 0x7fff)
    return (rand() << 15) | rand();
  return rand();
}

inline int GetRandMax() {
  if (RAND_MAX == 0x7fff)
    return 0x3fffffff;
  return RAND_MAX;
}

// this function computes the workload undertaken by each thread.
void init_work_by_thread()
{
    total_num_updates=0;
}
void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}
void printNumUpdates(double *timer__, double seqtime, int max_k)
{
    // printf("# visited k-cliques = %d , n = %d,  kdensity = %lf \n", total_num_updates, 
    //     curr_mc.size(), (double)total_num_updates/curr_mc.size() );
    printf("# visited k-cliques = %d , ratio to Ck(G)= %lf \n", 
        total_num_updates, (double)total_num_updates/cliq_cnt[max_k]);

    // int num_threads = omp_get_max_threads();

    // double sum = 0;
    // double max = -1.0*INT_MAX ; 
    // double sum_time = 0 ; 
    // double max_time = -1.0*INT_MAX ; 
    // for(int i=0;i<num_threads;i++)
    // {
    //     // printf("\tthread %d, work = %lld, time = %lf\n", i, work_by_thread[i], timer__[i]);
    //     sum+=work_by_thread[i];
    //     sum_time += timer__[i]; 
    //     max = max > work_by_thread[i] ? max : work_by_thread[i] ; 
    //     max_time = max_time > timer__[i] ? max_time : timer__[i] ; 
    // }
    // double mean = sum / num_threads ; 
    // double meant = sum_time / num_threads ; 

    // // the standard deviation is not a useful indicator 
    // // printf("work = %lld, mean = %lf, max/mean = %lf \n", (long int)sum, mean, max/mean); 
    // printf("time = %lf,  mean = %lf, max/mean = %lf \n", sum_time, meant, max_time/meant); 

    // // if(num_threads>1) printf("over head = %lf , over head span = %lf\n", sum_time - seqtime, (sum_time - seqtime)/num_threads);

    // if(num_threads>1)  printf("speedup = %lf , numthreads = %d \n", seqtime / max_time, num_threads );
    // if(num_threads==1) printf("sequential time = %lf \n", seqtime  );

    // it actually makes more sense to compute max/mean.
}

// when finished the work at leave node, effectively remove the last added vertex. 
// recurse back to previous state. 
void update_length( treeNode *root)
{
    if(root->id>=0 ){

        if(root->label==0){
            hset.pop_back();
        }else{
            pset.pop_back();
        }

    }
}

// for each clique not contained in the MC, distribute its weight evenly to non-MC vertices. 

void updateRho(int curr_itr, long long unsigned int *rho, int r)
{
    int data[r];
    int *visited = (int *) Calloc(1, sizeof(int));
    *visited = 0;
    int argminhold = -1; 
    int minhold = INT_MAX*1.0;
    for(auto ver:hset){
        if(rho[ver]<minhold){
            argminhold = ver; 
            minhold = rho[ver]; 
        }
    }
    if(GRAPH_REDUCE){
        updateRhoHelper(curr_itr, visited, rho, data, 0, pset__.size()-1, 0, r, &argminhold, &minhold);
    }else{
        updateRhoHelper(curr_itr, visited, rho, data, 0, pset.size()-1,   0, r, &argminhold, &minhold);
    }
    Free(visited);
}

void updateRhoHelper(int curr_itr, int *visited, long long unsigned int *rho, 
        int *data, int start, int end, int index, int r, int *argminhold, int *minhold)
{
    if(curr_itr!=-1 && curr_itr==*visited) 
        return ; 
    if (index == r)
    {
        int arg_min = *argminhold ; 
        double min_rho = *minhold ;
        // printf("cliq: \n");
        // for(auto v:hset) printf("h: %d\n", v);
        // for (int j=0; j<r; j++)printf("p: %d\n", data[j]);
        // printf("\n");
        // search the pivot vertices   
        for (int j=0; j<r; j++){
            if(rho[data[j]]<min_rho){
                arg_min = data[j]; 
                min_rho = rho[data[j]]; 
            }
        }
        rho[arg_min]++;
        total_num_updates++;

        if(curr_itr!=-1) // only track this number when curr_itr is a positive number
            (*visited)++; 

        return;
    }
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        if(GRAPH_REDUCE){
            data[index] = pset__[i];
        }else{
            data[index] = pset[i];
        }
        updateRhoHelper(curr_itr, visited, rho, data, i+1, end, index+1, r, argminhold, minhold);
    }
    return;
}

void kcliquePlusPlus(int n, long long unsigned int *rho, treeNode *root, int depth, int max_k)
{

    if(root==NULL) 
        return ;

    if(root->id>=0){

        // KCC-based graph reduction.
        if(KCC && (!KCC_lookup[find_func(root->id)])) 
            return;  

        if(root->label==0){

            // if(hset.size()+1>max_k)return;// we do not need too many hold vertices. 

            if(GRAPH_REDUCE && Ck[root->id]< ceil(sub_den)) // if one of the hold vertex should be discarded, the whole branch is discarded
                return ; 

            hset.push_back(root->id); 
        }else{
            pset.push_back(root->id); 
        }
    }
    // recurse on its children if any
    if(root->child_node.size()>0){
        // we only shuffle the child_node vector when depth > 0 .
        // if(RAMDOM && depth>0){
        if(RAMDOM ){
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::shuffle(std::begin(root->child_node), std::end(root->child_node), std::default_random_engine(seed));
        }
        for(auto child : root->child_node)
        {
            assert(child !=NULL);
            if(MAX_DEPTH && child->max_depth<max_k) 
                continue;
            kcliquePlusPlus(n, rho, child, depth+1, max_k);
        }
        update_length(root);
        return; 
    }

    // at leave node. we can get the k-cliques: 
    if(GRAPH_REDUCE){
        // construct pset__
        pset__.clear();
        for(auto x:pset) if(Ck[x] >= ceil(sub_den)) pset__.push_back(x);

        // shuffle pset__
        if(RAMDOM){
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::shuffle(std::begin(pset__), std::end(pset__), std::default_random_engine(seed));
        }
        // updating this in every iteration incurs big overhead when kclique number is small. 

        // when the search scope is small enough, do not update this. 
        for(auto v : hset)   
            Ck_next[v] += nCr[pset__.size()][max_k - hset.size()];    
        for(auto v : pset__) 
            Ck_next[v] += nCr[pset__.size()-1][max_k - hset.size()-1]; 

    }else{
        // make a copy of pset. 
        // pset__ = pset; 
        // if(RAMDOM){
        //     unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //     std::shuffle(std::begin(pset__), std::end(pset__), std::default_random_engine(seed));
        // }
    }

    if(SAMPLE){
        num_to_enumerate = (int) nCr[pset__.size()][max_k - hset.size()] * sampling_ratio; 
    }

    if( hset.size() >= max_k){
        if(hset.size() == max_k)
        { 
            int arg_min = -1; 
            double min_rho = INT_MAX*1.0; 
            for(auto ver:hset){
                if(rho[ver]<min_rho){
                    arg_min = ver; 
                    min_rho = rho[ver]; 
                }                
            }
            rho[arg_min]++; 
            total_num_updates++;
        }
        update_length(root);
        return ; 
    } 

    // whether or not we apply graph reduction, we always operate on pset__, while pset is kept intact. 
    // this creates big overhead. 

    // at this point hset.size() < max_k.
    if(GRAPH_REDUCE){
        // operate on pset__ 
        if(hset.size()+pset__.size() >= max_k)
        {
            if(hset.size()+pset__.size() == max_k)
            {
                int arg_min = -1; 
                double min_rho = INT_MAX*1.0;
                for(auto ver: hset){
                    if(rho[ver] <min_rho){
                        arg_min = ver;  min_rho = rho[ver]; 
                    }            
                }
                for(auto ver: pset__){
                    if(rho[ver] <min_rho){
                        arg_min = ver; min_rho = rho[ver]; 
                    }            
                }
                rho[arg_min]++;
                total_num_updates++;
                update_length(root);
                return ; 
            }
        }
    }else{
        // operate on pset
        if(hset.size()+pset.size() >= max_k)
        {
            if(hset.size()+pset.size() == max_k)
            {
                int arg_min = -1; 
                double min_rho = INT_MAX*1.0;
                for(auto ver: hset){
                    if(rho[ver] <min_rho){
                        arg_min = ver;  
                        min_rho = rho[ver]; 
                    }            
                }
                for(auto ver: pset){
                    if(rho[ver] <min_rho){
                        arg_min = ver; 
                        min_rho = rho[ver]; 
                    }            
                }
                rho[arg_min]++;
                total_num_updates++;
                update_length(root);
                return ; 
            }
        }
    }

    printf("at leaf level \n");
    if(!BATCH_UPDATE){
        updateRho( num_to_enumerate, rho, max_k - hset.size());
    }else{
        // currently, we always turn on graph reduction and batch processing together 
        if(!GRAPH_REDUCE) {
            pset__=pset;
            if(RAMDOM){
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::shuffle(std::begin(pset__), std::end(pset__), std::default_random_engine(seed));
            }
        }
        batch_rho_update(0,num_to_enumerate, nCr, rho, max_k);
    }
    update_length(root);
    return;     
}

void show_arrays(int tmp, long long unsigned int *rho, int *hold_set, int *hold_len, int *pivot_set, int *pivot_len)
{
    // printf("p=%d, h=%d \n", *pivot_len, *hold_len);
    printf("hold: "); 
    for (int j=0; j<*hold_len; j++)  printf("%d ", hold_set[j] );
    printf("; "); 
    printf("pivot: "); 
    for (int j=0; j<*pivot_len; j++) printf("%d ", pivot_set[j] );
    printf("\n");     
}

void even_distribute_weight(double updates, int count1, long long unsigned int *rho, int *argmin_h_idxs)
{
    if(count1==1){
        int id_ = hset[argmin_h_idxs[0]];
        rho[id_]+=updates;
        return;
    }
    int mod_remainder = (int) updates % count1; 
    int increment = (int) updates / count1;
    int itr__ =0;
    for(int k=0;k<count1;k++){
        int id_ = hset[argmin_h_idxs[k]];
        if(itr__<mod_remainder){
            rho[id_]+=increment+1;
            itr__++;
        }else{
            rho[id_]+=increment;
        }      
    }
}
// total_num_updates+=count1

void pt(int xx){
    for(int yy = 0 ; yy<xx; yy++) printf("\t");
}

// always operate on pset__
void batch_rho_update(int depth, int remainder, double nCr[1001][401], long long unsigned int *rho, int max_k)
{
    if(depth > RECURSION_LIMIT){
        // setting a max depth for recursion
        // this needs to be re-examined!!! 
        updateRho( -1, rho, max_k - hset.size());
        return;
    }

    double gap ;
    double first = INT_MAX*1.0;     // the minimum
    double second = INT_MAX*1.0;    // the second smallest 
    int argmin=-1; 
    int min_label =-1; 
    // 
    double maxh = -INT_MAX*1.0, maxp = -INT_MAX*1.0; 
    //
    double minh=INT_MAX*1.0, minp=INT_MAX*1.0, secondh=INT_MAX*1.0, secondp=INT_MAX*1.0; 
    int argmin_h=-1, argmin_h_idx, argmin_p=-1, argmin_p_idx; 
    int argmin_h_idxs[100];
    // search the hold vertices.
    int count1 = 1, count2 = 1; 

    // these linear scans can be parralleized 
    pt(depth);printf("BatchUpdate\n");

    pt(depth); printf("hold: \n");
    for (int j=0; j<hset.size(); j++)
    {
        pt(depth); printf("v: %d  weight: %llu \n", hset[j], rho[ hset[j]]);
        maxh = maxh > rho[ hset[j]] ? maxh : rho[ hset[j]] ; 
        if(rho[ hset[j]] < minh){    // found a new minimum.
            secondh = minh ; 
            minh = rho[ hset[j]]; 
            count1 = 1;              // reset the number of minimum to 0.
            argmin_h_idxs[0] = j;
            argmin_h = hset[j];
            argmin_h_idx = j; 
        }else if(rho[hset[j]] < secondh && rho[hset[j]]!=minh  ){
            secondh = rho[hset[j]]; 
        }else if(rho[hset[j]]==minh){
            argmin_h_idxs[count1]=j;
            count1++; // used to check if the minimum is unique.
        }
    }
    pt(depth); printf("pivot: \n");
    for (int j=0; j<pset__.size(); j++){
        pt(depth); printf("v: %d  weight: %llu \n", pset__[j], rho[ pset__[j]]);
        maxp = maxp > rho[ pset__[j]] ? maxp : rho[ pset__[j]] ; 
        if(rho[ pset__[j]] < minp){
            secondp = minp ; 
            minp = rho[ pset__[j]]; 
            count2 = 1;
            argmin_p = pset__[j];
            argmin_p_idx = j; 
        }else if(rho[pset__[j]] < secondp && rho[pset__[j]]!=minp  ){
            secondp = rho[pset__[j]]; 
        }else if(rho[pset__[j]]==minp){
            count2++;
        }
    }  
    // printf("minh = %lf, minp = %lf \n", minh, minp);
    // count1  = |argminh|, count2 = |argminp|
    if(minh < minp){
        first = minh;
        second =  (secondh==INT_MAX) ? minp : (secondh < minp ? secondh : minp);
        argmin = argmin_h; 
        min_label = 0;
        gap = count1 * (second - minh);
    }
    if(minh >= minp){
        first = minp;
        second =  (secondp==INT_MAX) ? minh : (secondp < minh ? secondp : minh);
        argmin = argmin_p; 
        min_label = 1;
        if(count2>1) second = first ;
        gap = second - first +1;
    }
    // specialy handlings
    if(hset.size()+pset__.size() == max_k)
    {
        if(minh <  minp) rho[argmin_h]++;
        if(minh >= minp) rho[argmin_p]++;
        total_num_updates++;
        return;
    }
    if(min_label==1 && pset__.size() ==1)
    {
        rho[argmin] ++; 
        total_num_updates++;    
        return;   
    }
    if(remainder==-1){     
        if(min_label==0){
            pt(depth);printf("Case 1 \n");
            double updates = nCr[pset__.size()][max_k - hset.size()] ; 
            if(updates<=gap){
                // printf("// case 1.1 \n");
                even_distribute_weight(updates, count1, rho, argmin_h_idxs); 
                total_num_updates+=count1;
            }else{
                // printf("// case 1.2 \n"); 
                even_distribute_weight(gap, count1, rho, argmin_h_idxs); 
                total_num_updates+=count1;
                batch_rho_update(depth+1, updates - gap, nCr, rho, max_k);
            }
            return ;
        }
        if(min_label==1 && pset__.size()>1){
            pt(depth);printf("Case 2 \n");
            double partial_updates = nCr[pset__.size()-1][max_k - hset.size()-1] ; 
            if(partial_updates<=gap){
                // printf("// case 2.1 \n");
                rho[argmin] += partial_updates; 
                total_num_updates++; 
                
                if(hset.size()+pset__.size()-1>=max_k)
                {
                    // remove argmin_p from pivots 
                    bool switched = false ; 
                    if(argmin_p_idx!=pset__.size()-1){
                        int swapped = pset__[argmin_p_idx] ; 
                        pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                        pset__[pset__.size()-1] = swapped ; 
                        switched = true;
                    }
                    
                    // remove the last element 
                    int to_be_added_back = pset__[pset__.size()-1];
                    pset__.pop_back();
                    
                    batch_rho_update(depth+1, -1, nCr, rho, max_k);
                    pset__.push_back(to_be_added_back);
                    if(switched){
                        int swapped = pset__[argmin_p_idx] ; 
                        pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                        pset__[pset__.size()-1] = swapped ; 
                    }
                }
            }else{
                // printf("// case 2.2\n");
                rho[argmin] += gap;
                total_num_updates++;

                // insert argmin_p to hold: 
                hset.push_back(pset__[argmin_p_idx]); 

                // remove argmin_p from pivots 
                bool switched = false ; 
                if(argmin_p_idx!=pset__.size()-1){
                    int swapped = pset__[argmin_p_idx] ; 
                    pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                    pset__[pset__.size()-1] = swapped ; 
                    switched = true;
                }
                int to_be_added_back = pset__[pset__.size()-1];
                pset__.pop_back();

                // recursive call 1.
                batch_rho_update(depth+1, partial_updates - gap, nCr, rho, max_k);

                hset.pop_back();// then remove argmin at the end.

                // recursive call 2.
                if(hset.size()+pset__.size()>=max_k) batch_rho_update(depth+1, -1, nCr, rho, max_k);
                
                // add back the argmin_p:
                pset__.push_back(to_be_added_back); 

                if(switched){
                    // if switched, switch back
                    int swapped = pset__[argmin_p_idx] ; 
                    pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                    pset__[pset__.size()-1] = swapped ; 
                }
            }    
        }
    }
    if(remainder>0){
        // process remainder many cliques.
        if(min_label==0){
            pt(depth);printf("Case 3 \n");
            double updates = nCr[pset__.size()][max_k - hset.size()] ; 
            if(updates <= gap && updates <= remainder){
                // printf("// case 3.1 \n");
                even_distribute_weight(updates, count1, rho, argmin_h_idxs); 
                total_num_updates+=count1;
                batch_rho_update(depth+1, remainder-updates, nCr, rho, max_k);
            }
            if(updates <= gap && updates >  remainder){
                // printf("// case 3.2 \n");
                even_distribute_weight(remainder, count1, rho, argmin_h_idxs); 
                total_num_updates+=count1;
            }
            if(gap < updates && updates <= remainder ){
                // printf("// case 3.3 \n");
                even_distribute_weight(gap, count1, rho, argmin_h_idxs); 
                total_num_updates+=count1;
                batch_rho_update(depth+1, remainder-gap, nCr, rho, max_k);
            }
            if(gap <= remainder && remainder < updates ){
                // printf("// case 3.4 \n");
                even_distribute_weight(gap, count1, rho, argmin_h_idxs); 
                total_num_updates+=count1;
                batch_rho_update(depth+1, remainder-gap, nCr, rho, max_k);
            }
            if(gap < updates && remainder < gap ){
                // printf("// case 3.5 \n");
                even_distribute_weight(remainder, count1, rho, argmin_h_idxs); 
                total_num_updates+=count1;
            }
        }
        if(min_label==1 && pset__.size()>1)
        {
            pt(depth);printf("Case 4 \n");
            double partial = nCr[pset__.size()-1][max_k - hset.size()-1] ; 
            if(partial <= gap && partial <= remainder)
            {
                // printf("// case 3.6 \n");
                rho[argmin] += partial; 
                total_num_updates++;
                
                if(pset__.size()+hset.size()-1>=max_k)
                {
                    // remove argmin from pivots
                    bool switched = false ; 
                    if(argmin_p_idx!=pset__.size()-1){
                        int swapped = pset__[argmin_p_idx] ; 
                        pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                        pset__[pset__.size()-1] = swapped ; 
                        switched = true;
                    }
                    int to_be_added_back = pset__[pset__.size()-1];
                    pset__.pop_back();

                    batch_rho_update(depth+1, remainder-partial, nCr, rho, max_k);
                    // add back the argmin_p:
                    pset__.push_back(to_be_added_back); 
                    
                    if(switched){
                        int swapped = pset__[argmin_p_idx] ; 
                        pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                        pset__[pset__.size()-1] = swapped ; 
                    }
                }
            }
            if(partial <= gap && partial > remainder){
                // printf("// case 3.7 \n");
                rho[argmin]+= remainder; 
                total_num_updates++;
            }
            if(gap < partial && partial <= remainder )
            {
                // printf("// case 3.8 \n");
                rho[argmin] += gap; 
                total_num_updates++;
                // Now process partial-gap cliques containig argmin:

                // move argmin from pivots to holds;
                hset.push_back(pset__[argmin_p_idx]); 

                // remove argmin_p from pivots 
                bool switched = false ; 
                if(argmin_p_idx!=pset__.size()-1){
                    int swapped = pset__[argmin_p_idx] ; 
                    pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                    pset__[pset__.size()-1] = swapped ; 
                    switched = true;
                }
                int to_be_added_back = pset__[pset__.size()-1];
                pset__.pop_back();
                
                batch_rho_update(depth+1, partial-gap, nCr, rho, max_k);

                hset.pop_back();// remove argmin from holds;

                if(pset__.size()+hset.size()-1>=max_k){
                    batch_rho_update(depth+1, remainder-partial, nCr, rho, max_k);
                }
                pset__.push_back(to_be_added_back);

                if(switched){
                    int swapped = pset__[argmin_p_idx] ; 
                    pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                    pset__[pset__.size()-1] = swapped ; 
                }
            }
            if( gap <= remainder && remainder < partial)
            {
                // printf("// case 3.9 \n");
                rho[argmin] += gap; 
                total_num_updates++;
                // process remainder - gap cliques containing argmin.

                // move argmin from pivots to holds;
                hset.push_back(pset__[argmin_p_idx]);

                // remove argmin_p from pivots 
                bool switched = false ; 
                if(argmin_p_idx!=pset__.size()-1){
                    int swapped = pset__[argmin_p_idx] ; 
                    pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                    pset__[pset__.size()-1] = swapped ; 
                    switched = true;
                }
                int to_be_added_back = pset__[pset__.size()-1];
                pset__.pop_back();

                batch_rho_update(depth+1, remainder-gap, nCr, rho, max_k);
                hset.pop_back();
                
                // add back 
                pset__.push_back(to_be_added_back); 
                if(switched){
                    int swapped = pset__[argmin_p_idx] ; 
                    pset__[argmin_p_idx] = pset__[pset__.size()-1]; 
                    pset__[pset__.size()-1] = swapped ; 
                }
            }
            if(gap < partial && remainder<gap)
            {
                // printf("// case 3.10 \n");
                assert(remainder<gap);// need to process remainder many cliques containing partial.
                rho[argmin]+= remainder;
                total_num_updates++;
            } 
        }
    }
}

int get_sampled_cliques(int remainder, int r)
{
    int data[r];
    int *visited = (int *) Calloc(1, sizeof(int));
    *visited = 0;
    int result = get_sampled_cliques_helper(remainder, visited, data, 0, pset__.size()-1, 0, r);
    Free(visited);
    return result;
}

int get_sampled_cliques_helper(int remainder, int *visited, int *data, int start, int end, int index, int r)
{
    if(remainder==*visited) 
        return 0; 

    if(sampled_kcliques.size() >= 1.1*sample_size) 
        return 0; 
    
    if (index == r)
    {
        // sampled_kcliques
        vector<int> curr_cliq = hset; 
        for(int j=0;j<r;j++) curr_cliq.push_back(data[j]);
        sampled_kcliques.push_back(curr_cliq);

        // while visiting this k-clique. use it to update the rho array. 
        
        // assert(cliq_idx = sampled_kcliques.size());
        (*visited)++; 
        return 1;
    }
    int tmp=0;
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = pset__[i];
        tmp+=get_sampled_cliques_helper(remainder, visited, data, i+1, end, index+1, r);
    }
    return tmp;
}


int printCombination_exact(int r)
{
    int data[r];
    return combinationUtil_exact(data, 0, pset.size()-1, 0, r);
}

int combinationUtil_exact(int *data, int start, int end, int index, int r)
{
    // Current combination is ready to be printed, print it
    if (index == r)
    {
        // now operating on the cliq_idx clique.
        int max_k = r + hset.size(); 
        // for each entry need to record the vertex. 
        int node_idx = cliq_idx * max_k;

        for(auto ver:hset){
            ck_exact[node_idx] = ver ; 
            node_idx++;
        }
        for(int j=0;j<r;j++){
            ck_exact[node_idx] = data[j] ; 
            node_idx++;
        }
        cliq_idx++; 
        return 1;
    }
    int tmp=0;
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = pset[i];
        tmp = tmp+combinationUtil_exact(data, i+1, end, index+1, r);
    }
    return tmp;
}

// combination , this function can be used to compute Ck(MCS). 
int printCombination(int r)
{
    int data[r];
    return combinationUtil(data, 0, pset.size()-1, 0, r);
}
 
int combinationUtil(int *data, int start, int end, int index, int r)
{
    // Current combination is ready to be printed, print it
    if (index == r)
    {
        // printf("Clique: ");
        // for(auto ver:hset) printf("%d ", ver);
        // printf(" ; ");
        // for (int j=0; j<r; j++) printf("%d ", data[j]);
        // printf("\n");
        vector<int> tmp;
        
        int num_dominated = 0;
        
        for(auto ver:hset){
            if(MC_lookup[ver]){
                num_dominated++;
            }
            tmp.push_back(ver);
        }
        for(int j=0;j<r;j++){
            if(MC_lookup[data[j]]) {
                num_dominated++;
            }
            tmp.push_back(data[j]);
        }
        if(num_dominated==hset.size()+r){
            total_num_updates++;
        }

        // printf("num dominated = %d \n", num_dominated);
        // assert(num_dominated==hset.size()+r || num_dominated==0);
        
        all_mc.push_back(tmp);

        return 1;
    }

    // otherwise make a recursive call.
    int tmp=0;
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = pset[i];
        tmp = tmp+combinationUtil(data, i+1, end, index+1, r);
    }
    return tmp;
}

// compute yi for the cliques in MC
void listYi_MC(int r, int *level, int *y)
{
    int data[r];
    combineYi_MC( data, 0, curr_mc.size()-1, 0, r, level, y);
}
void combineYi_MC(int *data, int start, int end, int index, int r, int *level, int *y)
{
    // Current combination is ready to be printed, print it
    if (index == r)
    {
        // go through these vertices, find the first one with the max level
        printf("clique: "); 
        for (int j=0; j<r; j++){
            printf(" %d ", data[j]);
        }
        printf("\n"); 

        int arg_max = -1; 
        int max_level = -1;;
        // just search the pivot vertices        
        for (int j=0; j<r; j++){
            if(level[data[j]]>max_level){
                arg_max = data[j]; 
                max_level = level[data[j]]; 
            }
        }
        // at this point 
        y[arg_max]++;
        return;
    }
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = curr_mc[i];
        combineYi_MC(data, i+1, end, index+1, r,level, y);
    }
}

// compute yi:
void listForYi(int r, int *level, unsigned long long int *y)
{
    int data[r];
    int argmaxhold = hset[0];
    int maxhold = level[argmaxhold];
    for (int j=1; j<hset.size(); j++){
        if(level[hset[j]] > maxhold ){
            argmaxhold = hset[j]; 
            maxhold = level[hset[j]]; 
        }
    } 
    if(GRAPH_REDUCE){
        combineForYi( data, 0, pset__.size()-1, 0, r, level, y, &argmaxhold, &maxhold);
    }else{
        combineForYi( data, 0, pset.size()-1, 0, r, level, y, &argmaxhold, &maxhold);
    }
}
void combineForYi(int *data, int start, int end, int index, int r, int *level, unsigned long long int *y, int *argmaxhold, int *maxhold)
{
    // Current combination is ready to be printed, print it
    if (index == r)
    {
        // go through these vertices, find the first one with the max level
        int arg_max = *argmaxhold ; 
        int max_level = *maxhold;
        // search the pivot vertices        
        for (int j=0; j<r; j++){
            if(level[data[j]]>max_level){
                arg_max = data[j]; 
                max_level = level[data[j]]; 
            }
        }
        y[arg_max]++;
        yyy++;
        return;
    }
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        if(GRAPH_REDUCE){
            data[index] = pset__[i];
        }else{
            data[index] = pset[i];
        }
        combineForYi(data, i+1, end, index+1, r,level, y, argmaxhold, maxhold);
    }
}
// (1) when graph reduction is applied only visit the kcliques in G'
// this function should have a sampling switch. 
// when sampling, graph reduce is set to false.
void compute_yi(treeNode *root, int depth, int max_k, int *level, unsigned long long int *y)
{
    if(root==NULL) 
        return ;

    if(root->max_depth<max_k) 
        return ; 

    // add the current node (and its label) to the container. 
    if(root->id>=0){
        if(root->label==0){

            if(GRAPH_REDUCE && Ck[root->id]< ceil(sub_den))
                return ;

            if(SAMPLE && (!ver_exists_in_approx[root->id]))
                return ;

            hset.push_back(root->id); 
        }else{
            pset.push_back(root->id); 
        }
    }
    // recurse on children if any 
    if(root->child_node.size()>0){
        for(auto child : root->child_node) 
            compute_yi(child, depth+1, max_k,level,y);
    }
    else{ 

        if(GRAPH_REDUCE){
            pset__.clear();
            for(auto x:pset) if(Ck[x] >= ceil(sub_den)) pset__.push_back(x);
        }

        if(SAMPLE){
            // pset__.clear();
            int pset_size = 0;
            for(auto x:pset) if(ver_exists_in_approx[x]) pset_size++;
            yyy+= nCr[pset_size][max_k - hset.size()]; 
        }

        // at leave node
        if( !SAMPLE && hset.size() <= max_k) 
            listForYi(max_k - hset.size(), level, y);

    }
    update_length(root);
}

void list_enumeration_cliques_in_approximate_solution(int r)
{

    if(hset.size()+pset__.size() < r+hset.size()) return ;

    // printf("h:");
    // for(auto v:hset) printf("%d ", v);
    // printf("\n");
    // printf("p:");
    // for(auto v:pset__) printf("%d ", v);
    // printf("\n\n");

    assert(hset.size()+pset__.size() >= r+hset.size()); 

    int data[r];
    combine_enumeration_cliqs_N_approx( data, 0, pset__.size()-1, 0, r);
}
void combine_enumeration_cliqs_N_approx(int *data, int start, int end, int index, int r)
{
    // Current combination is ready to be printed, print it
    if (index == r)
    {
        // go through these vertices, find the first one with the max level
        m__++;
        // printf("%ld clique visited \n", m__ );
        
        Network::Vertex v = curr_network->AddVertex();
        for(auto ver:hset){
            curr_network->AddEdge(v, R[level[ver]], num_nodes);
        }
        // for(auto ver:pset__){
        for(int j=0;j<r;j++){
            int ver = data[j];
           curr_network->AddEdge(v, R[level[ver]], num_nodes);
        }
        curr_network->AddEdge(*source_ptr, v, num_nodes);
        
        return;
    }
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = pset__[i];
        combine_enumeration_cliqs_N_approx(data, i+1, end, index+1, r);
    }
}

void enumerate_cliques_in_approximate_solution(treeNode *root, int num_nodes, int max_k, int *level)
{
    if(root==NULL) 
        return ;

    if(root->max_depth<max_k) 
        return ; 

    if(root->id>=0){
        if(root->label==0){
            if( level[root->id] >= num_nodes )
                return ;
            hset.push_back(root->id); 
        }else{
            pset.push_back(root->id); 
        }
    }
    // recurse on children if any 
    if(root->child_node.size()>0){
        for(auto child : root->child_node) enumerate_cliques_in_approximate_solution(child, num_nodes, max_k,level);
    }
    else{ 
        // at leave node
        pset__.clear();
        for(auto x:pset) if(level[x] < num_nodes) pset__.push_back(x);

        if( hset.size() <= max_k)
            list_enumeration_cliques_in_approximate_solution(max_k - hset.size());
    }
    update_length(root);
}

void destroyTreeFromRoot(treeNode *root)
{
    if(root==NULL) return;

    // printf("destorying %d\n", root->id);
    if(root->child_node.size()>0){
        for(auto child : root->child_node){
            destroyTreeFromRoot(child); 
        }
        root->child_node.clear();
    }
    Free(root);
    return;
}
// cut_point by default is root. 
// status is default true. one a pushed vertex is not in the MC, then status becomes false and returns.

int find_func(int i)
{
    // naive implementation
    // if(group[i] == i) return i;
    // return find_func(group[i]);
    // path compression
    if(group[i]!=i){
        group[i] = find_func(group[i]); 
    }
    return group[i];
}

void union_func(){
    // union the roots pointed by hset and pset.
    // go through all vertices in hset and pset and find the root of each subtree.
    vector<int> roots;

    // go through all vertices in hset and pset. find the tree with maximum rank. 
    int argmax = -1; 
    int maxrank = 0;

    for(auto v:hset){
        int root_v = find_func(v);
        if(uf_rank[root_v] > maxrank){
            argmax = root_v;
            maxrank = uf_rank[root_v];
        }
        roots.push_back(root_v); 
    }

    for(auto v:pset){
        int root_v = find_func(v);
        if(uf_rank[root_v] > maxrank){
            argmax = root_v;
            maxrank = uf_rank[root_v];
        }
        roots.push_back(root_v); 
    }
    // attach all other trees to the tree rooted at argmax.
    for(auto root : roots){
        if(root!=argmax){
            group[root] = argmax; 
            if(uf_rank[root]==uf_rank[argmax]){
                uf_rank[argmax]++;
            }
        }
    }
}


void UFKclique(treeNode *root, int depth, int max_k)
{

    if(root==NULL) return ;

    if(root->max_depth<max_k ) return ; 

    // add the current node (and its label) to the container. 
    if(root->id>=0){
        if(root->label==0){
            hset.push_back(root->id); 
        }else{
            pset.push_back(root->id); 
        }
    }

    if(root->child_node.size()>0){
        // recurse on its children if any
        for(auto child : root->child_node){
            UFKclique(child, depth+1, max_k);
        }

    }else{ 
        // at leave node. we can get the k-cliques: 
        if( hset.size() <= max_k) {

            /*
            // merge the vertices in this branch.
            int groupid = group[hset[0]]; 
            // reset lookup table
            fill(MC_lookup.begin(),MC_lookup.end(), false); 
            for(auto xx : hset) MC_lookup[group[xx]] = true;
            for(auto yy : pset) MC_lookup[group[yy]] = true; 
            for(int i=0;i<MC_lookup.size();i++){
                if( MC_lookup[ group[i]  ] ){
                    group[i] = groupid ; 
                }
            }
            */


            // go though hset and pset, find their groups
            // for all vertices within these groups. give them the same group id. 
            union_func();
        }
    }
    if(root->id>=0){
        if(root->label==0){
            hset.resize(hset.size()-1); 
        }else{
            pset.resize(pset.size()-1); 
        }
    }
}

// when sampling is on. we use this to collect kcliques. 
// max_k is the queried k 
// base_k is the smallest clique we find kclique in. 
// (1) number of cliques sampled capped at 1.1 * 10^7. 
// (2) added rounding
long int ListRootToLeavePath(treeNode *root, int depth, int max_k, int base_k)
{

    if(root==NULL) 
        return 0;
    
    if(sampled_kcliques.size() >= 1.1*sample_size)
        return 0;

    // when only visiting depth >= SCT_K. 
    // since we unavoidally build some trees with depth < SCT_K. 
    // we can go through these shallow branches when listing
    // kq-cliques. 
    if(root->max_depth<max_k ) 
        return 0; 

    // add the current node (and its label) to the container. 
    if(root->id>=0){
        if(root->label==0){

            if(GRAPH_REDUCE && Ck[root->id]< ceil(sub_den)) {
                // exit(0);
                return 0;
            }
            hset.push_back(root->id); 
        }else{
            pset.push_back(root->id); 
        }
    }
    long int tmp = 0 ; 
    if(root->child_node.size()>0){
        // recurse on its children if any
        // when sampling. we need not to shuffle the branches because we need to visit all branches anyways 
        for(auto child : root->child_node){
            tmp += ListRootToLeavePath(child, depth+1, max_k, base_k);
        }
    }else{ 
        // at leave node. we can get the k-cliques: 
        if( hset.size() <= max_k) {
            long int curr=0 ; 
            if(GRAPH_REDUCE){
                // make pset_ array 
                pset__.clear();
                for(auto x:pset) if(Ck[x] >= ceil(sub_den)) pset__.push_back(x);
                // update Ck in hset and pset_.
                for(auto v : hset)   Ck_next[v] += nCr[pset__.size()][max_k - hset.size()];    
                for(auto v : pset__) Ck_next[v] += nCr[pset__.size()-1][max_k - hset.size()-1]; 
            }
            if(SAMPLE){
                pset__ = pset; 
                double total_work = nCr[pset__.size()][max_k - hset.size()]; 
                
                // int num_to_enumerate=   (int) (total_work* sampling_ratio); 
                int num_to_enumerate=   (int) (total_work* sampling_ratio+0.5); // rounding

                // solution: 
                // (1) if num_to_enumerate == 0, round it up to 1. 
                // (2) visit the deeper trees first. 

                if(num_to_enumerate>0){
                    // if the number we sample is positive, then we shuffle the pset_ array
                    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                    std::shuffle(std::begin(pset__), std::end(pset__), std::default_random_engine(seed));
                    curr = get_sampled_cliques(num_to_enumerate, max_k - hset.size()); 
                }else{
                    curr = 0; 
                }

                // printf("total = %lf \t to-be-sampled = %d  \t actually-sampled = %d \n", total_work, num_to_enumerate, curr); 

                // if(curr==0){
                //     printf("sampled cliques = %d \n", sampled_kcliques.size());
                // }
                assert(curr <= nCr[pset__.size()][max_k - hset.size()]);
            }
            tmp += curr ; 
        }
    }
    // root->num_cliques = tmp; // record the number of cliques into the tree.
    update_length(root);
    return tmp;
}

void MC_summary(){

    for(auto cliq : all_mc){
        
        printf("cliq: ");
        
        sort(cliq.begin(), cliq.end());

        for(auto v : cliq) printf("%d ", v);
        
        printf("\n\n");
    }

}


void showTree(treeNode *root, int depth, int max_k, long long unsigned int *rho, bool *lookup)
{
    if(root==NULL) return ; 

    if(root->max_depth<max_k) return; 


    for(int i=0;i<depth ; i++) printf("\t"); 

    if(root->label==1){
        printf("node %d: pivot\n", root->id); 
    }else{
        printf("node %d: hold \n",  root->id); 
    }

    // if(root->id>=0) {
    //     // printf("vertex %d\n", root->id);     
    //     if(rho)    rho[root->id] = 1;
    //     if(lookup && root->max_depth == 102 ) lookup[root->id]=true;
    // }

    if(root->child_node.size()==0){
        // at leave node.
        // if(root->label==1){

            treeNode *tmpParent = root->parent ; 

            // printf("vertex %d at leaf, parent = %d, position = %d \n", root->id, root->parent->id, pos);

            // this is creating problems
            // destroyTreeFromRoot(root);
            // to properly delete this subtree. 
            // root=NULL;

            // printf("after deletion, check the child list of parent\n");
            // for(auto x : tmpParent->child_node){
            //     if(x!=NULL){
            //         printf("sibling = %d \n", x->id);
            //     }else{
            //         printf("empty branch \n");
            //     }
            // }
            // tmpParent->child_node.clear();

        // }

        return ;
    }else{
        for(auto child : root->child_node){
        // for(int k=0; k < root->child_node.size(); k++){
            // treeNode* child = root->child_node[k];
            assert(child!=NULL);
            showTree(child,depth+1, max_k, rho,lookup); 
        }
        // if ons of the subtrees rooted at its children is deleted, resize the child_node array here.
    }

}

void showemptybranches(){
    printf("num of empty branches = %d \n", num_of_empty_branches);
    printf("size of pruned tree = %d \n", size_of_pruned_tree); 
}

// Maxflow related functions
Network::Network(): g(), rev(get(edge_reverse, g)) {

}

Network::Vertex Network::AddVertex() {
  return add_vertex(g);
}

void Network::AddEdge(Vertex &v1, Vertex &v2, const long capacity) {
  Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
  Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
  put(edge_capacity, g, e1, capacity);
  rev[e1] = e2;
  rev[e2] = e1;
}

unsigned long long Network::MaxFlow(Vertex &s, Vertex &t) {
  return push_relabel_max_flow(g, s, t);
}

/*
bool stable_verify(int n, double * rho_init){

    double min__ = INT_MAX*1.0; 
    double max__ = -1.0*INT_MAX ; 

    double sum___ = 0; 
    double sum_all = 0; 
    for(int i=0;i<n;i++){
        if(MC_lookup[i]){ // covered 
            min__ = min__ < rho_init[i] ? min__ : rho_init[i]; 
            // printf("covered rho : % lf \n", rho_init[i]);
            sum___+=  rho_init[i] ; 
        }else{ // uncovered 
            max__ = max__ > rho_init[i] ? max__ : rho_init[i]; 
        }
        if(rho_init[i]>=0) sum_all +=  rho_init[i] ; 
        // printf(" rho : % lf \n", rho_init[i]);
    }

    printf("min rho = %lf (density = %lf ),  max rho = %lf  , sum = % lf \n",  min__, sum___/stable_set.size(),  max__, sum___) ; 

    printf("sum all rho = %lf \n", sum_all);

    if(min__ > max__){
        return true;
    }else{
        return false;
    }
}


void stable_combine(int r, double * rho_init)
{
    int data[r];
    stable_combine_until(data, 0, pset.size()-1, 0, r, rho_init);
}
 
void stable_combine_until(int *data, int start, int end, int index, int r, double * rho_init)
{

    if (index == r)
    {
        // printf("Clique: ");
        // for(auto ver:hset) printf("%lf ", rho_init[ver]);
        // printf(" ; ");
        // for (int j=0; j<r; j++) printf("%lf ", rho_init[ data[j]]);
        // printf("\n");

        double min__ = INT_MAX*1.0; 
        int argmin = -1; 

        // double min_covered = INT_MAX*1.0; 
        // int argmin_covered = -1; 

        int num_covered = 0; 

        vector<int> not_covered_ver ; 

        for(auto ver:hset){
            if(!MC_lookup[ver]){
                if(rho_init[ver] < min__){
                    min__ = rho_init[ver];
                    argmin = ver; 
                }
                not_covered_ver.push_back(ver);
            }else{
                // if(rho_init[ver] < min__){
                //     min_covered = rho_init[ver];
                //     argmin_covered = ver; 
                // }
                num_covered++;
            }
        }
        for(int j=0;j<r;j++){
            int ver = data[j];
            if(!MC_lookup[ver]){
                if(rho_init[ver] < min__){
                    min__ = rho_init[ver];
                    argmin = ver; 
                }
                not_covered_ver.push_back(ver);
            }else{
                // if(rho_init[ver] < min__){
                //     min_covered = rho_init[ver];
                //     argmin_covered = ver; 
                // }
                num_covered++;
            }
        }
        // printf("\n");
        double cliquesize =hset.size()+r ;  
        if(num_covered==cliquesize) {
            // need to distribute according to 
            for(auto ver:hset)      rho_init[ver] += 1.0/cliquesize; 
            for(int j=0;j<r;j++)    rho_init[data[j]] += 1.0/cliquesize; 

            // rho_init[argmin_covered] += 1.0;

        }else{
            // for(auto ver : not_covered_ver) rho_init[ver] += 1.0/not_covered_ver.size();  // this is the worst 

            double denominator = 0;
            // for(auto ver : not_covered_ver) denominator+= 1.0/( per_ver_cnt[ver][cliquesize]  ); 
            // for(auto ver : not_covered_ver) rho_init[ver] += 1.0/( per_ver_cnt[ver][cliquesize]*denominator );

            for(auto ver : not_covered_ver){
                double cnt_factor = per_ver_cnt[ver][cliquesize]+1; 
                if(ver == argmin){
                    denominator+= 2.0/( cnt_factor *(rho_init[ver]+1)  ); 
                }else{
                    denominator+= 1.0/( cnt_factor *(rho_init[ver]+1)  ); 
                }   
            }
            double incre_sum = 0 ;
;
            for(auto ver : not_covered_ver){
                double incre ;
                double cnt_factor = per_ver_cnt[ver][cliquesize]+1; 
                if(ver == argmin){
                    incre = 2.0/( cnt_factor*(rho_init[ver]+1)*denominator ); 
                }else{
                    incre = 1.0/( cnt_factor*(rho_init[ver]+1)*denominator ); 
                }   
                assert(incre >=0 );
                rho_init[ver] += incre;
                incre_sum+=incre;
            }

            // printf("incre sum = %lf \n", incre_sum);
            // (1.0/per_ver_cnt[ver][cliquesize])/denominator; 
            

            // make these as even as possible. 
            // rho_init[argmin] += 1.0;
        }
        // printf("num covered = %d \n", num_covered);

        return ;
    }

    // otherwise make a recursive call.
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = pset[i];
        stable_combine_until(data, i+1, end, index+1, r, rho_init);
    }
}

void stable_init(treeNode *root, int depth, int max_k, double * rho_init)
{

    if(root==NULL) return ;

    if(root->max_depth<max_k ) return ; 

    // add the current node (and its label) to the container. 
    if(root->id>=0){

        if(!KCC_lookup[find_func(root->id)]) return; 
        
        if(root->label==0){
            hset.push_back(root->id); 
        }else{
            pset.push_back(root->id); 
        }
    }

    if(root->child_node.size()>0){
        // recurse on its children if any
        for(auto child : root->child_node){
            stable_init(child, depth+1, max_k, rho_init);
        }
    }else{ 
        // at leave node. we can get the k-cliques: 
        // printf("at leave node\n");
        if( hset.size() <= max_k) {

            int num_covered = 0 ;
            for(auto xx : hset) if(MC_lookup[xx]) num_covered++;
            for(auto xx : pset) if(MC_lookup[xx]) num_covered++;

            // printf("If all or none, simply distribute weight.\n"); 
            if( num_covered== hset.size()+pset.size() || num_covered == 0 ){
                for(auto xx : hset ) rho_init[xx] += nCr[pset.size()][max_k - hset.size()]/max_k ; 
                for(auto xx : pset ) rho_init[xx] += nCr[pset.size()-1][max_k - hset.size()-1]/max_k ; 
            }else{
                stable_combine(max_k - hset.size(), rho_init);  
            }
        }
    }

    if(root->id>=0){
        if(root->label==0){
            hset.resize(hset.size()-1); 
        }else{
            pset.resize(pset.size()-1); 
        }
    }
}
*/
/*! \brief

    \param list an input graph, represented as an array of linked lists of integers

    \param size the number of vertices in the graph

    \return the degeneracy of the input graph.
*/

int computeDegeneracy(LinkedList** list, int size)
{
    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    LinkedList** verticesByDegree = (LinkedList**) Calloc(size, sizeof(LinkedList*));

    // array of lists of vertices, indexed by degree
    Link** vertexLocator = (Link**) Calloc(size, sizeof(Link*));

    int* degree = (int*) Calloc(size, sizeof(int));

    for(i = 0; i < size; i++)
    {
        verticesByDegree[i] = createLinkedList();
    }

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<size; i++)
    {
        degree[i] = length(list[i]);
        vertexLocator[i] = addFirst(verticesByDegree[degree[i]], ((int) i));
    }
    
    int currentDegree = 0;

    int numVerticesRemoved = 0;

    while(numVerticesRemoved < size)
    {
        if(!isEmpty(verticesByDegree[currentDegree]))
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int vertex = (int)(size_t) getFirst(verticesByDegree[currentDegree]);

            deleteLink(vertexLocator[vertex]);

            degree[vertex] = -1;

            LinkedList* neighborList = list[vertex];

            Link* neighborLink = neighborList->head->next;

            while(!isTail(neighborLink))
            {
                int neighbor = (int)(size_t) neighborLink->data;

                if(degree[neighbor]!=-1)
                {
                    deleteLink(vertexLocator[neighbor]);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        vertexLocator[neighbor] = 
                            addFirst(verticesByDegree[degree[neighbor]], 
                                     (int) neighbor);
                    }
                }

                neighborLink = neighborLink->next;
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }

    }

    for(i = 0; i<size;i++)
    {
        destroyLinkedList(verticesByDegree[i]);
    }

    Free(vertexLocator);
    Free(verticesByDegree);
    Free(degree);

    return degeneracy;
}

/*! \brief

    \param list an input graph, represented as an array of linked lists of integers

    \param size the number of vertices in the graph

    \return an array of NeighborLists representing a degeneracy ordering of the vertices.

    \see NeighborList
*/

NeighborList** computeDegeneracyOrderList(LinkedList** list, int size)
{

#ifdef DEBUG
    printf("degeneracy is %d\n", computeDegeneracy(list, size));
#endif

    NeighborList** ordering = (NeighborList**)Calloc(size, sizeof(NeighborList*));

    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    LinkedList** verticesByDegree = (LinkedList**) Calloc(size, sizeof(LinkedList*));

    // array of lists of vertices, indexed by degree
    Link** vertexLocator = (Link**) Calloc(size, sizeof(Link*));

    int* degree = (int*) Calloc(size, sizeof(int));

    for(i = 0; i < size; i++)
    {
        verticesByDegree[i] = createLinkedList();
        ordering[i] = (NeighborList*)Malloc(sizeof(NeighborList));
        ordering[i]->earlier = createLinkedList();
        ordering[i]->later = createLinkedList();
    }

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree

    for(i=0; i<size; i++)
    {
        degree[i] = length(list[i]);
        //printf("degree[%d] = %d\n", i, degree[i]);
        vertexLocator[i] = addFirst(verticesByDegree[degree[i]], ((int) i));
    }
    
    int currentDegree = 0;

    int numVerticesRemoved = 0;

    while(numVerticesRemoved < size)
    {
        if(!isEmpty(verticesByDegree[currentDegree]))
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int vertex = (int)(size_t) getFirst(verticesByDegree[currentDegree]);

            deleteLink(vertexLocator[vertex]);

            ordering[vertex]->vertex = vertex;
            ordering[vertex]->orderNumber = numVerticesRemoved;

            degree[vertex] = -1;

            LinkedList* neighborList = list[vertex];

            Link* neighborLink = neighborList->head->next;

            while(!isTail(neighborLink))
            {
                int neighbor = (int)(size_t) neighborLink->data;
                //printf("Neighbor: %d\n", neighbor);

                if(degree[neighbor]!=-1)
                {
                    deleteLink(vertexLocator[neighbor]);
                    addLast(ordering[vertex]->later, (int)neighbor);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        vertexLocator[neighbor] = 
                            addFirst(verticesByDegree[degree[neighbor]], 
                                     (int) neighbor);
                    }
                }
                else
                {
                    addLast(ordering[vertex]->earlier, (int) neighbor);
                }

                neighborLink = neighborLink->next;
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }

    }

    for(i = 0; i<size;i++)
    {
        destroyLinkedList(verticesByDegree[i]);
    }

    Free(vertexLocator);
    Free(verticesByDegree);
    Free(degree);

    return ordering;
}

/*! \brief

    \param list an input graph, represented as an array of linked lists of integers

    \param size the number of vertices in the graph

    \return an array of NeighborListArrays representing a degeneracy ordering of the vertices.

    \see NeighborListArray
*/

NeighborListArray** computeDegeneracyOrderArray(LinkedList** list, int size, int max_k)
{
    NeighborList** ordering = (NeighborList**)Calloc(size, sizeof(NeighborList*));

    int i = 0;

    int degeneracy = 0;
    
    // array of lists of vertices, indexed by degree
    LinkedList** verticesByDegree = (LinkedList**) Calloc(size, sizeof(LinkedList*));

    // array of lists of vertices, indexed by degree
    Link** vertexLocator = (Link**) Calloc(size, sizeof(Link*));

    int* degree = (int*) Calloc(size, sizeof(int));

    for(i = 0; i < size; i++)
    {
        verticesByDegree[i] = createLinkedList();
        ordering[i] = (NeighborList*)Malloc(sizeof(NeighborList));
        ordering[i]->earlier = createLinkedList();
        ordering[i]->later = createLinkedList();
    }

    // fill each cell of degree lookup table
    // then use that degree to populate the 
    // lists of vertices indexed by degree
    // fprintf(stderr, "Ordering is:\n" );
    for(i=0; i<size; i++)
    {
        degree[i] = length(list[i]);
        // fprintf(stderr, "i=%d, degree[i]=%d\n",i, degree[i] );
        vertexLocator[i] = addFirst(verticesByDegree[degree[i]], ((int) i));
    }
    // fprintf(stderr, "Ordering is:\n" );
    int currentDegree = 0;

    int numVerticesRemoved = 0;

    // track corenumbers
    while(numVerticesRemoved < size)
    {
        if(!isEmpty(verticesByDegree[currentDegree]))
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int vertex = (int)(size_t) getFirst(verticesByDegree[currentDegree]);

            corenumber[vertex] = degeneracy ; 
            // printf("vertex=%d, degeneracy=%d\n", vertex, degeneracy); 

            deleteLink(vertexLocator[vertex]);

            ordering[vertex]->vertex = vertex;
            ordering[vertex]->orderNumber = numVerticesRemoved;
            // printf("vertex = %d, degree[vertex]=%d\n", vertex, degree[vertex]);
            degree[vertex] = -1;

            LinkedList* neighborList = list[vertex];

            Link* neighborLink = neighborList->head->next;

            while(!isTail(neighborLink))
            {
                int neighbor = (int)(size_t) (neighborLink->data);
                // printf("neighbor = %d, degree[neighbor] = %d, ordering[neighbor]->orderNumber=%d\n",neighbor, degree[neighbor], ordering[neighbor]->orderNumber);
                if(degree[neighbor]!=-1)
                {
                    deleteLink(vertexLocator[neighbor]);
                    addLast(ordering[vertex]->later, (int)neighbor);

                    // printf("In. neighbor = %d, degree[neighbor] = %d, ordering[neighbor]->orderNumber=%d\n",neighbor, degree[neighbor], ordering[neighbor]->orderNumber);
                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        vertexLocator[neighbor] = 
                            addFirst(verticesByDegree[degree[neighbor]], 
                                     (int) neighbor);
                    }
                }
                else
                {
                    addLast(ordering[vertex]->earlier, (int) neighbor);
                }

                neighborLink = neighborLink->next;
            }

            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }

    }
    
    NeighborListArray** orderingArray = (NeighborListArray**)Calloc(size, sizeof(NeighborListArray*));

    for(i = 0; i<size;i++)
    {
        // only allocate for the degeneracy fullfilling vertices. 
        // if(CN_PRUNE && corenumber[ordering[i]->vertex]+1 < max_k){         
        //     continue;
        // }
        // printf("vertex = %d, cn = %d\n",ordering[i]->vertex,  corenumber[ordering[i]->vertex]);

        orderingArray[i] = (NeighborListArray*)Malloc(sizeof(NeighborListArray));
        orderingArray[i]->vertex = ordering[i]->vertex;
        orderingArray[i]->orderNumber = ordering[i]->orderNumber;
        // if (orderingArray[i]->vertex == 175421) fprintf(stderr, "175421, orderNumber = %d\n", ordering[i]->orderNumber );
        // if (orderingArray[i]->vertex == 573) fprintf(stderr, "573, orderNumber = %d\n", ordering[i]->orderNumber );
        orderingArray[i]->laterDegree = length(ordering[i]->later);
        orderingArray[i]->later = (int *)Calloc(orderingArray[i]->laterDegree, sizeof(int));
       
        int j=0;
        Link* curr = ordering[i]->later->head->next;
        while(!isTail(curr))
        {
            orderingArray[i]->later[j++] = (int)(size_t)(curr->data);
            curr = curr->next;
        }

        orderingArray[i]->earlierDegree = length(ordering[i]->earlier);
        orderingArray[i]->earlier = (int *)Calloc(orderingArray[i]->earlierDegree, sizeof(int));
        
        j=0;
        curr = ordering[i]->earlier->head->next;
        while(!isTail(curr))
        {
            orderingArray[i]->earlier[j++] = (int)(size_t)(curr->data);
            curr = curr->next;
        }
    }

    for(i = 0; i<size;i++)
    {
        // fprintf(stderr, "vertex = %d, orderNumber = %d, laterdeg = %d, earlierdeg = %d\n", orderingArray[i]->vertex, orderingArray[i]->orderNumber, orderingArray[i]->laterDegree, orderingArray[i]->earlierDegree );
        Free(ordering[i]);
        destroyLinkedList(verticesByDegree[i]);
    }
    Free(ordering);
    Free(vertexLocator);
    Free(verticesByDegree);
    Free(degree);
    return orderingArray;
}

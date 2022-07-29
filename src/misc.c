#include<assert.h>
#include<stdio.h>
#include<time.h>
#include<sys/resource.h>
#include<stdlib.h>
#include<string.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"degeneracy_helper.h"

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

using namespace std; // trying to use vectors and stuff 

Network *curr_network ; 
long int m__ ;
int num_nodes ;
Network::Vertex *source_ptr; 

int *reordered__; 
// int *ck ;  
unsigned long long int count_; 

extern bool GRAPH_REDUCE; 
extern bool SAMPLE ;
extern bool BATCH_UPDATE ;

bool SAMPLE_REDUCE = true ; 

long long unsigned int *rho ; 

int *reordered; 
int *level ; 

unsigned long long int *y ; 

int num_iterations = 2;  // this is the total number of iterations
int iter_num;

// sampling related 
int sample_size = 1000000; 
                
extern double sampling_ratio ; 
extern vector<int> ck_sampling;
vector<vector<int>> sampled_kcliques; 
vector<bool> ver_exists_in_approx; 

bool HYBRID = false;
int min_k = 43;

int cliquesize_ ; 
int graph_num_vertices;

double sub_den; 
double last_sub_den ; 
double last_checked_density=0; 
int last_search_scope; 
double fk;
double LB;

extern vector<int> group; 
extern vector<int> uf_rank;
extern vector<int> curr_mc ;
extern vector<bool> MC_lookup; 
extern vector<bool> KCC_lookup; 
extern vector<double> cliq_cnt; 
extern vector<vector<double>> per_ver_cnt ; 
extern vector<int> max_k_per_ver ; 
extern long int total_num_updates ;
extern double yyy ;

vector<Network::Vertex> R;

vector<double> Ck ; 
vector<double> Ck_next; 
vector<int> corenumber; 

// parallel related switches 
bool DISTRIBUTE_DYNAMIC   = false ;                      // looks like this switch is the best
bool DISTRIBUTE_EVENLY    = true ; 
bool PLAIN = false ; 

bool found_KDDS ;

// parallel related variables 
unsigned long long int total_work;
unsigned long long int individual_work;
unsigned long long int max_local; 
unsigned long long int max_load = 0 ;
int num_threads; 
// newary chunks[1000000] ;// the size of chunks is limited 
int chunks_len; 
std::vector< std::vector<int>> distribute_helper; 

double nCr[1001][401];
// 
bool KCC =   false  ; // only turn on for big-w graphs
bool EXACT = false; 
extern bool KCC_COUNT;

static int CDF_NodeCmp(const void *a, const void *b) 
{
  double x = rho[*(const unsigned *)a];
  double y = rho[*(const unsigned *)b];
  if (x > y) return -1;
  if (x < y) return 1;
  return 0;
}

bool compareTreeNode(treeNode* x, treeNode* y){return (x->max_depth > y->max_depth);}

void runAndPrintStatsCliques(  LinkedList** adjListLinked,
                               int n, const char * gname, 
                               char T, int max_k, int flag_d)
{
    //
    fflush(stderr);
    cliquesize_ = max_k;
    graph_num_vertices  = n ; 

    clock_t start = clock();

    // double totalCliques = 0;
    int deg = 0, m = 0;

    // printf("initialize corenumbers. \n");
    corenumber.resize(n); 
    fill(corenumber.begin(), corenumber.end(),0);

    // printf("core decomposition process. \n");
    NeighborListArray** orderingArray = computeDegeneracyOrderArray(adjListLinked, n, max_k);

    // printf("done decomp \n");
    for(int i=0;i<n;i++)
    {
        if(orderingArray[i]==NULL) 
            continue ;
        if (deg < orderingArray[i]->laterDegree) 
            deg = orderingArray[i]->laterDegree;
        m += orderingArray[i]->laterDegree;

        // printf("vertex: %d , orderNumber: %d \n", orderingArray[i]->vertex, orderingArray[i]->orderNumber); // 
    }
    // printf("done listing order numbers \n");

    // if (max_k == 0) max_k = deg + 1; // this is important

    // printf("initialize the root node and pass it in the recursion. \n"); 
    treeNode *root = (treeNode *)Calloc(1, sizeof(treeNode)); 
    root->id = -1 ; 
    root->label = 0; 
    root->max_depth = 0; 
    root->num_cliques = 0; 
    root->parent = NULL; 
    root->pos = -1; 

    // pass the root to build the SCT index. 
    printf("Building the SCT index\n");
    // max depth of each node computed along the way

    // pass in the MC look up table in the SCT build function
    init_MC(n, adjListLinked); 

    int SCT_k = (HYBRID == true) ? min_k : max_k ;  
    listAllCliquesDegeneracy_A(root, orderingArray, n, deg + 1, SCT_k); // max_k is actually current k.

    // in this process, record the vertices that we did not build the tree for but has a legitimate outneighrborhood. 
    clock_t end = clock();
    printf("SCT build time: %lf\n", (double)(end-start)/(double)(CLOCKS_PER_SEC));
    // when building the tree, I track the max reachable depth for each node, which is useful for later pruning. 
    printf("Max depth of tree is %d \n", root->max_depth); 

    // sort the subtrees by their max_depth. 
    // sort(root->child_node.begin(), root->child_node.end(), compareTreeNode); // seems like this is trivial. 

    int MAX_TREE_DEPTH = root->max_depth+1; 

    int size_reordered = 0 ;
    
    double seq_time = 0 ; 

    int inputk = max_k;
    for(int i=0;i<1;i++){

        printf("k=%d \n", max_k);
        if(EXACT) found_KDDS = false;
        // distribute_helper.clear();
        init_work_by_thread();
        double start__=omp_get_wtime(); 

        num_threads = 1 ; 
        sub_den =  nCr[curr_mc.size()][max_k] / curr_mc.size() ; 
        printf("size of maximum clique = %d, MC-kclique-density = %lf  \n", curr_mc.size(), sub_den); 

        KCC = curr_mc.size() >= 100; 
        if(KCC)
        {
            printf("KCC-based analysis is turned on! \n");
            clock_t ttt = clock();
            // long int total____ = ListRootToLeavePath(root, 0, max_k, 0);
            // printf("total____ = %ld, valid num_updates = %ld \n", total____, total_num_updates);
            // exit(0);

            total_work = 0 ;

            int size_of_NB_graph = 0 ; 
            vector<bool> visited;
            visited.resize(n); 
            fill(visited.begin(), visited.end(), false); 

            // vector<int> curr_mc = get_curr_mc();
            // printf("size of maximum clique = %d \n", curr_mc.size());
            int num_cc = 0; 

            // Compute a kclique-isolating partition
            // to begin with, all vertices flag is just vertex id. 
            group.resize(n) ;
            uf_rank.resize(n) ;
            for(int i=0;i<n;i++) group[i] = i ; 
            fill(uf_rank.begin(), uf_rank.end(),1); // initialize all ranks to be 1 

            UFKclique(root, 0, max_k); 
            // during this process, we can compute the kclique density of a preset suboptimal solution.
            
            // only visit the vertices in Gk. 
            vector<int> community_size; 
            community_size.resize(n);
            fill(community_size.begin(), community_size.end(), 0);
            // a similar array can be used to store the maximumn kclique count per vertex in each KCC.

            vector<double> max_ver_cnt ; 
            max_ver_cnt.resize(n);
            fill(max_ver_cnt.begin(), max_ver_cnt.end(), 0);

            KCC_lookup.resize(n);
            fill(KCC_lookup.begin(), KCC_lookup.end(), false);

            for(int i=0; i<n ; i++){
                if(per_ver_cnt[i][max_k]>0) {
                    int group_id = find_func(i); 
                    community_size[ group_id ]++; 
                    // for each group. keep track of the maximum per vetex count 
                    max_ver_cnt[group_id] = max_ver_cnt[group_id] > per_ver_cnt[i][max_k] ? max_ver_cnt[group_id] : per_ver_cnt[i][max_k] ; 
                }
            }

            // We only visit each group once.
            int total_size = 0; 
            int discarded_kcc = 0 ;
            for(int i=0;i<n;i++)
            {
                // for each vertex i , its group id is find_func(i), the size of the group is community_size[find_func(i)].
                // actually, for each vertex i, its group id is find_func(v); 
                // the max_ver_cnt is max_ver_cnt[group id] 
                int group_id = find_func(i); 

                if(community_size[group_id]==0 || visited[group_id]) continue; 

                visited[group_id] = true;

                // printf("KCC %d: size = %d, max cnt = %0.f \n", group_id, community_size[group_id], max_ver_cnt[group_id]); 
                if(max_ver_cnt[group_id]>=sub_den)
                {
                    // if(max_ver_cnt[group_id]>=  max_k*sub_den){ // I think this may be wrong.
                    printf("KCC %d: size = %d, max cnt = %0.f \n", group_id, community_size[group_id], max_ver_cnt[group_id]); 

                    KCC_lookup[group_id] = true; // this should be updated at the end of each iteration.

                    // check if MC is unique
                    if(cliq_cnt[curr_mc.size()]==1){

                        printf("\tIn this graph the MC is unique. \n");
                        
                        if( find_func(curr_mc[0]) == group_id){
                            printf("\tMC resides in this KCC! \n");

                            // if size of MC and this KCC are the same, then this KCC itself is the KDDS. 
                            if(curr_mc.size() == community_size[i]){
                                printf("\tBy theorem 1, the MC is the KDDS of the graph.\n");
                                found_KDDS = true;
                                break;
                            }
                            else{
                                printf("\tCompute max cliq count per vertex outside MC!\n");
                                // in this cases, check the maximum clique count in KCC\MC. 
                                double max_outside_cnt = 0 ;
                                for(int j=0;j<n;j++){
                                    if(find_func(j) == group_id && !MC_lookup[j]){
                                        max_outside_cnt = max_outside_cnt > per_ver_cnt[j][max_k] ? max_outside_cnt : per_ver_cnt[j][max_k]; 
                                    }
                                }
                                if(sub_den > max_outside_cnt){
                                    printf("\tdensity(MC) = %0.f > max cnt out side = %0.f\n", sub_den, max_outside_cnt);
                                    printf("\tBy theorem 1, the MC is the KDDS of the graph\n");
                                    found_KDDS = true;
                                    break;
                                }
        
                                // control the number of iterations this while-loop runs. 
                                vector<bool> cliques_lookup = MC_lookup ; 
                                vector<bool> nb_cliques_lookup = cliques_lookup ; 
                                vector<int> curr_clique = curr_mc;

                                int max_iter_ = 0;
                                while(sub_den <= max_outside_cnt){
                                    // while(max_k*sub_den <= max_outside_cnt){ // can I put a factor k here? 
                                    printf("\tTheorem 1 not satisfied! \n");
                                    printf("\tdensity(MC) = %0.f <= max cnt out side = %0.f\n", sub_den, max_outside_cnt);

                                    // compute NB(curr_clique), update nb_cliques_lookup
                                    for(auto v_ : curr_clique){
                                        Link* curr = adjListLinked[v_]->head ; 
                                        curr = curr->next ;
                                        while(curr != adjListLinked[v_]->tail)
                                        {
                                            int u_ = curr->data ; 
                                            if(!nb_cliques_lookup[u_] && find_func(u_)==group_id )  // only consider the neighbors in Gk. (in the current KCC)
                                                nb_cliques_lookup[u_] = true; 
                                            curr = curr->next ; 
                                        }
                                    }
                                    // find the largest clique not in NB(cliques).

                                    vector<vector<int>> clique_summary; 
                                    clique_summary.resize(curr_mc.size()+1); 

                                    vector<int> large_clique; 
                                    int mk_outside_NB=-1; 
                                    for(int j=0;j<n;j++){
                                        if(!nb_cliques_lookup[j] && find_func(j) == group_id){
                                            clique_summary[max_k_per_ver[j]].push_back(j); 
                                            mk_outside_NB = mk_outside_NB > max_k_per_ver[j] ? mk_outside_NB : max_k_per_ver[j]; 
                                        }
                                    }
                                    // printf("size of 102 ver = %d \n", clique_summary[102].size());
                                    // printf("mk_outside_NB = %d, number of vertices with this cnt = %d \n", mk_outside_NB, large_clique.size());
                                    for(int xx = mk_outside_NB; xx >=3; xx--){
                                        if(clique_summary[xx].size() == xx){
                                            large_clique = clique_summary[xx];
                                            break;
                                        }
                                    }
                                    printf("leaving out %d-clique \n", large_clique.size());
                                    curr_clique = large_clique;
                                    for(auto xx : curr_clique) cliques_lookup[xx] = true ;
                                    max_outside_cnt = 0 ;
                                    for(int j=0;j<n;j++){
                                        if(find_func(j)==group_id && !cliques_lookup[j]){
                                            max_outside_cnt = max_outside_cnt > per_ver_cnt[j][max_k] ? max_outside_cnt : per_ver_cnt[j][max_k]; 
                                        }
                                    }
                                    max_iter_ ++;
                                    if(max_iter_==6) break;
                                }

                                if(sub_den > max_outside_cnt){
                                    printf("\tdensity(MC) = %0.f > max cnt out side = %0.f\n", sub_den, max_outside_cnt);
                                    printf("\tBy theorem 1, the MC is the KDDS of the graph\n");
                                    found_KDDS = true;
                                    break;
                                }else{
                                    printf("need more iterations to leave more cliques out \n");
                                }
                            }
                        }else{
                            printf("MC is not in this KCC\n");
                        }
                    }

                }else{
                    discarded_kcc++;
                }
                num_cc++;
                total_size+=community_size[group_id];
            }

            if(!found_KDDS){
                printf("number of KCC = %d (%d of them pruned by MC-density), total= %d \n", num_cc, discarded_kcc, total_size);
            }
        }
        
        if(found_KDDS) break;

        if((!found_KDDS) && EXACT){
            printf("Phase 1: Sampling 10^7 kcliques. use them to update rho. \n");
            // Obtain a sub-optimal k-clique density. 
            rho = (long long unsigned int *) Calloc(n, sizeof(long long unsigned int));
            for(int i=0;i<n;i++) rho[i] = 0.0;
            reordered = (int *) Calloc(n, sizeof(int));
            level = (int *) Calloc(n, sizeof(int));
            y = (unsigned long long int *) Calloc(n, sizeof(unsigned long long int));

            // sample 10^7 k-cliques 
            sampled_kcliques.clear();
            sampling_ratio = (sample_size / cliq_cnt[max_k]) ; 
            printf("Sampling_ratio = %lf \n",sampling_ratio);
            if(sampling_ratio>=1) 
                sampling_ratio = 1 ; 
            // in the first iteration collect the kclique. In this process, compute the initial values of Ck(v, G'). G' is the sampled graph
            GRAPH_REDUCE = false;
            SAMPLE = true;
            long int cliq_tmp = ListRootToLeavePath(root, 0, max_k, 0);

            printf("cliq_tmp = %ld, sampling_ratio = %lf  \n", cliq_tmp, sampling_ratio);

            // shuffle the sampled kcliques 
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::shuffle(std::begin(sampled_kcliques), std::end(sampled_kcliques), std::default_random_engine(seed));         
            
            // if(sampling_ratio<1){
            printf("cliq_tmp = %ld, sample_size = %d, sampled_kcliques.size() = %d \n", cliq_tmp, sample_size, sampled_kcliques.size());
                // assert(cliq_tmp == sample_size && sample_size ==sampled_kcliques.size());
            

            // go through the sampled cliques for one iteration.
            printf("Go through the sampled cliques for one iteration, sampled_kcliques.size() = %d \n", sampled_kcliques.size());
            for(unsigned ic = 0 ; ic<sampled_kcliques.size(); ic++){
                unsigned node_index = 0;
                for (unsigned i = 1; i < max_k; ++i) {
                    if (rho[ sampled_kcliques[ic][node_index] ] > rho[sampled_kcliques[ic][i]]) node_index = i;
                }
                rho[sampled_kcliques[ic][node_index]]++;
            }

            // recover the sub-optimal density in the original graph. 
            printf("recover the sub-optimal density in the original graph. \n");
            
            bool smart_sampling = true ; 
            
            for (int i = 0; i < n; ++i) 
                reordered[i]  = i; 

            qsort(reordered, n, sizeof(int), CDF_NodeCmp); 

            // #pragma omp parallel for
            for (int i = 0; i < n; ++i){
                level[reordered[i]] = i;
                y[i]=0;
            }

            if(smart_sampling){
                // when sampling, we go through the sampled cliques to get yi.
                for(unsigned ic = 0 ; ic<sampled_kcliques.size(); ic++){
                    unsigned node_index = 0;
                    for (unsigned i = 1; i < max_k; ++i) {
                        if (level[ sampled_kcliques[ic][node_index] ] < level[sampled_kcliques[ic][i]]){
                            node_index = i;
                        }
                    }
                    y[sampled_kcliques[ic][node_index]]++;
                }  
            }
            else{
                yyy=0;
                SAMPLE = false ; 
                compute_yi(root, 0, max_k, level, y); 
                SAMPLE = true ; 
                printf("compute yi updates = %ld \n", yyy);
            }

            unsigned long long int sum = 0;
            num_nodes = 1; 
            double density = -1; 
            for(int i_=0; i_<n; ++i_)// for(int i_=0; i_< n; ++i_)
            {
                sum += y[reordered[i_]];
                if ((double)sum / (i_ + 1) >= density) { // should we use >= or >
                    density = (double)sum / (i_ + 1);
                    num_nodes = i_ + 1; // info.n is the sweet spot
                }
            }  

            if(smart_sampling ){
                printf("Approx density on sampled graph = %lf \n", density);
                printf("In the last iteration of SCTList-sampling, recover the density on the original graph \n"); 
                ver_exists_in_approx.resize(n);
                fill(ver_exists_in_approx.begin(), ver_exists_in_approx.end(),false);
                for(int xx=0; xx< num_nodes ; xx++)
                    ver_exists_in_approx[reordered[xx]]=true;

                yyy=0;
                compute_yi(root, 0, max_k, level, y);
                printf("# k-cliques in approx solution = %lf \n", yyy);
                density = yyy*1.0/num_nodes;
            }       
            printf("%d nodes, k-clique density= %lf \n", num_nodes, density);

            // break;
                    
            sub_den = density; 
            
            printf("Phase 2: Iteratively apply graph reduction until the search scope does not change.\n");
            int search_scope = last_search_scope = 0;
            // initialize Ck to be Ck(v, G).
            Ck.resize(n);
            Ck_next.resize(n);
            for(int i=0;i<n;i++) {
                Ck[i] = per_ver_cnt[i][max_k]; // we need this 
                if(Ck[i]>= ceil(sub_den)) 
                    search_scope++;
            }
            printf("Search scope (vertices with per ver cnt >= ceil(MC-density) = %d \n", search_scope);
            GRAPH_REDUCE = true; 
            SAMPLE = false; 
            BATCH_UPDATE = false;
            while(last_search_scope != search_scope){
                fill(Ck_next.begin(),Ck_next.end(),0);
                ListRootToLeavePath(root, 0, max_k, 0);
                Ck = Ck_next ;
                last_search_scope = search_scope; // search scopes updating
                search_scope = 0 ; 
                for(int xxx = 0 ;xxx<n ; xxx++)
                {
                    if(Ck[xxx]>=ceil(sub_den) ) {
                        search_scope++;
                    }
                    else{ 
                        rho[xxx] = 0 ; // setting rho to zero
                    }
                }
                printf("search scope = %d \n", search_scope);
            }

            // reset all rho values. // this is stupid 
            for(int xx= 0;xx<n;xx++) rho[xx] = 0; 
            
            // check goldberg.

            // Phase 3: Doubling the number of iterations until DDS is found. 
            int T__ = 10; 
            int num_iter_run = 0 ;
            while(!found_KDDS){

                // run approximate for T_ iterations. 
                printf("running %d iterations\n", T__); 
                for(iter_num=0; iter_num < T__; iter_num++){
                    if( (iter_num+1) % (int)(T__/10)==0) printf("\t itr %d \n", iter_num+1); 
                    kcliquePlusPlus(n, rho, root, 0, max_k);
                }
                // obtain approximate solution, inheriting the rho values. 
                if(1==1){
                    // sort vertices by rho values 
                    for (int i = 0; i < n; ++i) 
                        reordered[i]  = i; 

                    qsort(reordered, n, sizeof(int), CDF_NodeCmp); 

                    for (int i = 0; i < n; ++i){
                        level[reordered[i]] = i;
                        y[i]=0;
                    }

                    yyy=0;
                    compute_yi(root, 0, max_k, level, y); 
                    printf("compute yi updates = %ld \n", yyy);  

                    // Find the approximate maximum density
                    unsigned long long int sum = 0;
                    num_nodes = 1; 
                    double density = -1; 
                    for(int i_=0; i_<(GRAPH_REDUCE ? search_scope : n); ++i_)// for(int i_=0; i_< n; ++i_)
                    {
                        sum += y[reordered[i_]];
                        if ((double)sum / (i_ + 1) >= density) { // should we use >= or >
                            density = (double)sum / (i_ + 1);
                            num_nodes = i_ + 1; // info.n is the sweet spot
                        }
                    }
                    printf("approx: %d nodes, k-clique density= %lf\n",num_nodes, density);
                    // printf("sol: ");
                    // for(int i_=0; i_<num_nodes; i_++ ){
                    //     printf("%d ", reordered[i_]);
                    // }
                    // printf("\n");

                    last_sub_den = sub_den;
                    sub_den = sub_den > density ? sub_den : density ;
                    // check optimality.
                    bool check_optimal =  (last_checked_density!=density) && last_sub_den==sub_den ;
                    if(check_optimal)
                    {
                        // Checking 1 // the rho values should be density vectors.
                        num_iter_run+=T__; 
                        if(goldberg_check_optimal(density,num_iter_run)) {
                            printf("The current solution is optimal by Goldberg's condition! (tighter UB)\n"); 
                            printf("%d nodes, k-clique density= %lf \n",num_nodes, density);
                            found_KDDS=true;
                            break;
                        }
                        else{
                            printf("Not optimal by Goldberg's condition.\n"); 
                        }
                        // Checking 2         
                        // when faced with the same solution in different order, the maxflow returns differenet results.                
                        if(maxflow_check_optimal(root)){
                            printf("The current solution is optimal by Maxflow!\n"); 
                            printf("%d nodes, k-clique density= %lf \n",num_nodes, density);
                            found_KDDS=true;
                            break;
                        }else{
                            // printf("The current solution is not optimal...\n"); 
                            printf("Not optimal by Maxflow\n"); 
                            // if(density>8065.9) 
                            //     exit(0);
                        }
                        
                        last_checked_density = density; 
                        /*
                        if(search_scope == num_nodes){
                            printf("When current search scope and the approx solution are the same, quickly check if it is KDDS of itself. \n");
                            // best density is density. 
                            vector<double> Ck_sorted = Ck ;
                            sort(Ck_sorted.begin(), Ck_sorted.end(), greater<double>());

                            // for(auto ele:Ck_sorted) if(ele>0) printf("CK: %lf \n", ele);

                            bool dominating = true; 
                            
                            double sum_ck=0;

                            for(int xxx=0;xxx<curr_mc.size();xxx++)
                                sum_ck += Ck_sorted[xxx];
                            
                            for(int size__ = curr_mc.size()+1; size__ < num_nodes; size__++){
                                // compute its upper bound. 
                                double ub__;

                                if(size__ < num_nodes-1){
                                    sum_ck +=  Ck_sorted[size__-1];
                                    ub__ = sum_ck/max_k/size__;
                                }else{
                                    ub__ = ( density*num_nodes - Ck_sorted[num_nodes-1]  )/size__;
                                }
                                if(density < ub__){
                                    printf("\t\tsize = %d, density = %lf, ub = %lf \n", size__, density, ub__ );
                                    dominating = false; 
                                    break;
                                }else{
                                    printf("(dominated) size = %d, density = %lf, ub = %lf \n", size__, density, ub__ );
                                }
                            }
                            if(dominating) {
                                printf("The current solution is optimal by Jerry's condition! (tighter UB)\n"); 
                                printf("%d nodes, k-clique density= %lf \n",num_nodes, density);
                                found_KDDS=true;
                                break;
                            }
                            else{
                                printf("Not optimal by Jerry's condition.\n"); 
                                // break;
                            }          

                        }
                        */
                    } 
                    // 
                }
                T__ = (int)(T__* 1.1);
            }

            Free(rho);
            Free(reordered);
            Free(level);
            Free(y);
        }
        // (sampling) approximate algorithm
        if(!EXACT)
        {
            // these are only needed when graph reduction is switched on. 
            int search_scope ;
            if(GRAPH_REDUCE )
            {
                Ck.resize(n);
                Ck_next.resize(n);
                search_scope=0;
                for(int i=0;i<n;i++) {
                    Ck[i] = per_ver_cnt[i][max_k];
                    if(Ck[i]>= ceil(sub_den)) 
                        search_scope++;
                }
                // need to compute the number of edges 
                // int num_edges =0 ;
                // for(int v_ = 0; v_ < n; v_++){
                //     Link* curr = adjListLinked[v_]->head ; 
                //     curr = curr->next ;
                //     // printf("xxx:%d", v_); 
                //     if(Ck[v_]>= ceil(sub_den)) {
                //         while(curr != adjListLinked[v_]->tail)
                //         {
                //             if(Ck[curr->data]>= ceil(sub_den)){
                //                 num_edges++;
                //             }
                //             curr = curr->next ; 
                //         }
                //     }
                // }
                // printf("num edges = %d \n", num_edges);
                printf("Search scope (vertices with per ver cnt >= ceil(MC-density) = %d \n", search_scope);
                
            }

            if(SAMPLE )
            {
                GRAPH_REDUCE = false;
                Ck.resize(n);
                fill(Ck.begin(), Ck.end(),0);
                // num_iterations = 2;
            }
        
            rho = (long long unsigned int *) Calloc(n, sizeof(long long unsigned int));
            for(int i=0;i<n;i++) rho[i] = 0.0;
            reordered = (int *) Calloc(n, sizeof(int));
            level = (int *) Calloc(n, sizeof(int));
            y = (unsigned long long int *) Calloc(n, sizeof(unsigned long long int));

            // double start__=omp_get_wtime(); 
            for(iter_num=0; iter_num < num_iterations; iter_num++)
            {   
                printf("\nrunning %d iteration: \n", iter_num+1); 

                init_work_by_thread();
        
                // when the sampling switch is on, we simply pick and choose which subtrees we explore. 
                if(GRAPH_REDUCE){
                    fill(Ck_next.begin(),Ck_next.end(),0); // recompute this.
                }
                    
                // We sample once and run for many iterations
                if(SAMPLE){
                    // in the first iteration, collect sigma k-cliques. 
                    if(iter_num == 0 ){                    
                        sampled_kcliques.clear();

                        sampling_ratio = sample_size / cliq_cnt[max_k] ; 

                        printf("max_k = %d, Sampling_ratio = %lf, cliq_cnt[max_k] = %lf \n",
                            max_k, sampling_ratio, cliq_cnt[max_k]);

                        if(sampling_ratio>=1) sampling_ratio = 1 ; 
                        
                        // ideally we want to increase the sampling ratio a bit. 

                        // in the first iteration collect the kclique. In this process, compute the initial values of Ck(v, G'). G' is the sampled graph
                        long int cliq_tmp = ListRootToLeavePath(root, 0, max_k, 0);

                        // shuffle the sampled kcliques 
                        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                        std::shuffle(std::begin(sampled_kcliques), std::end(sampled_kcliques), std::default_random_engine(seed));
                                            

                        printf("cliq_tmp = %d \n",cliq_tmp);
                        printf("sample_size = %d \n",sample_size);
                        printf("sampled_kcliques.size() = %d \n", sampled_kcliques.size());

                        // if(sampling_ratio<1){
                        //     assert(cliq_tmp == sample_size && sample_size ==sampled_kcliques.size());
                        // }
                        

                        for(unsigned ic = 0 ; ic<sampled_kcliques.size(); ic++){
                            unsigned node_index = 0;
                            
                            Ck[sampled_kcliques[ic][node_index]]++;
  
                            // go through the vertices of each \kclique
                            for (unsigned i = 1; i < max_k; ++i) {
                                if (rho[ sampled_kcliques[ic][node_index] ] > rho[sampled_kcliques[ic][i]]) node_index = i;
                                Ck[sampled_kcliques[ic][i]]++;
                            }
                            rho[sampled_kcliques[ic][node_index]]++;
                            total_num_updates++;
                        }
                    }
                    else if(iter_num == 1 ){ 

                        // fill(Ck_next.begin(), Ck_next.end(),0);
                        vector<vector<int>> reduced_sampled_cliques;

                        for(unsigned ic = 0 ; ic<sampled_kcliques.size(); ic++){

                            unsigned node_index = 0;
                            bool dominated = false;
                            if(SAMPLE_REDUCE) dominated = Ck[sampled_kcliques[ic][node_index]] < ceil(sub_den);

                            // Ck_next[sampled_kcliques[ic][node_index]]++;
                            for (unsigned i = 1; i < max_k; ++i) {
                                if (rho[ sampled_kcliques[ic][node_index] ] > rho[sampled_kcliques[ic][i]]) 
                                    node_index = i;
                                if(SAMPLE_REDUCE && Ck[sampled_kcliques[ic][i]] < ceil(sub_den)) 
                                    dominated = true;
                                // Ck_next[sampled_kcliques[ic][i]]++;
                            }
                            if(!dominated){
                                rho[sampled_kcliques[ic][node_index]]++;
                                total_num_updates++;
                                reduced_sampled_cliques.push_back(sampled_kcliques[ic]);
                            }
                        }    
                        sampled_kcliques = reduced_sampled_cliques;
                        // Ck = Ck_next;    
                    }
                    else{
                        // shuffle the sampled kcliques more often is not useful. 
                        for(unsigned ic = 0 ; ic<sampled_kcliques.size(); ic++){
                            // only proceed when all vertices satisfy Ck[ver] >= ceil(sub_den); 

                            unsigned node_index = 0;
                            for (unsigned i = 1; i < max_k; ++i) {
                                if (rho[ sampled_kcliques[ic][node_index] ] > rho[sampled_kcliques[ic][i]]) node_index = i;
                            }
                            rho[sampled_kcliques[ic][node_index]]++;
                            total_num_updates++;
                        }              
                    }
                    // in each iteration, visit all sampled kcliques.    
                }else{
                    kcliquePlusPlus(n, rho, root, 0, max_k);
                }
                

                if(GRAPH_REDUCE){
                    Ck = Ck_next ; 
                }

                // double xxx2=omp_get_wtime(); 


                // if sampling and in the first iteration. compute the approximate denstiy
                if(SAMPLE && iter_num == 0){
                    for (int i = 0; i < n; ++i) 
                        reordered[i]  = i;    

                    qsort(reordered, n, sizeof(int), CDF_NodeCmp); 

                    for (int i = 0; i < n; ++i){
                        level[reordered[i]] = i;
                        y[i]=0;
                    }
                    // compute yi;
                    for(unsigned ic = 0 ; ic<sampled_kcliques.size(); ic++){
                        unsigned node_index = 0;
                        for (unsigned i = 1; i < max_k; ++i) {
                            if (level[ sampled_kcliques[ic][node_index] ] < level[sampled_kcliques[ic][i]]){
                                node_index = i;
                            }
                        }
                        y[sampled_kcliques[ic][node_index]]++;
                    }
                    // find the best solution under this ordering // vertices 0 <= i < num_nodes. ver = reordered. 
                    unsigned long long int sum = 0;
                    int num_nodes = 1; 
                    double density = -1; 
                    for (int i = 0; i < n; ++i) 
                    {
                        sum += y[reordered[i]];
                        if ((double)sum / (i + 1) >= density) { // should we use >= or >
                            density = (double)sum / (i + 1);
                            num_nodes = i + 1; 
                        }
                    }
                    printf("Approx density on sampled graph = %lf \n", density);
                    sub_den = density; 
                }

                // this dictates how often we obtain the sub-optimal density
                if( (GRAPH_REDUCE &&  search_scope> 100 && (iter_num+1)%(num_iterations/2) ==1  ) || iter_num+1 == num_iterations  ) 
                { 
                    bool smart_sampling = true ; 

                    if(GRAPH_REDUCE){
                        printf("current search scope larger than 100\n");
                        last_search_scope = search_scope;
                        search_scope = 0 ; 
                        printf("sub_den = %lf \n",sub_den);
                        for(int xxx = 0 ;xxx<n ; xxx++)
                        {
                            if(Ck[xxx]>=ceil(sub_den) ) {
                                search_scope++;
                            }
                            // else{
                            //     rho[xxx] = 0 ;
                            // }
                        }
                    }
                    // sort vertices by rho values 
                    for (int i = 0; i < n; ++i) 
                        reordered[i]  = i; 

                    qsort(reordered, n, sizeof(int), CDF_NodeCmp); 

                    // #pragma omp parallel for
                    for (int i = 0; i < n; ++i){
                        level[reordered[i]] = i;
                        y[i]=0;
                    }
                    
                    if(smart_sampling && SAMPLE){
                        // when sampling, we go through the sampled cliques to get yi.
                        for(unsigned ic = 0 ; ic<sampled_kcliques.size(); ic++){
                            unsigned node_index = 0;
                            for (unsigned i = 1; i < max_k; ++i) {
                                if (level[ sampled_kcliques[ic][node_index] ] < level[sampled_kcliques[ic][i]]){
                                    node_index = i;
                                }
                            }
                            y[sampled_kcliques[ic][node_index]]++;
                        }  
                    }else{
                        yyy=0;
                        bool old_label = SAMPLE; 
                        SAMPLE = false;
                        compute_yi(root, 0, max_k, level, y); 
                        SAMPLE = old_label; 
                        printf("compute yi updates = %ld \n", yyy);
                    }

                    // Find the approximate maximum density
                    unsigned long long int sum = 0;
                    num_nodes = 1; 
                    double density = -1; 
                    for(int i_=0; i_<(GRAPH_REDUCE ? search_scope : n); ++i_)// for(int i_=0; i_< n; ++i_)
                    {
                        sum += y[reordered[i_]];
                        if ((double)sum / (i_ + 1) >= density) { // should we use >= or >
                            density = (double)sum / (i_ + 1);
                            num_nodes = i_ + 1; // info.n is the sweet spot
                        }
                    }

                    if(smart_sampling && SAMPLE){
                        printf("Approx density on sampled graph = %lf \n", density);
                        printf("In the last iteration of SCTList-sampling, recover the density on the original graph \n"); 
                        ver_exists_in_approx.resize(n);
                        fill(ver_exists_in_approx.begin(), ver_exists_in_approx.end(),false);
                        for(int xx=0; xx< num_nodes ; xx++)
                            ver_exists_in_approx[reordered[xx]]=true;

                        yyy=0;
                        compute_yi(root, 0, max_k, level, y);
                        printf("# k-cliques in approx solution = %lf \n", yyy);
                        density = yyy*1.0/num_nodes;
                    }

                    last_sub_den = sub_den;
                    // update sub_den
                    sub_den = sub_den > density ? sub_den : density ;

                    // Compute the upper bound
                    double upperbound = 0; 
                    sum=0;
                    int T = iter_num+1; 
                    double ip1ck = 1; // (i + 1) choose k
                    // for (unsigned i = 0; i < n; ++i) {
                    for(int i_=0; i_<(GRAPH_REDUCE ? search_scope : n); ++i_){
                    // for(int i_=0; i_<n; ++i_){
                        sum += rho[reordered[i_]];
                        if (i_ + 1 == max_k)
                            ip1ck = 1;
                        if (i_ + 1 > max_k)
                            ip1ck = (ip1ck * (i_ + 1)) / (i_ + 1 - max_k);

                        if (ip1ck < (double)sum / T)
                            upperbound = ip1ck / (i_ + 1);
                        else {
                            if (upperbound < (double)sum / T / (i_ + 1)){
                                upperbound = (double)sum / T / (i_ + 1);
                            }
                            break;
                        }
                    }
                    if(SAMPLE) 
                        upperbound = upperbound / sampling_ratio / (1 - sqrt(6 * log(n) / upperbound));
                    // use criteron B 
                    // double num_cliq = pow(num_nodes/(max_k-1), max_k-2);
                    // LB = num_cliq / num_nodes;
                    // assert(density<=upperbound);
                    printf("%d nodes, k-clique density= %lf, UB = %lf, error = %lf\n",num_nodes, density, upperbound, (upperbound-density)/density);
                }
                if(GRAPH_REDUCE){
                    printf("search scope = %d \n", search_scope);
                    // int num_edges =0 ;
                    // for(int v_ = 0; v_ < n; v_++){
                    //     Link* curr = adjListLinked[v_]->head ; 
                    //     curr = curr->next ;
                    //     if(Ck[v_]>= ceil(sub_den)) {
                    //         while(curr != adjListLinked[v_]->tail)
                    //         {
                    //             if(Ck[curr->data]>= ceil(sub_den)) num_edges++;
                    //             curr = curr->next ; 
                    //         }
                    //     }
                    // }
                    // printf("num edges = %d \n", num_edges);
                }
                    
                printNumUpdates(NULL, 0, max_k); 

                // for(int xx =0; xx < n ; xx ++) printf("v: %d,  weight = %llu \n", xx, rho[xx]);
                
            }
            Free(rho);
            Free(reordered);
            Free(level);
            Free(y);
        }
        // printNumUpdates(NULL, 0, max_k); 
        double end__ = omp_get_wtime(); 
        printf("SCTList took %f seconds\n", end__ - start__);
    }

    // Free(reordered__) ; 
    // this function shows how many rho updates are conducted. 

    clock_t end3 = clock();
    printf("KCList++ time: %lf\n", (double)(end3-end)/(double)(CLOCKS_PER_SEC));
    printf("Total time: %lf\n", (double)(end3-start)/(double)(CLOCKS_PER_SEC));

    destroyTreeFromRoot(root); // this function is problematic
    root = NULL ; 
    // show the clique counts: 
    // Free(cliqueCounts);
    // 
    for(int i = 0; i<n; i++)
    {
        Free(orderingArray[i]->later);
        Free(orderingArray[i]->earlier);
        Free(orderingArray[i]);
    }
    Free(orderingArray);
}

void populate_nCr()
{
    // printf("populate_nCr is called \n"); 
    FILE *infile;
    infile = fopen("./src/nCr.txt","r");
    double d=0;
    if(infile==NULL)
    {
        printf("file could not be opened\n");
        exit(1);
    }
    for(int row = 0; row < 1001; ++row)
    {
        for (int col = 0; col < 401; ++col)
        {
            if (!fscanf(infile,"%lf,",&d)) 
                fprintf(stderr, "Error\n");
            // fprintf(stderr, "%lf\n", d);
            nCr[row][col] = d;
        }
    }
    fclose(infile);
}

int nodeComparator(int node1, int node2)
{
    if ((int)(size_t)node1 < (int)(size_t)node2)
        return -1;
    if((int)(size_t)node1 > (int)(size_t)node2)
        return 1;

    return 0;
}

int sortComparator(int node1, int node2)
{
    if (*(int*)node1 < *(int*)node2)
        return -1;
    if(*(int*)node1 > *(int*)node2)
        return 1;

    return 0;
}

int qsortComparator(const void * node1, const void * node2)
{
    return ( *(int*)node1 - *(int*)node2 );
}

void printArray(int* array, int size)
{
    int i = 0;
    while(i<size)
        printf("%d ", array[i++]);
    printf("\n");
}

void printArrayOfLinkedLists(LinkedList** listOfLists, int size)
{
    // list graph contents

    int i=0;

    while(i<size)
    {
        if(!isEmpty(listOfLists[i]))
        {
            printf("%d:", i);
            // printListAbbv(listOfLists[i], &printInt);
        }
        i++;
    }
}

void printClique(int* clique)
{
    int i = 0;
    while(clique[i]!=-1)
    {
        printf("%d", clique[i]);
        if(clique[i+1]!=-1)
            printf(" ");
        i++;
    }
    printf("\n");
}


void printInt(int integer)
{
    printf("%d", (int)(size_t)integer);
}

void destroyCliqueResults(LinkedList* cliques)
{
    Link* curr = cliques->head->next;
    while(!isTail(curr))
    {
        int* clique = (int*)curr->data;

#ifdef DEBUG
        int i=0;
        while(clique[i] != -1)
        {
            printf("%d", clique[i]);
            if(clique[i+1] != -1)
                printf(" ");
            i++;
        }
        printf("\n");
#endif
        Free(clique);
        curr = curr->next;
    } 

    destroyLinkedList(cliques); 
}


LinkedList** readInGraphAdjList(int* n, int* m)
{
    int u, v; // endvertices, to read edges.

    if(scanf("%d", n)!=1)
    {
        fprintf(stderr, "problem with line 1 in input file\n");
        exit(1);
    }

    if(scanf("%d", m)!=1)
    {
        fprintf(stderr, "problem with line 2 in input file\n");
        exit(1);
    }

#ifdef DEBUG
    printf("Number of vertices: %d\n", *n);
    printf("Number of edges: %d\n", *m);
#endif
    
    LinkedList** adjList = (LinkedList**)Calloc(*n, sizeof(LinkedList*));

    int i = 0;
    while(i < *n)
    {
        adjList[i] = createLinkedList();
        i++;
    }

    i = 0;

    while(i < *m)
    {
        if(scanf("%d,%d", &u, &v)!=2)
        {
            printf("problem with line %d in input file\n", i+2);
            exit(1);
        }
        assert(u < *n && u > -1);
        assert(v < *n && v > -1);
        if(u==v)
            printf("%d=%d\n", u, v);
        assert(u != v);

        addLast(adjList[u], (int)v);
        // addLast(adjList[v], (int)u);

        i++;
    }

#ifdef DEBUG
    printArrayOfLinkedLists(adjList, *n);
#endif

    return adjList;
}


LinkedList** readInGraphAdjListToDoubleEdges(int* n, int* m, char *fpath)
{
    int u, v; // endvertices, to read edges.

    FILE *fp;
    fp = fopen (fpath,"r");
    if (!fp) 
    {
        fprintf(stderr, "Could not open input file.\n");
        exit(1);
    }

    if(fscanf(fp, "%d %d", n, m)!=2)
    {
        fprintf(stderr, "Number of vertices: %d\n", *n);
        fprintf(stderr, "Number of edges: %d\n", *m);
        fprintf(stderr, "problem with line 1 in input file\n");
        exit(1);
    }

    LinkedList** adjList = (LinkedList**)Calloc(*n, sizeof(LinkedList*));

    int i = 0;
    while(i < *n)
    {
        adjList[i] = createLinkedList();
        i++;
    }

    i = 0;
    // double maxv = 0;
    // int sizeofgraph=0;

    while(i < *m)
    {
        if (fscanf(fp, "%d %d\n", &u, &v)!=2)
        {
            printf("problem with line %d in input file, u=%d, v=%d\n", i+2, u, v);
            exit(1);
        }

        // only record this edge if u < v

        // if ((u>= *n) || (v >= *n)) printf("u = %d, v = %d \n", u, v);
        // if ((double) u > maxv) maxv = (double) u;
        // if ((double) v > maxv) maxv = (double) v;
        assert(u < *n && u > -1);

        if(v>=*n){
            printf("v = %d , n = %d\n", v, *n);
            exit(1);
        }
        
        assert(v < *n && v > -1);
        
        assert(u != v);
        if(u<v){ // this is for handling directed graphs.
            addLast(adjList[u], (int)v);
            addLast(adjList[v], (int)u);
            // printf("e:%d %d\n", u, v);
        }
        // printf("%d,%d\n", u, v);
        // else{
        //     printf("e:%d %d\n", v, u);
        // }
        i++;
    }
    /*
    int visited =0 ;
    // visit the adjacency list: 
    printf("xxx:%d\n", *n); 
    for(int v_ = 0; v_ < *n; v_++){
        Link* curr = adjList[v_]->head ; 
        curr = curr->next ;
        printf("xxx:%d", v_); 
        while(curr != adjList[v_]->tail)
        {
            // printf("%d %d\n", v_, curr->data);
            printf(" %d",curr->data);
            // visited++;
            // if(visited==100) exit(0);
            curr = curr->next ; 
        }
        visited++;
        // if(visited==100) exit(0);
        printf("\n"); 
    }
    exit(0);
    */

    *m = (*m) * 2;

    fclose(fp);
    return adjList;
}

int findBestPivotNonNeighborsDegeneracyCliques( int** pivotNonNeighbors, int* numNonNeighbors,
                                                int* vertexSets, int* vertexLookup,
                                                int** neighborsInP, int* numNeighbors,
                                                int beginX, int beginP, int beginR)
{
    // printf("check: %d %d %d \n", beginX, beginP, beginR);
    int pivot = -1;
    int maxIntersectionSize = -1;

    // this is vital.
    // iterate over each vertex in P 
    int j = beginP; 

    // iterate over each vertex in X union P
    // int j = beginX;

    while(j<beginR)
    {
        int vertex = vertexSets[j];
        int numPotentialNeighbors = min(beginR - beginP, numNeighbors[vertex]); //bug resolved by Shweta
        // printf("vertex = %d, numPotentialNeighbors  =%d \n", vertex , numPotentialNeighbors);
        int numNeighborsInP = 0;

        int k = 0;
        while(k<numPotentialNeighbors)
        {
            int neighbor = neighborsInP[vertex][k];
            int neighborLocation = vertexLookup[neighbor];

            if(neighborLocation >= beginP && neighborLocation < beginR)
            {
                numNeighborsInP++;
            }
            else
            {
                break;
            }
            k++;
        }

        if(numNeighborsInP > maxIntersectionSize)
        {
            pivot = vertex;
            
            // if(j<beginP){// check if vertex is coming from X, then this branch is not maximal.

            // }

            maxIntersectionSize = numNeighborsInP;
        }

        j++;
    }
    // printf("pivot = %d \n", pivot);


    // compute non neighbors of pivot by marking its neighbors
    // and moving non-marked vertices into pivotNonNeighbors.
    // we must do this because this is an efficient way
    // to compute non-neighbors of a vertex in 
    // an adjacency list.

    // we initialize enough space for all of P; this is
    // slightly space inefficient, but it results in faster
    // computation of non-neighbors.
    *pivotNonNeighbors = (int *)Calloc(beginR-beginP, sizeof(int));
    memcpy(*pivotNonNeighbors, &vertexSets[beginP], (beginR-beginP)*sizeof(int));

    // we will decrement numNonNeighbors as we find neighbors
    *numNonNeighbors = beginR-beginP;

    int numPivotNeighbors = min(beginR - beginP, numNeighbors[pivot]); //bug resolved by Shweta
  
    // printf("numPivotNeighbors = %d \n", numPivotNeighbors); 
    // mark the neighbors of pivot that are in P.
    j = 0;
    while(j<numPivotNeighbors)
    {
        int neighbor = neighborsInP[pivot][j];
        int neighborLocation = vertexLookup[neighbor];

        if(neighborLocation >= beginP && neighborLocation < beginR)
        {
            (*pivotNonNeighbors)[neighborLocation-beginP] = -1;
        }
        else
        {
            break;
        }

        j++;
    }

    // move non-neighbors of pivot in P to the beginning of
    // pivotNonNeighbors and set numNonNeighbors appriopriately.
    // printf("*numNonNeighbors = %d \n", *numNonNeighbors); 
    // if a vertex is marked as a neighbor, the we move it
    // to the end of pivotNonNeighbors and decrement numNonNeighbors.
    j = 0;
    while(j<*numNonNeighbors)
    {
        int vertex = (*pivotNonNeighbors)[j];

        if(vertex == -1)
        {
            (*numNonNeighbors)--;
            (*pivotNonNeighbors)[j] = (*pivotNonNeighbors)[*numNonNeighbors];
            continue;
        }

        j++;
    }
    // printf("numNonNeighbors = %d \n", *numNonNeighbors);
    return pivot; 
}


void fillInPandXForRecursiveCallDegeneracyCliques( int vertex, int orderNumber,
                                                   int* vertexSets, int* vertexLookup, 
                                                   NeighborListArray** orderingArray,
                                                   int** neighborsInP, int* numNeighbors,
                                                   int* pBeginX, int *pBeginP, int *pBeginR, 
                                                   int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{
    // not sure if this is correct
    if(orderingArray[orderNumber]==NULL) return;

    int vertexLocation = vertexLookup[vertex];

    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;
    
    *pNewBeginR = *pBeginR;
    *pNewBeginP = *pBeginR;

    //printf("Before 1st while\n");
    // swap later neighbors of vertex into P section of vertexSets
    int j = 0;
    while(j < orderingArray[orderNumber]->laterDegree)
    {
        int neighbor = orderingArray[orderNumber]->later[j];
        int neighborLocation = vertexLookup[neighbor];

        (*pNewBeginP)--;

        vertexSets[neighborLocation] = vertexSets[*pNewBeginP];
        vertexLookup[vertexSets[*pNewBeginP]] = neighborLocation;
        vertexSets[*pNewBeginP] = neighbor;
        vertexLookup[neighbor] = *pNewBeginP;

        j++; 
    }

    *pNewBeginX = *pNewBeginP;

// add sth here:
    /*
    // swap earlier neighbors of vertex into X section of vertexSets
    j = 0;
    while(j<orderingArray[orderNumber]->earlierDegree)
    {
        int neighbor = orderingArray[orderNumber]->earlier[j];
        int neighborLocation = vertexLookup[neighbor];

        (*pNewBeginX)--;
        vertexSets[neighborLocation] = vertexSets[*pNewBeginX];
        vertexLookup[vertexSets[*pNewBeginX]] = neighborLocation;
        vertexSets[*pNewBeginX] = neighbor;
        vertexLookup[neighbor] = *pNewBeginX;

        Free(neighborsInP[neighbor]);
        neighborsInP[neighbor] = (int*)Calloc(min(*pNewBeginR-*pNewBeginP,orderingArray[neighbor]->laterDegree), sizeof(int));
        numNeighbors[neighbor] = 0;

        // fill in NeighborsInP
        int k = 0;
        while(k<orderingArray[neighbor]->laterDegree)
        {
            int laterNeighbor = orderingArray[neighbor]->later[k];
            int laterNeighborLocation = vertexLookup[laterNeighbor];
            if(laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR)
            {
                neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                numNeighbors[neighbor]++;
            }

            k++;
        }

        j++; 

    }
    */
// 

    // reset numNeighbors and neighborsInP for this vertex
    j = *pNewBeginP;
    //printf("Before 2nd while\n");
    while(j<*pNewBeginR)
    {
        int vertexInP = vertexSets[j];
        //printf("vertexInP = %d, numNeighbors[vertexInP]=%d\n", vertexInP, numNeighbors[vertexInP] );
        //printf("Address being freed: %p\n", neighborsInP[vertexInP]);
        numNeighbors[vertexInP] = 0;
        Free(neighborsInP[vertexInP]);
        //printf("Allocating %d space for neighborsInP[vertexInP].\n", min( *pNewBeginR-*pNewBeginP, 
                                           //  orderingArray[vertexInP]->laterDegree 
                                           //+ orderingArray[vertexInP]->earlierDegree));
        neighborsInP[vertexInP]= (int *)Calloc( min( *pNewBeginR-*pNewBeginP, 
                                             orderingArray[vertexInP]->laterDegree 
                                           + orderingArray[vertexInP]->earlierDegree), sizeof(int));

        j++;
    }

    // count neighbors in P, and fill in array of neighbors
    // in P
    j = *pNewBeginP;
    //printf("Before 3rd while\n");
    while(j<*pNewBeginR)
    {
        int vertexInP = vertexSets[j];

        int k = 0;
        while(k<orderingArray[vertexInP]->laterDegree)
        {
            int laterNeighbor = orderingArray[vertexInP]->later[k];
            int laterNeighborLocation = vertexLookup[laterNeighbor];

            if(laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR)
            {
                neighborsInP[vertexInP][numNeighbors[vertexInP]] = laterNeighbor;
                numNeighbors[vertexInP]++;
                neighborsInP[laterNeighbor][numNeighbors[laterNeighbor]] = vertexInP;
                numNeighbors[laterNeighbor]++;
            }

            k++;
        }

        j++;
    }
}


void moveToRDegeneracyCliques( int vertex, 
                               int* vertexSets, int* vertexLookup, 
                               int** neighborsInP, int* numNeighbors,
                               int* pBeginX, int *pBeginP, int *pBeginR, 
                               int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{

    int vertexLocation = vertexLookup[vertex];
    
    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;

    // this is not a typo, initially newX is empty
    *pNewBeginX = *pBeginP;
    *pNewBeginP = *pBeginP;
    *pNewBeginR = *pBeginP;

    int sizeOfP = *pBeginR - *pBeginP;

    int j = (*pBeginP);
    while(j<(*pBeginR))
    {
        int neighbor = vertexSets[j];
        int neighborLocation = j;

        int numPotentialNeighbors = min(sizeOfP, numNeighbors[neighbor]); 
        int k = 0;
        while(k<numPotentialNeighbors)
        {
            if(neighborsInP[neighbor][k] == vertex)
            {
                vertexSets[neighborLocation] = vertexSets[(*pNewBeginR)];
                vertexLookup[vertexSets[(*pNewBeginR)]] = neighborLocation;
                vertexSets[(*pNewBeginR)] = neighbor;
                vertexLookup[neighbor] = (*pNewBeginR);
                (*pNewBeginR)++;
            }

            k++;
        }

        j++;
    }

    j = (*pNewBeginP);

    while(j < *pNewBeginR)
    {
        int thisVertex = vertexSets[j];

        int numPotentialNeighbors = min(sizeOfP, numNeighbors[thisVertex]); 

        int numNeighborsInP = 0;

        int k = 0;
        while(k < numPotentialNeighbors)
        {
            int neighbor = neighborsInP[thisVertex][k];
            int neighborLocation = vertexLookup[neighbor];
            if(neighborLocation >= *pNewBeginP && neighborLocation < *pNewBeginR)
            {
                neighborsInP[thisVertex][k] = neighborsInP[thisVertex][numNeighborsInP];
                neighborsInP[thisVertex][numNeighborsInP] = neighbor;
                numNeighborsInP++;
            }
            k++;
        }

        j++;
    }
}

void moveFromRToXDegeneracyCliques( int vertex, 
                                    int* vertexSets, int* vertexLookup, 
                                    int* pBeginX, int* pBeginP, int* pBeginR )
{
    int vertexLocation = vertexLookup[vertex];

    //swap vertex into X and increment beginP and beginR
    vertexSets[vertexLocation] = vertexSets[*pBeginP];
    vertexLookup[vertexSets[*pBeginP]] = vertexLocation;
    vertexSets[*pBeginP] = vertex;
    vertexLookup[vertex] = *pBeginP;

    *pBeginP = *pBeginP + 1;
    *pBeginR = *pBeginR + 1;
}

int findNbrCSC(int u, int v, int *CSCindex, int *CSCedges)
{
    int index = -1;

    int first = CSCindex[u], last = CSCindex[u+1] - 1;
    int middle = (first+last)/2;

    while (first <= last) 
    {
        if (CSCedges[middle] < v)
            first = middle + 1;    
        else if (CSCedges[middle] == v) 
        {
            index = middle;
            break;
        }
        else
            last = middle - 1;
 
        middle = (first + last)/2;
    }
    
    return index;
}


#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#include<unistd.h> 
#include<libgen.h>
    
#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"



int main(int argc,char** argv)
{

    if (argc != 9)
    {
        printf("Incorrect number of arguments.\n");
        printf("./degeneracy_cliques -i <file_path> -t <type> -k <max_clique_size> -d <data_flag>\n");
        printf("file_path: path to file\n");
        printf("type: A/V/E. A for just k-clique information, V for per-vertex k-cliques, E for per-edge k-cliques\n");
        printf("max_clique_size: max_clique_size. If 0, calculate for all k.\n");
        printf("data_flag: 1 if information is to be output to a file, 0 otherwise.\n");
        return 0;
    }

    int n; // number of vertices
    int m; // 2x number of edges

    // char *opt = NULL;
    int opt;
    char *fpath = (char *)Calloc(1000, sizeof(char));
    char t;
    int flag_d;
    int max_k = 0;

    while((opt = getopt(argc, argv, ":i:t:k:d:")) != -1)  
    {  
        switch(opt)  
        {  
            case 'i':  
                // printf("In case i. optarg = %s\n", optarg);
                // strcpy(fpath, optarg);
                fpath = optarg;
                // printf("fpath = %s\n", fpath);
                break;
            case 't':  
                t = *optarg;
                if ((t != 'A') && (t != 'V') && (t != 'E'))
                {
                    printf("Incorrect type. Type should be A, V or E.\n");
                    return 0;
                }
                break;
            case 'k':  
                max_k = atoi(optarg);
                break;
            case 'd': 
                flag_d = atoi(optarg);
                if ((flag_d < 0) || (flag_d > 2))
                {
                    printf("Incorrect flag for data. Shoudld be 0, 1 or 2.\n");
                    return 0;
                }
                break;
            default:
                printf("In default case.\n");
                abort ();
        }
    }
    

    printf("Parsed all arguments. t = %c, max_k = %d, flag_d = %d. About to get graph.\n", t, max_k, flag_d);
    LinkedList** adjacencyList = readInGraphAdjListToDoubleEdges(&n, &m, fpath);
    printf("readInGraphAdjListToDoubleEdges done.\n"); 

    int i;

    char *gname = basename(fpath);

    char *lastdot = strrchr (gname, '.');
    if (lastdot != NULL)
        *lastdot = '\0';


    populate_nCr();
    printf("about to run runAndPrintStatsCliques \n"); 
    runAndPrintStatsCliques(adjacencyList, n, gname, t, max_k, flag_d);

    i = 0;
    while(i<n)
    {
        destroyLinkedList(adjacencyList[i]);
        i++;
    }

    Free(adjacencyList); 

    return 0;
}


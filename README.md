# sourcecode

This folder contains the source code for the submited paper "Scaling Up 𝑘-Clique Densest Subgraph Detection" to SIGMOD with paper ID 240. All algorithms are implemented in C++. 

To compile:
    make clean && make

To run:
    ./bin/degeneracy_cliques -i $data -t A -d 0 -k $K
    
$K is the parameter k, the number of vertices in a k-clique. 

$data is the path to the graph file. The first line of the graph file has to integers: the number of vertices (n) and the number of edges (m). This line is followed by m lines, each with two integers indicating the endpoints of an edge. 

For example, a triangle is represented by:
3 3
0 1
0 2
1 2

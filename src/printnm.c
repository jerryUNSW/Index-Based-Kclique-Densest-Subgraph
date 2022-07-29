#include<assert.h>
#include<stdio.h>
#include<stdlib.h>

int main()
{
    int n,m;

    if(scanf("%d", &n)!=1)
        exit(1);
    fprintf(stderr, "Number of vertices: %d\n", n);
    if(scanf("%d", &m)!=1)
        exit(1);
    fprintf(stderr, "Number of edges: %d\n", m/2);

    return 0;
}

#ifndef _DJS_MEMORYMANAGER_H_
#define _DJS_MEMORYMANAGER_H_
#include<stdlib.h>






#ifdef MEMORY_DEBUG
#include<stdio.h>


static void* MallocWithCheck(size_t x)
{
    #ifdef ALLOW_ALLOC_ZERO_BYTES
    void* retvalue = malloc(x);
    #else
    void* retvalue = malloc(max(x,1));
    #endif

    if(retvalue==NULL)
    {
        fprintf(stderr, "ERROR, malloc returned null pointer, that means we probably ran out of memory...\n");
        exit(1);
    }

    return retvalue;
};

/*! \brief Call calloc, and ensure that it returns non-NULL.

    \param x, the number of array slots to allocate

    \param x, the number of bytes to allocate in each slot

    \return a pointer to the allocated memory. If it is null, exit the program with error.

*/

static void* CallocWithCheck(size_t x, size_t y)
{
    #ifdef ALLOW_ALLOC_ZERO_BYTES
    void* retvalue = calloc(x,y); 
    #else
    void* retvalue = calloc(max(x,1),max(y,1)); 
    #endif

    if(retvalue==NULL)
    {
        fprintf(stderr, "ERROR, calloc returned null pointer, that means we probably ran out of memory...\n");
        exit(1);
    }

    return retvalue;
};

#define Malloc(x) MallocWithCheck(x)
#define Calloc(x,y) CallocWithCheck(x,y)
#define Free(x) free(x)

#else

    #ifdef ALLOW_ALLOC_ZERO_BYTES
    #define Malloc(x) malloc(x)
    #define Calloc(x,y) calloc(x,y)
    #define Free(x) free(x)

    #else

    #define Malloc(x) malloc(max(x,1))
    #define Calloc(x,y) calloc(max(x,1),max(y,1))
    #define Free(x) free(x)

    #endif // ALLOW_ALLOC_ZERO_BYTES
#endif // MEMORY_DEBUG

#endif


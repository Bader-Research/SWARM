#include "inputs_select.h"
#include "swarm_random.h"

void fill_same(int M, int *arr)
{
    register int
        i;

    for (i=0 ; i<M ; i++)
        arr[i] = i;

} 

void fill_linear(int M, int *arr, THREADED)
{
    register int
        i;

    for (i=0 ; i<M ; i++)
        arr[i] = (MYTHREAD*M) + i;
}

void fill_random(int M, int *arr, THREADED)
{
    register int
        i;

    for (i=0 ; i<M ; i++)
        arr[i] = SWARM_random(TH);
}
    


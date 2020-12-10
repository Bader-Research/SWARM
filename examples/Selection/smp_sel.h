#ifndef _SMP_SEL
#define  _SMP_SEL

#include "swarm.h"
#include "swarm_random.h"

#define TEST_MACH_LIMIT 0x10000   /*size of array */
#define INP_I_I       1
#define INP_I_S       2
#define INP_I_R       3
#define INP_I_N       4
#define INP_I_W       5

void partition_with_two_fr_i(int *A, int *B, int M, int s0, int s1,
                          int *a0, int *a1);  /*defined in select_fr.c*/
                          
void smp_test(int i, int *A,int inp,int *global_A,THREADED);
int smp_select_median_i(int M, int *global_A,THREADED);
int smp_select_i(int M, int total_n, int i,int *global_A,
	THREADED);

int locate_k(int c0,int *B_size,int B_max,int *k);
int gather(int *temp_array,int *l,int *r, int *global_A);
#endif

#ifndef _RADIXSORT_H
#define _RADIXSORT_H

#include "swarm.h"

void countsort_swarm(long q,
		     int *lKey,
		     int *lSorted,
		     int R,
		     int bitOff, int m,
		     THREADED);

#define radixsort_swarm(a,b,c,d)   radixsort_swarm_s3(a,b,c,d)
void radixsort_swarm_s3(long q,
			int *lKeys,
			int *lSorted,
			THREADED);
void radixsort_swarm_s2(long q,
			int *lKeys,
			int *lSorted,
			THREADED);
void radixsort20_swarm_s1(long q,
			  int *lKeys,
			  int *lSorted,
			  THREADED);
void radixsort20_swarm_s2(long q,
			  int *lKeys,
			  int *lSorted,
			  THREADED);

void radixsort_check(long q,
		     int *lSorted);

#endif


#include <stdio.h>

#include "radixsort.h"

#define bits(x,k,j) ((x>>k) & ~(~0<<j))

/****************************************************/
void countsort_swarm(long q,
		     int *lKey,
		     int *lSorted,
		     int R,
		     int bitOff, int m,
		     THREADED)
/****************************************************/
/* R (range)      must be a multiple of SMPS */
/* q (elems/proc) must be a multiple of SMPS */
{
    register int
	j,
	k,
        last, temp,
	offset;
    
    int *myHisto,
        *mhp,
        *mps,
        *psHisto,
        *allHisto;

    long x;

    myHisto  = (int *)SWARM_malloc(THREADS*R*sizeof(int), TH);
    psHisto  = (int *)SWARM_malloc(THREADS*R*sizeof(int), TH);

    mhp = myHisto + MYTHREAD*R;

    for (k=0 ; k<R ; k++)
      mhp[k] = 0;
    
    pardo(x, 0, q, 1)
      mhp[bits(lKey[x],bitOff,m)]++;

    SWARM_Barrier();

    pardo(k, 0, R, 1) {
      last = psHisto[k] = myHisto[k];
      for (j=1 ; j<THREADS ; j++) {
	temp = psHisto[j*R + k] = last + myHisto[j*R +  k];
	last = temp;
      }
    }

    allHisto = psHisto+(THREADS-1)*R;

    SWARM_Barrier();
    
    offset = 0;

    mps = psHisto + (MYTHREAD*R);
    for (k=0 ; k<R ; k++) {
      mhp[k]  = (mps[k] - mhp[k]) + offset;
      offset += allHisto[k];
    }

    SWARM_Barrier();
    
    pardo(x, 0, q, 1) {
      j = bits(lKey[x],bitOff,m);
      lSorted[mhp[j]] = lKey[x];
      mhp[j]++;
    }

    SWARM_Barrier();

    SWARM_free(psHisto, TH);
    SWARM_free(myHisto, TH);
}

/****************************************************/
void radixsort_check(long q,
		     int *lSorted)
/****************************************************/
{
  long i;

  for (i=1; i<q ; i++) 
    if (lSorted[i-1] > lSorted[i]) {
      fprintf(stderr,
	      "ERROR: q:%ld lSorted[%6ld] > lSorted[%6ld] (%6d,%6d)\n",
	      q,i-1,i,lSorted[i-1],lSorted[i]);
    }
}

/****************************************************/
void radixsort_swarm_s3(long q,
			int *lKeys,
			int *lSorted,
			THREADED)
/****************************************************/
{
  int *lTemp;

  lTemp = (int *)SWARM_malloc_l(q*sizeof(int), TH);
			
  countsort_swarm(q, lKeys,   lSorted, (1<<11),  0, 11, TH);
  countsort_swarm(q, lSorted, lTemp,   (1<<11), 11, 11, TH);
  countsort_swarm(q, lTemp,   lSorted, (1<<10), 22, 10, TH);

  SWARM_free(lTemp, TH);
}

/****************************************************/
void radixsort_swarm_s2(long q,
			int *lKeys,
			int *lSorted,
			THREADED)
/****************************************************/
{
  int *lTemp;

  lTemp = (int *)SWARM_malloc_l(q*sizeof(int), TH);
			
  countsort_swarm(q, lKeys,   lTemp,   (1<<16),  0, 16, TH);
  countsort_swarm(q, lTemp,   lSorted, (1<<16), 16, 16, TH);

  SWARM_free(lTemp, TH);
}

/****************************************************/
void radixsort20_swarm_s1(long q,
			  int *lKeys,
			  int *lSorted,
			  THREADED)
/****************************************************/
{
  countsort_swarm(q, lKeys,   lSorted, (1<<20),  0, 20, TH);
}

/****************************************************/
void radixsort20_swarm_s2(long q,
			     int *lKeys,
			     int *lSorted,
			     THREADED)
/****************************************************/
{
  int *lTemp;

  lTemp = (int *)SWARM_malloc_l(q*sizeof(int), TH);
			
  countsort_swarm(q, lKeys,   lTemp,   (1<<10),  0, 10, TH);
  countsort_swarm(q, lTemp,   lSorted, (1<<10), 10, 10, TH);

  SWARM_free(lTemp, TH);
}











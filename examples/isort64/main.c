#include "swarm.h"
#include "radixsort.h"
#include "create_input.h"

#define DEBUG              0
#define TIMING             1

#define TEST_WIDE          1

#if TEST_WIDE
#define ARR_SIZE_SORT   ((long)pow(2.0,33))
#else
#define ARR_SIZE_SORT   (1<<14)
#endif
#define MIN_TIME       0.000001



static void test_radixsort_swarm(long N1, THREADED) {
  int *inArr, *outArr;

#if TIMING
  double secs, tsec;
#endif

  inArr  = (int *)SWARM_malloc_l(N1 * sizeof(int), TH);
  outArr = (int *)SWARM_malloc_l(N1 * sizeof(int), TH);

  create_input_nas_swarm(N1, inArr, TH);

#if TIMING
  SWARM_Barrier();
  secs = get_seconds();
#endif

  radixsort_swarm(N1, inArr, outArr, TH);

#if TIMING
  secs = get_seconds() - secs;
  secs = max(secs,MIN_TIME);
  tsec = SWARM_Reduce_d(secs,MAX, TH);
  on_one {
    fprintf(stdout,"T: %3d n: %13ld  SSort: %9.6lf  (MB:%5ld)\n",
	    THREADS,N1,tsec,
	    ((long)ceil(((double)N1*2*sizeof(int))/(double)(1<<20))));
    fflush(stdout);
  }
#endif

  SWARM_Barrier();

  on_one radixsort_check(N1,  outArr);

  SWARM_Barrier();

  SWARM_free(outArr, TH);
  SWARM_free(inArr, TH);
}



static void *swarmtest(THREADED)
{
  long  i;

#if TEST_WIDE
#define TEST_INC (i<<1)
#else
#define TEST_INC (i+=8)
#endif
  
#if DEBUG
  fprintf(stdout,"PE%3d: SMP_main()\n",MYTHREAD);
  fflush(stdout);
#endif

  SWARM_Barrier();
  
  for (i = ((long)1<<12) ; i<=ARR_SIZE_SORT ; i = TEST_INC)
    test_radixsort_swarm(i, TH);

  SWARM_Barrier();

  /*******************************/
  /* End of program              */
  /*******************************/
  
  SWARM_done(TH);
}

int main(int argc, char **argv) 
{
  SWARM_Init(&argc,&argv);
  SWARM_Run((void *)swarmtest);
  SWARM_Finalize();
  return 0;
}

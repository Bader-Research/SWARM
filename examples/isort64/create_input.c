#include "create_input.h"
#include "nas_r.h"


void create_input_random_swarm(int myN, int *x, THREADED) {
  create_seq_random_swarm( 317*(MYTHREAD+17),
			    _NAS_MULT,   
			    myN,
			    x,
			    TH);   
}


void create_input_nas_swarm(int n, int *x, THREADED) {
  register int tsize, mynum, thtot;

  tsize = n / THREADS;
  mynum = MYTHREAD;
  thtot = THREADS;
  
  create_seq_swarm( find_my_seed( mynum,
				   thtot,
				   (n >> 2),
				   _NAS_SEED,    /* Random number gen seed */
				   _NAS_MULT),   /* Random number gen mult */
		     _NAS_MULT,                  /* Random number gen mult */
		     tsize,
		     x+(tsize*MYTHREAD),
		     TH);   

}




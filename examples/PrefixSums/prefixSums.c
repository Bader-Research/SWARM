#include "swarm.h"

#define DEBUG 0

static void input_gen (int **result,int k,int *n,THREADED) {
  
  int i;
  
  *n = (1<<k);

#if DEBUG
  printf("T%3d: input is %d\n",MYTHREAD,*n);
  fflush(stdout);
  SWARM_Barrier();
#endif
  
  *result = (int *) SWARM_malloc(*n*sizeof(int),TH);
  
  if (*result == NULL) {
    exit(-1);
  }
  on_one_thread {
    for(i=0 ; i<*n ; i++) {
      (*result)[i] = i+1;
    }
  }

#if DEBUG
  for (i=0 ; i<*n ; i++)
    fprintf(stdout,"T%3d: result[%5d]: %12d\n",MYTHREAD,i,(*result)[i]);
  fflush(stdout);
#endif
  
  SWARM_Barrier();
  return;
}


static void prefixSums_Tree(int *result, int k, int n, THREADED) {
  int h;
  int p;
  int j;
  int i;
  int x;
  int y;

#if DEBUG
  on_one_thread {
    fprintf(stdout,"k: %d\n", k); 
  }
#endif

  SWARM_Barrier();

  for (h=1; h<=k ; h++) { 
      
    p = (1<<(k-h));
    j = (1<<h); 

    pardo (x,1,p+1,1) {
      result[(j*x) -1] += result[((1<<(h-1))*((x<<1)-1)) -1];
    } 
       
    SWARM_Barrier();
  }       

  for (h=k-2 ; h>=0 ; h--) {

    p = (1<<(k-h));
    j=(1<<h); 

    pardo (y,3,p+1,2) {
      result[(j*y) -1] += result[(j*(y-1)) -1];
    }
      
    SWARM_Barrier();
  }
  
#if DEBUG
 on_one_thread {
    for (i=0 ; i<n; i++) {
      fprintf(stdout,"prefix_sum: %12d\n", result[i]);
    }
 }
#endif

 return;
} 

static void prefixSums_Block(int *result,int n,THREADED) {
  int j;
  int i;
  int *p;
  int r;
  int start;
  int end;
  int add_value;

  r = n / THREADS;

#if DEBUG
  fprintf(stdout,"T%3d: value of r: %5d",MYTHREAD, r);
#endif

  p = (int *)SWARM_malloc(NOSHARE(THREADS)*sizeof(int),TH);

#if DEBUG   
  if (p == NULL) {
    exit(-1);
  }

  if (r * (THREADS) != n) {
    fprintf(stderr,"error\n");
  }
#endif

  start =  MYTHREAD*r + 1;
  end   = (MYTHREAD+1)*r;
  
  for (j=start ; j<end ; j++) 
    result[j] += result[j-1];
  
  p[NOSHARE(MYTHREAD)] = result[end-1];

  SWARM_Barrier();
  
  on_one_thread {
    for (j=1 ; j<THREADS ; j++)
      p[NOSHARE(j)] += p[NOSHARE(j-1)];
  }
    
  SWARM_Barrier();

  if (MYTHREAD>0) {
    add_value=p[NOSHARE(MYTHREAD-1)];
    
    for (j=start-1 ; j<end ; j++)
      result[j] += add_value;
  }
  
  SWARM_Barrier();

#if DEBUG
  on_one {
    for (i=0 ; i<n ; i++)
      fprintf(stdout,"prefix_sums are: %12d\n",result[i]);
  }
#endif

}


static void *PrefixSums_main(THREADED) {
  int n;
  int k;   

  int *result;
  int *resultseq;
  int i;
  double t;

#if DEBUG
  fprintf(stdout,"THARGC = %d\n",THARGC);
  for (i=0;i<THARGC;i++)
    fprintf(stdout,"THARGV[%d]=%s\n",i,THARGV[i]);
  SWARM_Barrier();
#endif
 
  if (THARGC != 1) {
    fprintf(stderr,"ERROR: you must call the program with one argument\n");
    exit(-1);
  }
  
  k = atoi(THARGV[0]);
  
  input_gen(&result, k,&n,TH);

  SWARM_Barrier();
  
  t = get_seconds();
  prefixSums_Tree(result, k, n, TH);
  t = get_seconds() - t;

  t = SWARM_Reduce_d(t, MAX, TH);
  on_one_thread {
    fprintf(stdout,"pTree  P: %3d  k: %12d n: %12d Time: %9.6f\n",THREADS, k, n, t);
    fflush(stdout);
  }
  
  input_gen(&resultseq, k,&n,TH);

  SWARM_Barrier();
  
  on_one_thread {
    
    fprintf(stdout,"n: %d\n",n);
    fflush(stdout);
    for (i=1 ; i<n ; i++)
      resultseq[i] += resultseq[i-1];
    fprintf(stdout,"Checking the result:\n");
    fflush(stdout);
    for (i=0 ; i<n ; i++)
      if (result[i] != resultseq[i]) {
        fprintf(stdout,"ERROR: i = %12d  result[i]: %6d resultseq[i]: %6d\n",
                i, result[i], resultseq[i]);
      }
    fprintf(stdout,"Done Checking.\n");
    fflush(stdout);
  }
  
  SWARM_Barrier();

  input_gen(&result, k,&n,TH);

  SWARM_Barrier();

  t = get_seconds();
  prefixSums_Block(result, n, TH);
  t = get_seconds() - t;

  t = SWARM_Reduce_d(t, MAX, TH);
  on_one_thread {
    fprintf(stdout,"pBlock P: %3d  k: %12d n: %12d Time: %9.6f\n",THREADS, k, n, t);
    fflush(stdout);
  }
  
  input_gen(&resultseq, k,&n,TH);

  SWARM_Barrier();
  
  on_one_thread {
    
    fprintf(stdout,"n: %d\n",n);
    fflush(stdout);
    for (i=1 ; i<n ; i++)
      resultseq[i] += resultseq[i-1];
    fprintf(stdout,"Checking the result:\n");
    fflush(stdout);
    for (i=0 ; i<n ; i++)
      if (result[i] != resultseq[i]) {
        fprintf(stdout,"ERROR: i = %12d  result[i]: %6d resultseq[i]: %6d\n",
                i, result[i], resultseq[i]);
      }
    fprintf(stdout,"Done Checking.\n");
    fflush(stdout);
  }
  
  SWARM_Barrier();

  SWARM_free(result, TH);
  SWARM_free(resultseq, TH);

  SWARM_Barrier();
  
  SWARM_done(TH);
}

int main(int argc, char **argv) {
  SWARM_Init(&argc,&argv);
  SWARM_Run((void *)PrefixSums_main);
  SWARM_Finalize();
  return 0;
}

#include "smp_sel.h"

#include "inputs_select.h"

#include "select_randsamp.h"

#include "select_seq.h"

#define TRUE  1
#define FALSE 0

#define DEBUG1

static void *Selection_main(THREADED)
{
  int *A,*global_A;
  int A_size=TEST_MACH_LIMIT;

  SWARM_random_init(TH);
  SWARM_srandom(MYTHREAD+1, TH);

  /*Initialize SMP parameter such as PROCS*/
  SWARM_Barrier();
 
  /*malloc*/ 
  global_A = (int *)SWARM_malloc(THREADS*A_size*sizeof(int),TH);
  assert_malloc(global_A);
  /*PROCS=NODES*THREADS;*/

  A_size=A_size/THREADS;
  A=global_A+MYTHREAD*A_size;

  SWARM_random_init(TH);
  SWARM_srandom(MYTHREAD+1, TH);
       
  smp_test(A_size,A,INP_I_I,global_A,TH);
/*  
    smp_test(A_size,A,INP_I_S,global_A,TH);
    smp_test(A_size,A,INP_I_R,global_A,TH);
    smp_test(A_size,A,INP_I_N,global_A,TH);
*/

  SWARM_Barrier();
  SWARM_done(TH);
}

void smp_test(int i, int *A,int inp,int *global_A,THREADED) {
  double start_time,end_time;
  int result,j;

  fprintf(stdout,"MYTHREAD:%d\n",MYTHREAD);
  switch(inp) {
  case INP_I_I:  fill_same(i, A);    break;
  case INP_I_S:  fill_linear(i, A, TH);  break;
  case INP_I_R:  fill_random(i, A, TH);  break;
  default: fprintf(stderr,"Input class %d not known.\n",inp); exit(-1);
  }

#ifdef DEBUG2
  fprintf(stdout,"input array:\n");
  for( j=0;j<i;j++) fprintf(stdout,"%d,",A[j]);
  fprintf(stdout,"\n");
#endif

  /*barrier and get start time*/
   SWARM_Barrier();

  result=smp_select_median_i(i,global_A,TH);

  /*barrier and get end time, print the running time*/
   SWARM_Barrier();

  fprintf(stdout,"Result:%d\n",result);

   end_time=get_seconds();
  
  fprintf(stdout,"Running time for sorting %d words: %f seconds\n",
	i,end_time-start_time);

 
  return;
}


int smp_select_median_i(int M, int *global_A,THREADED) {
 int t; /*global size*/
 int median;

 t = THREADS*M; 
 median=(t>>1)+1;

 return smp_select_i(M, t, median,global_A,TH);
}


int smp_select_i(int M, int total_n, int i,int *global_A,THREADED)
{
  int *B,*C;
  int c_less,c_mid;
  int    *B_size,*my_Bsize;
  int  C_size,C_max,B_max,B_index;
  int  *l,*r;  /*local array:a[l..r]*/
  int *my_l,*my_r;
  int  n,ni; /*n:global_size,ni:local size*/
  int k[2]; /*the two pivots*/
  int c0,c1; /*index of the two pivots*/
  double Smag;
  int a[2];  /*in step 5, a[2]:local c_less,c_middle; */
  int flag; 
  int delta;
  int result;
  int *temp_array;
  int j;


#ifdef DEBUG1
  fprintf(stdout,"at the beginning of sel:TH:%d,M%d,i:%d,n:%d\n",
	  MYTHREAD,M,i,total_n);
#endif

  /*malloc for B[B_max]: sample array; C[C_max]:gatherd sample*/

  B_max = (int)ceil(M * pow((double)total_n, PICK_EPS)/total_n);
  C_max=THREADS*B_max;

  C=(int *)SWARM_malloc(sizeof(int)*C_max,TH);

  B_size=(int*)SWARM_malloc(sizeof(int)*THREADS,TH);
  my_Bsize=B_size+MYTHREAD;
    
  l=(int*)SWARM_malloc(sizeof(int)*THREADS,TH);
  my_l=l+MYTHREAD;

  r=(int*)SWARM_malloc(sizeof(int)*THREADS,TH);
  my_r=r+MYTHREAD;

  n = total_n;
  *my_l = M*MYTHREAD;
  *my_r = M-1+ (*my_l);
  
  temp_array=(int *)malloc(sizeof(int)*M);
  assert_malloc(temp_array);

  while (n > THREADS*THREADS) {
    /* Step 0 */
    ni = *my_r - *my_l + 1;

    /*step1: select sample */
    *my_Bsize=(int)ceil((double)ni * pow((double)n, PICK_EPS-1.0));
    B_index=MYTHREAD*(*my_Bsize);
    B=C+B_index;

    *my_Bsize = all_random_pick_exact_i(global_A+(*my_l), B, ni, n,TH);

#ifdef DEBUG2
    fprintf(stdout,"n:%d,my_Bsize:%d,B_index:%d\n",n,*my_Bsize,B_index);
#endif

    SWARM_Barrier();

    on_one_thread { 
      C_size=0;
      for(j=0;j<THREADS;j++)   C_size+=B_size[j];
    } 

    Smag = (double)C_size * (double)i / (double)n;
    delta=1;
    flag=TRUE;

    while(flag) { 
      /*step3: P0 work on the global C*/
      SWARM_Barrier();

      on_one_thread {
        c0   = (int)ceil(Smag - sqrt((double)C_size * delta));
        c1   = (int)ceil(Smag + sqrt((double)C_size * delta));

#ifdef DEBUG2
	fprintf(stdout,"ni:%d,smag:%f,c_size:%d,c0:%d,c1:%d\n",ni,Smag,C_size,c0,c1);
	fflush(stdout);
#endif

      	if (c0 < 0)       { flag=FALSE; c0=0;}   
      	if (c1 >= C_size)  { flag=FALSE;c1=C_size-1;}
	
        k[0]=select_mom_i(C,C_size,c0);
        k[1]=select_mom_i(C,C_size,c1);
      }

      /*step4: p0 bcastk1, k2*/
      SWARM_Barrier();
      flag=SWARM_Bcast_i(flag,TH);

      k[0]=SWARM_Bcast_i(k[0],TH);
      k[1]=SWARM_Bcast_i(k[1],TH);

#ifdef DEBUG2
      fprintf(stdout,"k:%d,%d\n",k[0],k[1]);
#endif
     
      /*step 5: partitopn array into three parts*/
      partition_with_two_fr_i(global_A+(*my_l), temp_array, ni, k[0],
			      k[1],&a[0],&a[1]);

#ifdef DEBUG2
      fprintf(stdout,"my_l:%d,a:%d,%d\n",*my_l,a[0],a[1]);
#endif
      SWARM_Barrier();

      /*step 6:Reduce, get c_less,c_mid*/
      c_less=SWARM_Reduce_i(a[0],SUM,TH);
      c_mid=SWARM_Reduce_i(a[1],SUM,TH);

#ifdef DEBUG2
      fprintf(stdout,"n:%d,i:%d,c_less:%d,c_mid:%d\n",n,i,c_less,c_mid);
      fflush(stdout);
#endif


      /*step7:check wether i is in (c_less,c_less+c_mid)*/
      if ((i > c_less) && (i<=c_less+c_mid)) {
        n  = c_mid;   
        *my_l= (*my_l)+a[0];
        (*my_r)  = (*my_l) + a[1] - 1;
        i -= c_less; 
	flag=FALSE;
      }
      else { if (flag) 
	delta*=PICK_MULT;
      else   
	if (i<=c_less) {
	  n=c_less;
	  (*my_r) =(*my_l)+a[0]-1;
	}
	else {
	  n=n-(c_less+c_mid);
	  (*my_l)  = (*my_l)+a[0]+ a[1];
	  i=i-(c_less+c_mid);
	}
      }

#ifdef DEBUG2
      fprintf(stdout,"n:%d,i:%d,my_l:%d,my_r:%d\n",n,i,*my_l,*my_r);
      fflush(stdout);
#endif 

    }
  }

  SWARM_Barrier();
#ifdef DEBUG2
  fprintf(stdout,"\nafter loop:MYTHREAD:%d,i:%d,my_l:%d,my_r:%d\n",
	  MYTHREAD,i,*my_l,*my_r);
  /*for(j=*my_l;j<=*my_r;j++) fprintf(stdout,"%d,",*(global_A+j));
    fprintf(stdout,"\n");*/
  fflush(stdout);
#endif

  free(temp_array);
  temp_array=(int *)malloc(sizeof(int)*n);
  assert_malloc(temp_array);

  /*step 8:Gather A[l,r]*/
  on_one_thread {
    n=gather(temp_array,l,r,global_A); 
  }
  SWARM_Barrier();

  /*step 9: P0 perform a sequential selection and broadcast the result*/
  on_one_thread { 
    result=select_mom_i(temp_array,n,i); 
  }

  SWARM_Barrier();
  result=SWARM_Bcast_i(result,TH);

  SWARM_free(C,TH);
  SWARM_free(B_size,TH);
  SWARM_free(l,TH);
  SWARM_free(r,TH);
  free(temp_array);

  return (result);
}

int locate_k(int c0,int *B_size,int B_max,int *k) {
  int temp_i,temp_index,j;

  temp_i=0;
  for(j=0;j<THREADS;j++) {
    temp_i+=B_size[j];
    if (c0 <temp_i) { 
      temp_index=j;   
      temp_i-=B_size[j];
      break;
    }
  }
  *k=B_max*temp_index+c0-temp_i;
      
  return 0;
}

int gather(int *temp_array,int *l,int *r, int *global_A) 
{
  int j,num,size=0;
  int *ptr=temp_array;

  for (j=0;j<THREADS;j++) {
    num=r[j]-l[j]+1;
    memcpy((char *) ptr,(char *)(global_A+l[j]),
	   num*sizeof(int));
    
    ptr=ptr+num;
    size+=num;
  }

  return size;
}


int main(int argc, char **argv) 
{
  SWARM_Init(&argc,&argv);
  SWARM_Run((void *)Selection_main);
  SWARM_Finalize();
  return 0;
}

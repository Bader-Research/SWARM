
#include "swarm_random.h"
#include "graph.h"
#include "stack.h"

/* try to do a balanced breadth-first search*/
void spanning_tree_breadth_B1(V* graph, int nVertices,THREADED)
{
#define S_POINTS (THREADS*THREADS*THREADS*2)
#define MYCOLOR (MYTHREAD+1)
#define GREY    (THREADS+1)
#define MAX_WALKS 10
#define MIN_WALKS 4
#define INDEX(a,b) ((a-1)*(THREADS)+(b-1))

  double start,end;
  double interval;

  int * color, first_time, work_to_steal,myroot,b,max, host, *stack,top,count,bottom=-1,i,j,n,root,walks,r,*starting_points;
  double power;
  int ** stack_M, **top_M, **bottom_M;
  int * finished;
  int * count_M;

  SWARM_random_init(TH);
  SWARM_srandom(MYCOLOR, TH);

  stack_M = SWARM_malloc(THREADS*sizeof(int *),TH);
  top_M = SWARM_malloc(THREADS*sizeof(int *),TH);
  bottom_M=SWARM_malloc(THREADS*sizeof(int *),TH);

  stack_M[MYTHREAD]=malloc(nVertices*sizeof(int));
  stack=stack_M[MYTHREAD];
  top_M[MYTHREAD]=&top;
  bottom_M[MYTHREAD]=&bottom;

  finished = SWARM_malloc(THREADS*sizeof(int),TH);
  finished[MYTHREAD]=0;
  count_M=SWARM_malloc(THREADS*sizeof(int),TH);

  color=SWARM_malloc(sizeof(int)*nVertices,TH); 

  pardo(i,0,nVertices,1){
    color[i]=0;
  }

  bottom=-1;
  top=-1;
  count=0;

  /*lets select a point to start in the graph*/
  on_one_thread {
    root=(SWARM_random(TH)%nVertices);
    power=1/(float)THREADS;
    walks=(int)pow(nVertices,power);
    if(walks<MIN_WALKS) walks=MIN_WALKS;
    if(walks>MAX_WALKS) walks=MAX_WALKS;
    color[root]=MYCOLOR;
    graph[root].parent=root;
  }

  root=SWARM_Bcast_i(root,TH);
  walks=SWARM_Bcast_i(walks,TH);
  starting_points=SWARM_malloc(sizeof(int)*S_POINTS,TH);
  
  /* printf("walks is %d, root is %d\n",walks,root);*/

  /*lets first generate some candidates for the threads to start their random walks*/

  myroot=root;
  on_one_thread{
    j=0;
    start=get_seconds();
    push(myroot,stack_M[j],top_M[j]);
    for(i=0;i<S_POINTS;i++)
      {
	starting_points[i]=myroot;
	r=SWARM_random(TH);
	if(r%2==0){
	  for(r=0;r<graph[myroot].n_neighbors;r++)
	    {
	      n=graph[myroot].my_neighbors[r];
	      if(color[n]==0){
		graph[n].parent=myroot;
		color[n]=MYCOLOR;
		myroot=n;
		count++;
		break;
	      }
	    }
	}
	else {
	  for(r=graph[myroot].n_neighbors-1;r>=0;r--)
	    {
	      n=graph[myroot].my_neighbors[r];
	      if(color[n]==0){
		graph[n].parent=myroot;
		color[n]=MYCOLOR;
		myroot=n;
		count++;
		break;
	      }      
	    }
	}
	push(myroot,stack_M[j],top_M[j]);
	j=(j+1)%THREADS;
	/*printf("seed %d\n",myroot);*/
      }
    end=get_seconds();
    interval=end-start;
    printf("Time used for generating starting points is %f\n",interval/1000000000);
  }
	       
  SWARM_Barrier();

  /*now each thread randomly picks a start point from the staring points*/
  r=SWARM_random(TH)%(THREADS*THREADS);
  myroot=starting_points[2*MYTHREAD*THREADS*THREADS+r];
    
  SWARM_Barrier();
  
  /* printf("my root is %d\n", myroot);*/
  push(myroot, stack, &top);
  count++;

  SWARM_Barrier();

  first_time=1;
  work_to_steal=0;

  while(first_time || work_to_steal)
    {
      /*if(work_to_steal) printf("stealing work\n");*/
      while(!is_empty( stack, &top, &bottom))
	{
	  n=pop(stack,&top,bottom);

	  if(n==-1){
	    printf("stack overflow\n");
	    top=-1;
	    bottom=-1;
	    break;
	  }

	  for(i=0;i<graph[n].n_neighbors;i++)
	    {
	      if(color[graph[n].my_neighbors[i]]==0) {/*found new frontier*/
		color[graph[n].my_neighbors[i]]=MYCOLOR;
		graph[graph[n].my_neighbors[i]].parent=n;
		push(graph[n].my_neighbors[i],stack,&top);
		count++;
	      }
	  
	    }
	}
      printf("Thread %d:done with this stack,current count is %d \n", MYTHREAD, count);

      if(first_time){
	first_time=0;
	finished[MYTHREAD]=1;
      }

      work_to_steal=0;
      max=-1;
      for(i=0;i<THREADS;i++)
	{
	  if(!finished[i]){
	    n=*(top_M[i]);
	    b=*(bottom_M[i]);
	    if(n<=0) continue;
	    if(n-b>max) {max=n-b;host=i;}
	  }
	}

      if(count>nVertices/THREADS) /* I already did my share*/
	  r=b+(n-b)/THREADS; /*I won't steal half of the work, that will put too much work on me. I take a 1/THREADS piece*/
      else  /* I didn't have chance to do my share, but lets make up.*/ 
	r=b+max((n-b)/THREADS,min(nVertices/THREADS-count,(n-b)/2));

      printf("Thread%d: thread %d's stack is %d tall\n",MYTHREAD,host,n-b);
      if(r<(*top_M[host])) {
	work_to_steal=1;		 
	(*bottom_M[host])=(r-1);
	while((r--)>b)
	  push(stack_M[host][r],stack,&top);
	break;
      }
    }

  count_M[MYTHREAD]=count;
  SWARM_Barrier();
  printf("Thread %d count is %d\n",MYTHREAD, count);  

  on_one_thread{
    int max=0, min=nVertices;
    for(i=0;i<THREADS;i++)
      {
	if(count_M[i]>max) max=count_M[i];
	if(count_M[i]<min) min=count_M[i];
      }
    printf("===span_breadth_B1:The difference between counts is %d\n",max-min);
  }
  SWARM_Barrier();

  SWARM_free(count_M,TH);
  SWARM_free(finished,TH);
  SWARM_free(color, TH);
  SWARM_free(starting_points, TH);
  SWARM_free(stack_M,TH);
  SWARM_free(top_M,TH);
  free(stack);
}





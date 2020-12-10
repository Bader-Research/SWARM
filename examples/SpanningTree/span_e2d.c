
#include "swarm_random.h"
#include "graph.h"
#include "stack.h"

/* Before we try to do a balanced breadth-first search, eliminate the connected d2 vertices*/
void spanning_tree_e2d(V* graph,int nVertices,THREADED)
{
#define S_POINTS (THREADS*THREADS*THREADS*2)
#define MYCOLOR (MYTHREAD+1)

  double start,end;
  double interval;

  int * color, first_time, work_to_steal,myroot,b,start_sr,start_sl, *stack,top,count,bottom=-1;
  int i,j,n,root,r,visited,neighbor,ret=0,n_neighbors;
  int ** stack_M, **top_M, **bottom_M;
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
  count_M=SWARM_malloc(THREADS*sizeof(int),TH);

  color=SWARM_malloc(sizeof(int)*nVertices,TH); 

  start=get_seconds();
  pardo(i,0,nVertices,1){
    color[i]=0;
  }
  end=get_seconds();
  interval=end-start;
  printf("The time used for setting up is %f \n",interval/1000000000);
 
  bottom=-1;
  top=-1;
  count=0;

  SWARM_Barrier();

  start=get_seconds();
  ret=eliminate_2d_vertices(graph,nVertices,&root,TH);
  SWARM_Barrier();
  end=get_seconds();
  interval=end-start;
  on_one_thread printf("METRICS:The time used for eliminating vertices is %f \n",interval/1000000000);
  SWARM_Barrier();
  
  if(ret==0) root=(SWARM_random(TH)%nVertices);
  /*lets select a point to start in the graph*/
  on_one_thread {
    color[root]=MYCOLOR;
    graph[root].parent=root;
    myroot=root;
    j=0;
    push(myroot,stack_M[j],top_M[j]);
    j=(j+1)%THREADS;
    i=0;
    for(i=0;i<S_POINTS;i++)
      {
	n_neighbors=graph[myroot].n_neighbors;
	r=SWARM_random(TH);
	if(r%2==0){
	  for(r=0;r<n_neighbors;r++)
	    {
	      visited++;
	      n=graph[myroot].my_neighbors[r];
	      if(color[n]==0){
		graph[n].parent=myroot;
		color[n]=MYCOLOR;
		myroot=n;
		count++;
		push(myroot,stack_M[j],top_M[j]);
		/*printf("sub tree : %d \n",myroot);*/
		j=(j+1)%THREADS;
		break;
	      } 
	    }
	  if(r==n_neighbors){
	    r=(SWARM_random(TH)%n_neighbors);
	    myroot=graph[myroot].my_neighbors[r];
	  }
	}
	else {
	  for(r=n_neighbors-1;r>=0;r--)
	    {
	      visited++;
	      n=graph[myroot].my_neighbors[r];
	      if(color[n]==0){
		graph[n].parent=myroot;
		color[n]=MYCOLOR;
		myroot=n;
		count++;
		push(myroot,stack_M[j],top_M[j]);
		/*printf("sub tree : %d \n",myroot);*/
		j=(j+1)%THREADS;
		break;
	      }     
	    }
	  if(r<0){
	    r=(SWARM_random(TH)%n_neighbors);
	    myroot=graph[myroot].my_neighbors[r];
	  }
	}		
      }
    end=get_seconds();
    interval=end-start;
  }
	       
  SWARM_Barrier();
  start=get_seconds();
 
  first_time=1;
  work_to_steal=0;
  start_sr=MYTHREAD; /*when I am out of work, where do i start to search. r towards right, l towards left*/
  start_sl=MYTHREAD;

  while(first_time || work_to_steal)
    {
      /*if(work_to_steal) printf("stealing work\n");*/
      while(!is_empty( stack, &top, &bottom))
	{
	  n=pop(stack,&top,bottom);
	  if(n==-1){
	    printf("THREAD %d:stack overflow\n",MYTHREAD);
	    bottom=-1;
	    top=-1;
	    break;
	  }
	  visited+=graph[n].n_neighbors;
	  for(i=0;i<graph[n].n_neighbors;i++)
	    {
	      neighbor=graph[n].my_neighbors[i];
	      if(color[neighbor]==0) {/*found new frontier*/
		color[neighbor]=1;
		graph[neighbor].parent=n;
		push(neighbor,stack,&top);
		count++;
	      }
	  
	    } 
	}
      /*printf("Thread %d:done with this stack, current count is %d \n", MYTHREAD, count);*/

      if(first_time) first_time=0;

      work_to_steal=0;
      if(MYTHREAD%2==0)
	{
	  for(j=0;j<THREADS;j++)
	    {
	      i=(j+start_sr)%THREADS; /*start searching in circular from my neighbor*/
	      if(i==MYTHREAD) continue;

	      n=*(top_M[i]);
	      b=*(bottom_M[i]);	
	      if(n-b<THRESHOLD) continue;	      	       
	      if(count>nVertices/THREADS) /*I did my share*/
		r=b+(n-b)/THREADS;
	      else r=b+max((n-b)/THREADS,min(nVertices/THREADS-count,(n-b)*3/4));
	      /*printf("Thread%d: thread %d's stack is %d tall\n",MYTHREAD,i,n-b);*/
	      if(r<(*top_M[i])) {
		(*bottom_M[i])=max(r-1,-1);
		/*printf("THREAD %d: I am taking %d elements\n",MYTHREAD, r-b);*/
		work_to_steal=1;		 
		while((r--)>b)
		  push(stack_M[i][r],stack,&top);
		break;		
	      }
	    }
	} else{
	  for(j=THREADS-1;j>0;j--)
	    {
	      i=(j+start_sl)%THREADS;
	      if(i==MYTHREAD) continue;
	      n=*(top_M[i]);
	      b=(*bottom_M[i]);
	      if(n-b<THRESHOLD) continue;	     
	      /*printf("Thread%d: thread %d's stack is %d tall\n",MYTHREAD,i,n-b);*/
	      if(count>nVertices/THREADS) /*I did my share*/
		r=b+(n-b)/THREADS;
	      else r=b+max((n-b)/THREADS,min(nVertices/THREADS-count,(n-b)*3/4));
	      if(r<(*top_M[i])) {
		(*bottom_M[i])=max(r-1,-1);
		/*printf("THREAD %d: I am taking %d elements\n",MYTHREAD, r-b);*/
		work_to_steal=1;
		while((r--)>b)
		  push(stack_M[i][r],stack,&top);
		break;
	      } 
	    }
	}
      start_sr=(start_sr+1)%THREADS;
      start_sl=(start_sl-1+THREADS)%THREADS;

    }
  count_M[MYTHREAD]=count;
  end=get_seconds();
  interval=end-start;
  printf("THREAD %d: I AM DONE..., time used: %f\n",MYTHREAD,interval/1000000000);
  SWARM_Barrier();
  printf("Thread %d count is %d, visited is %d\n",MYTHREAD, count,visited);  

  on_one_thread{
    int max=0, min=nVertices;
    for(i=0;i<THREADS;i++)
      {
	if(count_M[i]>max) max=count_M[i];
	if(count_M[i]<min) min=count_M[i];
      }
    printf("METRICS===span_e2d:The difference between counts is %d\n",max-min);
  }
  SWARM_Barrier();

  SWARM_free(count_M,TH);
  SWARM_free(color, TH);
  SWARM_free(stack_M,TH);
  SWARM_free(top_M,TH);
  free(stack);

}





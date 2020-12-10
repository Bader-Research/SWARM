
#include "swarm_random.h"
#include "graph.h"
#include "stack.h"

/*see if we can do a breadth-search by each thread, currently El and n_edges are not used*/
void spanning_tree_breadth(V* graph, int nVertices,THREADED)
{
#define S_POINTS (THREADS*THREADS*THREADS)
#define MYCOLOR (MYTHREAD+1)

  int * color, myroot, *stack,top,count,i,j,n,n_neighbors,root,r,bottom;
  int * count_M;
  int ** stack_M, **top_M, **bottom_M;

  SWARM_random_init(TH);
  SWARM_srandom(MYCOLOR, TH);

  count_M=SWARM_malloc(THREADS*sizeof(int),TH);
  color=SWARM_malloc(sizeof(int)*nVertices,TH);
  stack_M = SWARM_malloc(THREADS*sizeof(int *),TH);
  top_M = SWARM_malloc(THREADS*sizeof(int *),TH);
  bottom_M=SWARM_malloc(THREADS*sizeof(int *),TH);

  stack_M[MYTHREAD]=malloc(nVertices*sizeof(int));
  stack=stack_M[MYTHREAD];
  top_M[MYTHREAD]=&top;
  bottom_M[MYTHREAD]=&bottom;

  pardo(i,0,nVertices,1){
    color[i]=0;
  }

  SWARM_Barrier();

  top=-1;
  bottom=-1;
  count=0;

  /*lets select a point to start in the graph*/
  on_one_thread {
    root=(SWARM_random(TH)%nVertices);
    color[root]=MYCOLOR;
  }

  root=SWARM_Bcast_i(root,TH);
  printf("root is %d\n",root);

  /*lets first generate some candidates for the threads to start their random walks*/
	       
  SWARM_Barrier();
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
	      } ;
	    }
	  if(r==n_neighbors) {
	    r=(SWARM_random(TH)%n_neighbors);
	    myroot=graph[myroot].my_neighbors[r];
	  }
	}
	else {
	  for(r=n_neighbors-1;r>=0;r--)
	    {
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
    printf("generating sub tree done\n");
  }

  SWARM_Barrier();
  while(!is_empty( stack, &top,&bottom))
    {
      n=pop(stack,&top,bottom);
      if(n==-1) {
	printf("stack overflow\n");
	bottom=-1;
	top=-1;
	break;
	}
      if(MYTHREAD%2==0)
	{
	  for(i=0;i<graph[n].n_neighbors;i++)
	    {
	      if(color[graph[n].my_neighbors[i]]==0) {/*found new frontier*/
		color[graph[n].my_neighbors[i]]=MYCOLOR;
		push(graph[n].my_neighbors[i],stack,&top);
		count++;
	      }
	  
	    }
	} else{
	  for(i=graph[n].n_neighbors-1;i>=0;i--)
	    {
	      if(color[graph[n].my_neighbors[i]]==0) {/*found new frontier*/
		color[graph[n].my_neighbors[i]]=MYCOLOR;
		push(graph[n].my_neighbors[i],stack,&top);
		count++;	      
	      }
	    }
      
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
    printf("METRICS:=====span_breadth:The difference between counts is %d\n",max-min);
  }
  SWARM_Barrier();
  SWARM_free(count_M,TH);
  SWARM_free(stack_M,TH);
  SWARM_free(top_M,TH);
  SWARM_free(bottom_M, TH);
  SWARM_free(color, TH);
  free(stack);
}






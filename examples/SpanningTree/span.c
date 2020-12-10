#include <sys/types.h>
#include "graph.h"
#include "lock.h"

static V* G;
static V* T;
static E* El;
static int n_edge,max_d_node;

#define NANO 1000000000
#define DEBUG_BICONN 0 
#define DEBUG_VERIFY 0

static void *Span_main(THREADED)
{
  int i,t,n_vertices,max=0;
  double start,end;
  double interval,total=0;
  char * input_file;
  int * D;

  /*initialize graph from input file */  

  input_file = THARGV[0];

  start = get_seconds();
  on_one_thread{
	t = initialize_graph(input_file,&G,&n_vertices);
	if(t!=0) exit(0);
  }
  n_vertices=SWARM_Bcast_i(n_vertices,TH);
  end = get_seconds();
  interval=end-start;
   printf("Time for initialization is %f\n",interval/NANO);
  SWARM_Barrier();

  on_one_thread{
   max=0;
   for(i=0; i<n_vertices; i++)
   {
     if(G[i].n_neighbors>max){
	 max=G[i].n_neighbors; 
	 max_d_node=i;
     }
   }
  }
  SWARM_Barrier();
 
#if 0
  start = get_seconds();
  initialize_graph_edgelist(G, n_vertices,&El,&n_edge, TH);
  end = get_seconds();
  interval=end-start;
  printf("Time for initialization(graph edge list) is %f\n",interval/NANO);
  SWARM_Barrier();

  SWARM_Barrier();
  start = get_seconds();
  spanning_tree_CRCW(G,El,n_vertices,n_edge,TH);
  end = get_seconds();
  interval=end-start;
  total+=interval;
  printf("Time for spanning_tree_CRCW is %f\n",interval/NANO);
  SWARM_Barrier();

  pardo(i,0,n_edge,1)
  {
     El[i].workspace=0;
  }

  SWARM_Barrier();
  start = get_seconds();
  spanning_tree_random01(G,El,n_vertices,n_edge,TH);
  end = get_seconds();
  interval=end-start;
  total+=interval;
  printf("Time for spanning_tree_random01 is %f\n",interval/NANO);
  SWARM_Barrier();
#endif

  on_one_thread printf("METRICS:%d THREADS, file %s\n", THREADS, input_file);
  SWARM_Barrier();
  start = get_seconds();
  spanning_tree_CREW(G,n_vertices,&T,TH);
  end = get_seconds();
  interval=end-start;
  total+=interval;
  SWARM_Barrier();
  on_one_thread printf("METRICS:Time for spanning_tree_CREW is %f\n",interval/NANO);

  SWARM_Barrier();
  start = get_seconds();
  spanning_tree_breadth(G,n_vertices,TH);
  end = get_seconds();
  interval=end-start;
  total+=interval;
  SWARM_Barrier(); 
  on_one_thread printf("METRICS:Time for spanning_tree_breadth is %f\n",interval/NANO);

  SWARM_Barrier();
  start = get_seconds();
  spanning_tree_breadth_B2(G,n_vertices,max_d_node,TH);
  end = get_seconds();
  interval=end-start;
  total+=interval;
  SWARM_Barrier();
  on_one_thread printf("METRICS:Time for spanning_tree_breadth_B2 is %f\n",interval/NANO);

  SWARM_Barrier();
  start = get_seconds();
  spanning_tree_e2d(G,n_vertices,TH);
  end = get_seconds();
  interval=end-start;
  total+=interval;
  SWARM_Barrier();
  on_one_thread printf("METRICS:Time for spanning_tree_e2d is %f\n",interval/NANO);

#if DEBUG_VERIFY
  D=SWARM_malloc(sizeof(int)*n_vertices,TH);
  pardo(i,0,n_vertices,1) D[i]=G[i].parent;
  SWARM_Barrier();
  pardo(i,0,n_vertices,1)
    {
      while(D[i]!=D[D[i]]){
	D[i]=D[D[i]];
      }
    } 
  SWARM_Barrier();
  printf("pointer jumping done\n");


  pardo(i,0,n_vertices,1)
    {
      if(D[i]!=D[0]) {
	printf("error\n");
	break;
      }
    }
  printf("done verifying\n");
  SWARM_free(D,TH);
#endif

  on_one_thread{
	delete_graph(G,n_vertices);
	if(El) free(El);
  } 

  SWARM_Barrier();
  SWARM_done(TH);
}

int main(int argc, char **argv) 
{
  SWARM_Init(&argc,&argv);
  SWARM_Run((void *)Span_main);
  SWARM_Finalize();
	return 0;
}

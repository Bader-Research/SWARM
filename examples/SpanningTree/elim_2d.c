#include "graph.h"
#include <stdlib.h>

#define DEBUG_ELIM 1

static int in_my_range(int i, int n, THREADED)
{
  if(i>=(n/THREADS)*MYTHREAD && i<(n/THREADS)*MYTHREAD+(n/THREADS) ) return 1;
  else return 0;
}

/*eliminate degree 2 vertices. The basic strategy is to connect the two neighbors of the degree 2 vertex, thus roll out the degree 2 vertex. The n_neighbors of the degreee 2 vertex is set to 0. 
  Return value: 0 if there are no vertices eliminated, 1 there are vertices eliminated.
  In the case that there are vertices eliminated, root will be set to point to one of the vertex that is still in the graph ( not eliminated).
*/
int eliminate_2d_vertices(V* G, int n_vertices,int * root,THREADED)
{
  int i,j,n,m,p,w,my_starting_point, old_degree,index1,index2, end1, end2,count=0,total_count=0;
  int * starting_point;
  int  *count_matrix;
  count_matrix=SWARM_malloc(sizeof(int)*THREADS,TH);

  /* I still need to test if it is worth while to eliminate the d-2 vertices */
  pardo(i,0,n_vertices,1)
    {
      if(G[i].n_neighbors==2 && G[G[i].my_neighbors[0]].n_neighbors<=2 && G[G[i].my_neighbors[1]].n_neighbors<=2) count++;
    }
  count_matrix[MYTHREAD]=count;
  SWARM_Barrier();
  on_one_thread{
    for(i=1;i<THREADS;i++) count_matrix[0]+=count_matrix[i];
  }
  SWARM_Barrier();
  if(count_matrix[0]<n_vertices/THREADS){
    printf("Not enough d2 vertices to eliminate. ignore...\n");
    SWARM_free(count_matrix,TH);
    return (0);
  }
  
  starting_point=SWARM_malloc(sizeof(int),TH);
  
  n=(n_vertices/THREADS)*MYTHREAD;
  old_degree=G[n].n_neighbors;
  if(G[n].n_neighbors==2) G[n].n_neighbors=1; 
  /*This is specially for the whole graph is just a string. This way we will break the string into THREADS sub strings, hopefully equally sized, then remove all the intermediate d2 nodes. If we didn't do the break, then all the threads have to find the two ends of the whole string, which will give no speed up when we have more processors*/
  my_starting_point=n;

  SWARM_Barrier();

  /* for each non-degree-2 vertices who has neighbors that are degree-2, we want to follow the chain to the other end (where the degree is not 2), and link the two ends together hence remove all the in-between degree 2 vertices.*/


  printf("eliminate_2d_vertices:before pardo\n");

  pardo(i,0,n_vertices,1)
    {
      if(G[i].n_neighbors!=2) {	
	if(i==my_starting_point) m=old_degree;
	else m=G[i].n_neighbors;
	for(j=0;j<m;j++)
	  {
	    n=G[i].my_neighbors[j];
	    if(G[n].n_neighbors==2){
	      count=0;
	      end1=i; p=i; index1=j;
	      while(G[n].n_neighbors==2){
		count++;
		w=n;
		if(G[n].my_neighbors[0]==p) n=G[n].my_neighbors[1];
		else n=G[n].my_neighbors[0];
		p=w;
	      }
	      end2=n;
	      if(in_my_range(end2,n_vertices,TH) || end2>end1 ) {
		/*I will do the linking two ends if the two ends are in my charge, or the other end is beyond my charge. This is to avoid two processors try to link the two ends at the same time*/
		//printf("THREAD %d: In my range (%d,%d)\n", MYTHREAD,end1,end2);
		for(index2=0;G[end2].my_neighbors[index2]!=p;index2++);
		G[end1].my_neighbors[index1]=end2;
		G[end2].my_neighbors[index2]=end1;
		total_count+=count;
	      } else {};//printf("THREAD %d: Not in my range (%d,%d)\n", MYTHREAD,end1,end2);
	    }
	  }
      }
    }

  printf("removed vertices:%d\n",total_count);
#if DEBUG_ELIM
  count_matrix[MYTHREAD]=total_count;
  SWARM_Barrier();
  on_one_thread{
    for(i=1;i<THREADS;i++)
      count_matrix[0]+=count_matrix[i];
    printf("METRICS:the total removed vertices are:%d\n", count_matrix[0]);
  }
  SWARM_free(count_matrix,TH);
  printf("Leaving eliminate_2d_vertices\n"); 
#endif

  SWARM_Barrier();
  G[my_starting_point].n_neighbors=old_degree;
  (*starting_point)=my_starting_point;
  SWARM_Barrier();
  *(root)=(*starting_point);
  SWARM_free(starting_point,TH);
  return 1;
}







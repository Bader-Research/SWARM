#include "listrank.h"
#include <sys/types.h>

#define NANO 1000000000
static int nlist;

static void ListFB(list_t * list)
{ int i;

  for (i = 0; i < nlist; i++) list[i].succ = i + 1;
  list[nlist-1].succ = 0;
}


static void ListBF(list_t *list)
{ int i;

  list[0].succ = 0;
  list[1].succ = 0;
  for (i = 2; i <nlist; i++) list[i].succ = i - 1;
}


static void ListRandom(list_t *list)
{ 
  int i,j,t, *buf;
  double s;
  
  buf = malloc(sizeof(int)*nlist);
  
  for(i=0; i<nlist; i++)
  	buf[i]=i;


  for(i=0; i<nlist;i++)
  	{
		s = drand48();
        t = (int)(i+(s*nlist-i));
        if(t!=i)
        {
			j = buf[i];
			buf[i]=buf[t];
			buf[t]=j;
        }

	}
 
  for(i=0;i<nlist-1;i++)
  	list[buf[i]].succ=buf[i+1];

  list[buf[nlist-1]].succ=0;
  
  free(buf);
}

static void *Listrank_main(THREADED)
{
	list_t * list;
	int k = 8;

	if (THARGC != 1) {
	  on_one 
	    fprintf(stderr,"Usage: Listrank -t <number of threads> -- <number of elements in list>\n");
	  exit(-1);
	}
	on_one {
	  nlist = atoi(THARGV[0]);
	  //k = 8;
	}
	SWARM_Barrier();
	
	list = SWARM_malloc(sizeof(list_t)*nlist,TH);
	
	on_one{
		ListFB(list);
	}
	SWARM_Barrier();
	
	on_one fprintf(stdout,"Ranking ListFB (n=%d)... ",nlist);
	printf("k = %d\n", k);
	list_ranking(nlist, k, list, TH);
	SWARM_Barrier();
	on_one fprintf(stdout,"Done.\n");

    	on_one{
		ListBF(list);
	}
	SWARM_Barrier();

	on_one fprintf(stdout,"Ranking ListBF (n=%d)... ",nlist);
	list_ranking(nlist, k, list, TH);
	on_one fprintf(stdout,"Done.\n");
	SWARM_Barrier();
	
	on_one{
		ListRandom(list);
	}
	SWARM_Barrier();
	
	on_one fprintf(stdout,"Ranking ListRandom (n=%d)... ",nlist);
	list_ranking(nlist, k, list, TH);
	on_one fprintf(stdout,"Done.\n");
	SWARM_Barrier();

	SWARM_free(list,TH);
	SWARM_done(TH);
}

int main(int argc, char **argv) 
{
  SWARM_Init(&argc,&argv);
  SWARM_Run((void *)Listrank_main);
  SWARM_Finalize();
  return 0;
}

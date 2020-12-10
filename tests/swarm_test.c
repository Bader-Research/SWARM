
#include <swarm.h>
#include <swarm_random.h>

static void routine (THREADED)
{
     int i;
     
	#ifdef HAVE_SPRNG
     	on_one
       	printf("Using SPRNG:\n");
	#else
    	on_one
      	printf("Not using SPRNG:\n");
	#endif

     SWARM_Barrier();
     SWARM_random_init(TH);
     SWARM_srandom(MYTHREAD+1,TH);
     for (i=0; i<10; i++)
       printf("T%3d: sample %2d: %12ld\n",MYTHREAD,i,SWARM_random(TH));
}

int main (int argc, char **argv)
{
     SWARM_Init(&argc,&argv);

     SWARM_Run ((void *)routine);

     SWARM_Finalize();

     return 0;
}


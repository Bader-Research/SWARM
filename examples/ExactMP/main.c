/******************************************************************************
* main.c                                                                      *
* Checks arguments and forwards requests to "framework.c"                     *
******************************************************************************/

#include "header.h"

void *ExactMP_main(THREADED) {

   int i, initBestCost, printer, skipSize, lhsOffSet;
   time_t st, et;

   if (MYTHREAD==0 && (THARGC==0 || THARGC > 3)) {
    printf("USAGE: ExactMP [-t <# threads>] -- [-q] <input matrix> [desired # of arrangements for randomization]\n");
    exit(0);
   }

   if (MYTHREAD==0 && (strcmp(THARGV[0], "-q")==0)) {
     verbose=0;
     printf("\nQuiet Mode turned ON\n");
     fflush(stdout);
   } else if (MYTHREAD==0 && (strcmp(THARGV[0], "-q")!=0)) {
     verbose=1;
     printf("\nQuiet Mode turned OFF\n");
     fflush(stdout);
   }
   if (MYTHREAD==0) {
      st = time((time_t *)NULL);
      printf("\n\t\tSTARTING JOB %s WITH %d THREADS @ TIME %s",               \
                    THARGV[!verbose], THREADS, ctime(&st));
      NUM_TBR_TREES = 600;
      if (THARGC==3) {
        NUM_TBR_TREES = atoi(THARGV[2]);
      } else if (THARGC==2 && verbose) {
        NUM_TBR_TREES = atoi(THARGV[1]);
      }
   }

   setFrameWork(THARGV[!verbose],TH);

   if (matrix.num_taxa < 8) {
      if (MYTHREAD==0) {
        initBestCost = bestCost;
      	branchAndBound(TH);
      }
   } else {
        //if (MYTHREAD==0) bestCost = 1000;
        initBestCost = bestCost;
      	bnb(TH);
   }
   SWARM_Barrier();

   printer = 0;
   if (MYTHREAD==printer) {
       et = time((time_t *)NULL);
       /* printf("\n\t\t-o-o-o-o-\n"); */
       printf("\n\nPNI Sites: %d, Const Sites: %d, PI Sites: %d\n",           \
                        numPNISites, numConstSites, matrix.num_pars_inf_sites);
       printf("InitScore: %d  started  @ %s",initBestCost+pNiCost, ctime(&st));
       printf("BestScore: %d converged @ %s",bestCost+pNiCost, ctime(&et));
       printf("Number of best solutions retained: %d\n",solQCtr+1);
       printf("\n\t\t-o-o-o-o-\n\n");
   }
   SWARM_Barrier();
   printf("Thread %d decomposed %.0f nodes\n",MYTHREAD, loadDistribution[MYTHREAD]);
   fflush(stdout);
   SWARM_Barrier();

   if (MYTHREAD==printer) {
       printf("\n\t\t-o-o-o-o-\n\n");
       lhsOffSet = 4;
       skipSize = 5 + matrix.num_pars_inf_sites;
       for (i=0; i <= solQCtr; i++) {
         displayTreeInNexusFormat(solQ[i],lhsOffSet, skipSize, TH);
         printf("\n");
       }
       printf("\n\t\t-o-o-o-o-\n\n");
   }
   SWARM_Barrier();

   if (verbose) {
     printf("Thread %d terminating\n",MYTHREAD);
     fflush(stdout);
   }

   SWARM_done(TH);
}

int main(int argc, char **argv) 
{
  SWARM_Init(&argc,&argv);
  SWARM_Run((void *)ExactMP_main);
  SWARM_Finalize();
  return 0;
}

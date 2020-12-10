/******************************************************************************
* framework.c                                                                 *
* Sets the framework for the algorithm in the following manner:               *
* 	(1) Compute initial task with 3 taxa                                  *
* 	(2) Preprocess the input data and assign State Encoding to the taxa   *
* 	(3) Parsimonious Site Reordering                                      *
* 	(4) Get initial best cost using SG algorithm                          *
* 	(5) Decide the taxa addition order using maxMini algorithm            *
* 	(6) Push the first task to the task stack                             *
******************************************************************************/
#include "header.h"

void setFrameWork(char *argv, THREADED) {

  int i, skipSize, lhsOffSet, rhsOffSet;

  if (MYTHREAD==0) {

     if (verbose) {
       printf("Setting framework...\n");
       fflush(stdout);
     }

     pthread_mutex_init(&(lock), NULL);
     solQCtr = -1;

     loadDistribution = (double *)calloc(THREADS, sizeof(double));
     assert(loadDistribution);

     request = (int *)calloc(THREADS, sizeof(int));
     assert(request);

     taskID = (int *)calloc(THREADS, sizeof(int));
     assert(taskID);

     for (i=0; i < THREADS; i++) {
       request[i] = 0;
     }
     for (i=0; i < THREADS; i++) {
       taskID[i] = i;
     }
     nextID = THREADS;

     numBestSolutions = 0;
     taskSTop = NULL;

     taxaQHead = NULL;
     taxaQTail = NULL;

     tAddOrdQHead = NULL;
     tAddOrdQTail = NULL;

     bestCostStackHead = NULL;

     taxaQueue = (int *)calloc(matrix.num_taxa, sizeof(int));
     assert(taxaQueue);

     if (verbose) {
       printf("\nPreprocessing input...\n");
     }
     preProcess(argv, TH);

     if (verbose) {
       printf("\nReordering sites...\n");
     }
     reorderSites(TH);

     skipSize = 5 + matrix.num_pars_inf_sites,
     lhsOffSet = 4;
     rhsOffSet = 1+matrix.num_taxa;
     solQ = (int ***)calloc(keepTrees, sizeof(int **));
     assert(solQ);
     for (i=0; i < keepTrees; i++) {
         solQ[i] = (int **)calloc(2, sizeof(int *));
         assert(solQ[i]);
         solQ[i][0] =                                                         \
         (int *)calloc((2*matrix.num_taxa*skipSize)+lhsOffSet,sizeof(int));
         assert(solQ[i][0]);

         solQ[i][1]=(int *)calloc((2*matrix.num_taxa*skipSize)+1,sizeof(int));
         assert(solQ[i][1]);
     }

     /*printf("\nApplying NJ Algorithm...\n");
     nj(TH); */

     if (verbose) {
     printf("\nApplying modified greedy algorithm...");
     }
     getIstBestCost(TH);

     if (verbose) {
     printf("\t\t[BEST SCORE %d]\n",bestCost+pNiCost);
     }

     if (verbose) {
     printf("\nDeciding taxa addition order...");
     }
     taxaAddOrder(TH);
  
     if (verbose) {
     printf("\nApplying eckDay greedy algorithm...");
     }
     eckDay(TH);
     if (verbose) {
     printf("\t\t[BEST SCORE %d]\n",bestCost+pNiCost);
     }

     if (verbose) {
     printf("\nApplying TBR on current best solution...");
     }
     tbr(bestCostStackHead->bCNode, TH);
     if (verbose) {
     printf("\t[BEST SCORE %d]\n",bestCost+pNiCost);
     }
     if (verbose) {
     printf("\nApplying 3-level %d-arrng/level randomized algorithm...\n",   \
                                                            NUM_TBR_TREES);
     }
     randOptimize(TH);

     printf("\nPNI Sites: %d, Const Sites: %d, PI Sites: %d\n",               \
                numPNISites, numConstSites, matrix.num_pars_inf_sites);
     printf("\nApplying B&B with a score of %d...\n\n",bestCost+pNiCost);
  }
  SWARM_Barrier();
}

void printTaxaStatus (struct taskNode *srcTaskNode, THREADED) {
  struct taxaQ *tNode;
  /*printf("\n\t\t\t\t\t---------------------------------------\n");
  printf("\t\t\t\t\t    TAXA STATUS FOR ABOVE TREE \n");
  printf("\t\t\t\t\t---------------------------------------\n");*/
  printf("\n\t\t\t\t    REMAINING TAXA FOR ABOVE TREE: ");
  for (tNode=srcTaskNode->taxaRemHead; tNode != NULL; tNode = tNode->next) {
    printf("%d ",tNode->taxa);
  }
  printf("\n\n\t\t\t\t    TAXA DONE FOR ABOVE TREE: ");
  for (tNode=srcTaskNode->taxaDoneHead; tNode != NULL; tNode = tNode->next) {
    printf("%d ",tNode->taxa);
  }
  printf("\n");
}

void printTree (struct subTask *sourceSubTask, int pos, THREADED) {
  int lc, rc, pr, i, arrIndx;
  struct subTask *sT;
  
  sT = sourceSubTask;

  printf("\n\t\t\t\t\t---------------------------------------\n");
  printf("\t\t\t\t\t     TREE IN PREORDERED ARRAY FORM\n");
  printf("\t\t\t\t\t---------------------------------------\n");
  printf("\n\tHighest Int Node: %d\t",sT->highestInternalNode);
  printf("Last taxon added: %d\t",sT->lastTaxaAdded);
  if (pos != -1) {
    printf("in position: %d\t",pos);
  }
  printf("Size of array : %d\t",sT->sizeOfSubTaskArray);
  printf("Cost of sub task: %d\n",sT->costOfSubTask);
  /*printf("\n\t\tARRAY ELEMENTS:\t");
  for (i=0; i < sT->sizeOfSubTaskArray; i++) {
	printf("[%d] ",sT->node[i].id);
  }
  printf("\n\n");*/
  printf("\t\t\tLCHILD\tNODE\tRCHILD\tPARENT\tCOST AT THIS NODE\t{STATE VECTOR}\n");
  for(arrIndx=sourceSubTask->sizeOfSubTaskArray-1;arrIndx>=0;arrIndx=arrIndx-1){
      rc = sT->node[arrIndx].rChildIndx;
      lc = sT->node[arrIndx].lChildIndx;
      pr = sT->node[arrIndx].parentIndx;
      if (sT->node[arrIndx].id > 0) {
        if (sT->node[arrIndx].parentIndx == ROOT) {
          printf("\t\t\t null\t%2d\t null\t root\t\t%2d\t\t", \
                      sT->node[arrIndx].id,sT->node[arrIndx].costAtThisNode);
	}
	else {
          printf("\t\t\t null\t%2d\t null\t%3d\t\t%2d\t\t",  \
      sT->node[arrIndx].id, sT->node[pr].id,sT->node[arrIndx].costAtThisNode);
	}
      }
      else {
        if (sT->node[arrIndx].parentIndx == ROOT) {
          printf("\t\t\t%3d\t%2d\t%3d\t root\t\t%2d\t\t",sT->node[lc].id,   \
       sT->node[arrIndx].id, sT->node[rc].id,sT->node[arrIndx].costAtThisNode);
	}
	else {
          printf("\t\t\t%3d\t%2d\t%3d\t%3d\t\t%2d\t\t",sT->node[lc].id,   \
       sT->node[arrIndx].id, sT->node[rc].id, sT->node[pr].id,  \
       sT->node[arrIndx].costAtThisNode);
	}
      }
      printf("{");
      for (i=0; i < matrix.num_pars_inf_sites; i++) {
         printf("%d ",sT->node[arrIndx].stateVector[i]);
      }
      printf("}\n");
  }
  printf("\t\t\t%3d\t%2d\t%3d\t%3d\t\t%2d\t\t",                              \
    sT->node[sT->root.lChildIndx].id, sT->root.id,                           \
    sT->node[sT->root.rChildIndx].id, sT->root.id, sT->root.costAtThisNode);
  printf("{");
  for (i=0; i < matrix.num_pars_inf_sites; i++) {
    printf("%d ",sT->root.stateVector[i]);
  }
  printf("}\n");
}

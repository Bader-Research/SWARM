#include "header.h"

void branchAndBound(THREADED) {
  struct taskNode *parentTaskNode, *childTaskNode; 
  struct subTask *parentTask, *childTask; 
  int currTaxa, pos, numElements, purdLowBnd;
  unsigned long added = 0;
  unsigned long checked = 0;
  int initBestCost;

  initBestCost = bestCost;
  flushBestCostStack(TH);
                                   /* While there are more incomplete tasks */
  while (taskSTop != NULL) {
                                   /* Get task, & extract tree node from it */
    parentTaskNode = popTask(); 
    parentTask = parentTaskNode->partialTask; 
    numElements = parentTask->sizeOfSubTaskArray;
#ifdef BnB_VERY_VERBOSE
    printf("\n\t\t<><><><><><><><><><><><><><><><><><><><>");
    printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
    printf("\t\t\t   POPPED TASK (BC %d) ",bestCost);
    printf("WITH NUM TAXA (%d) ",numElements/2+1);
    printf("CURRENT COST W/O PNI COST(%d)\n",parentTask->costOfSubTask);
    printf("\n\t\t<><><><><><><><><><><><><><><><><><><><>");
    printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
#ifdef BnB_VERY_VERY_VERBOSE
    printTree(parentTask, -1, TH);
    printTaxaStatus (parentTaskNode, TH);
#endif
#endif
                                   /* While more taxa to be added, continue */
    if ((currTaxa = dequeueTaxaQ(parentTaskNode)) != -1) {
                                   /* For each possible position */
        for (pos=0; pos < numElements; pos++) { 
           if (parentTask->node[pos].id!= -1) { 
                                   /* Generate child task by adding taxon 
                                      in given position */
              childTask=genNextTask(parentTask, currTaxa, pos, TH);
	      ++checked;
                                   /* If cost of child task < best cost */
              if (childTask->costOfSubTask < bestCost) {
		  ++added;
                                   /* Check if this is leaf task => 
                                      this is current best solution */
                  if (abs(childTask->highestInternalNode)+2==matrix.num_taxa) {
                                   /* Flush best cost stack, add this afresh */
                     flushBestCostStack(TH);
		     pushBestCostNode(childTask);
		     bestCost = childTask->costOfSubTask;
	             printf("BETTER BC FOUND: %d\n",bestCost);
#ifdef BnB_VERY_VERBOSE
    printf("\n\t\t\t\t\t==================================================\n");
    printf("\t\t\t\t\t  BETTER BC FOUND, FLUSH STACK AND ADD THIS (BC %d)\n", \
                                bestCost);
    printf("\t\t\t\t\t==================================================\n");
#ifdef BnB_VERY_VERY_VERBOSE
    printTree(childTask, pos, TH);
#endif
#endif
		  }
                                    /* If not leaf task, test lower bound, 
                                       if passes stack it */
		  else {
                     childTaskNode = copyTaxaQ(parentTaskNode); 
                     childTaskNode->partialTask = childTask; 
                                    /* Add the current taxa added into the 
                                       done taxa list */
		     childTaskNode->taxaDoneTail->next =                     \
                                   (struct taxaQ*)malloc(sizeof(struct taxaQ));
		     childTaskNode->taxaDoneTail =                           \
                                   childTaskNode->taxaDoneTail->next;
		     childTaskNode->taxaDoneTail->taxa = currTaxa;
		     childTaskNode->taxaDoneTail->next = NULL;

                                    /* Apply Purdoms lower bounding function */
	             purdLowBnd = purdomsLowerBnd(childTaskNode);	
                                    /* Best sol from this task > bestCost, 
                                       solution impossible */
		     if (childTaskNode->partialTask->costOfSubTask+purdLowBnd \
                                                       > bestCost) {
#ifdef BnB_VERY_VERBOSE
    printf("\n\t\t\t======================================================");
    printf("=======================================\n");
    printf("\t\t\t    UNPRODUCTIVE TASK, FAILED LOW BND TEST,PRUNED ");
    printf("(BC %d), (PLB %d)",bestCost,purdLowBnd);
    printf(", (LBC %d)",childTaskNode->partialTask->costOfSubTask+purdLowBnd);
    printf(", (NUM TAXA %d)\n",                                              \
                        childTaskNode->partialTask->sizeOfSubTaskArray/2+1);
    printf("\n\t\t\t======================================================");
    printf("=======================================\n");
#ifdef BnB_VERY_VERY_VERBOSE
    printTree(childTask, pos, TH);
    printTaxaStatus (childTaskNode, TH);
#endif
#endif
                         freeStackNode(childTaskNode, TH);
		     }
                                     /* Best sol from this task < bestCost, 
                                        solution may be possible */
		     else {
                         pushTask(childTaskNode); 
#ifdef BnB_VERY_VERBOSE
   printf("\n\t\t\t=====================================================");
   printf("=============================\n");
   printf("\t\t\t    PUSHED FOLLOWING INCOMPLETE TASK ");
   printf("(BC %d), (PLB %d), ",bestCost, purdLowBnd);
   printf("(LBC %d), ",childTaskNode->partialTask->costOfSubTask+purdLowBnd);
   printf("(NUM TAXA %d)\n",childTaskNode->partialTask->sizeOfSubTaskArray/2+1);
   printf("\n\t\t\t=====================================================");
   printf("=============================\n");
#ifdef BnB_VERY_VERY_VERBOSE
   printTree(childTask, pos, TH);
   printTaxaStatus (childTaskNode, TH);
#endif
#endif
		     }
		  }
              }
#ifdef GET_ONE_BEST_SOLUTION
	      else {
		  freeTree(childTask);
	      }
#endif
#ifdef GET_ALL_SOLUTIONS
                                     /* If cost is equal to best cost */
              else if (childTask->costOfSubTask == bestCost) {
		  ++added;
                                     /* Check if this is leaf task => this 
                                        is another best solution */
                  if (abs(childTask->highestInternalNode)+2==matrix.num_taxa) {
		     pushBestCostNode(childTask);
#ifdef BnB_VERY_VERBOSE
    printf("\n\t\t\t\t\t=========================================\n");
    printf("\t\t\t\t\t  ANOTHER BEST SOLUTION, STACK IT (BC %d)\n",bestCost);
    printf("\t\t\t\t\t=========================================\n");
#ifdef BnB_VERY_VERY_VERBOSE
    printTree(childTask, pos, TH);
#endif
#endif
		  }
                                     /* Else if this is not leaf task => 
                                        stack it */
		  else {
                     childTaskNode = copyTaxaQ(parentTaskNode); 
                     childTaskNode->partialTask = childTask; 
                                     /* Add the current taxa added into the 
                                        done taxa list */
		     childTaskNode->taxaDoneTail->next =                     \
                                  (struct taxaQ*)malloc(sizeof(struct taxaQ));
		     childTaskNode->taxaDoneTail =                           \
                                   childTaskNode->taxaDoneTail->next;
		     childTaskNode->taxaDoneTail->taxa = currTaxa;
		     childTaskNode->taxaDoneTail->next = NULL;

                                     /* Apply Purdoms lower bounding function */
	             purdLowBnd = purdomsLowerBnd(childTaskNode);	
                                     /* Best sol from this task > bestCost, 
                                        solution impossible */
		     if (childTaskNode->partialTask->costOfSubTask+purdLowBnd \
                                             > bestCost) {
#ifdef BnB_VERY_VERBOSE
   printf("\n\t\t\t================================================");
   printf("===========================================\n");
   printf("\t\t\t UNPRODUCTIVE TASK, FAILED LOW BND TEST,PRUNED ");
   printf("(BC %d), (PLB %d), ",bestCost,purdLowBnd);
   printf("(LBC %d), ",childTaskNode->partialTask->costOfSubTask+purdLowBnd);
   printf("(NUM TAXA %d)\n",childTaskNode->partialTask->sizeOfSubTaskArray/2+1);
   printf("\n\t\t\t================================================");
   printf("===========================================\n");
#ifdef BnB_VERY_VERY_VERBOSE
   printTree(childTask, pos, TH);
   printTaxaStatus (childTaskNode, TH);
#endif
#endif
                         freeStackNode(childTaskNode, TH);
		     }
                                       /* Best sol from this task < bestCost, 
                                          solution may be possible */
		     else {
                         pushTask(childTaskNode); 
#ifdef BnB_VERY_VERBOSE
    printf("\n\t\t\t=================================================");
    printf("=====================\n");
    printf("\t\t\t    PUSHED FOLLOWINING INCOMPLETE TASK ");
    printf("(BC %d), (PLB %d), (LBC %d)\n",bestCost, purdLowBnd, \
                 childTaskNode->partialTask->costOfSubTask+purdLowBnd);
    printf("\n\t\t\t=================================================");
    printf("=====================\n");
#ifdef BnB_VERY_VERY_VERBOSE
    printTree(childTask, pos, TH);
    printTaxaStatus (childTaskNode, TH);
#endif
#endif
		     }
		  }
              }
                                       /* Else if cost of child task is more 
                                          than best cost */
              else if (childTask->costOfSubTask > bestCost) {
                                       /* Prune it, unproductive task */
          
#ifdef BnB_VERY_VERBOSE
    printf("\n\t\t\t================================");
    printf("======================================\n");
    printf("\t\t\t   UNPRODUCTIVE TASK PRUNED, ORIG COST ");
    printf("(%d) > (BC %d), (NUM TAXA %d)\n",childTask->costOfSubTask, \
                     bestCost, childTask->sizeOfSubTaskArray/2+1);
    printf("\n\t\t\t================================");
    printf("======================================\n");
#ifdef BnB_VERY_VERY_VERBOSE
    printTree(childTask, pos, TH);
#endif
#endif
		  freeTree(childTask);
	      }
#endif
           }
        }
    }
                                        /* Free parent task, all children 
                                           generated and tested */
    freeStackNode(parentTaskNode, TH);
  }
                                        /* Print best solutions */
  printBestSolutions(TH);

}


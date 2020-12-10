/******************************************************************************
* randOpt.c                                                                   *
* Computes better upper bound through randomized algorithm                    *
******************************************************************************/

#include "header.h"
void randOptimize(THREADED) {
	int ctr=0, **arr, numTbrTrees, randNum, oldBc, prntCtr;

                                       /* generate arbitrary addition order 
	                                  sequences, number # of such sequences
                                          can be specified in header file */
	arr = (int **)calloc(NUM_TBR_TREES, sizeof(int *));
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		arr[numTbrTrees]=(int *)calloc(matrix.num_taxa-3, sizeof(int));
		for (ctr=0; ctr < matrix.num_taxa-3; ctr++) {
			arr[numTbrTrees][ctr] = -1;
		}
	}

	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		ctr=0;
		/*srandom(numTbrTrees+1);*/
		srandom(random()+1);
		while (ctr < matrix.num_taxa-3) {
			randNum	= random()%(matrix.num_taxa-3);
			if (arr[numTbrTrees][randNum] == -1) {
				arr[numTbrTrees][randNum] = ctr+4;
				++ctr;
			}
		}
	}
                                        /* apply EDG, maxmini, IEDG algorithms 
	                                   on each of those addition orders 
                                           for better solution */
	oldBc = bestCost+pNiCost;
        if (verbose) {
	printf("\tRandomization level I: ");
	fflush(stdout);
        }
        prntCtr=0;
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		nQ4RandOpt(arr, numTbrTrees, TH);
		TBRSmartGreedy(TH);
                freeRandTaxaQ(TH);
                if (verbose) {
                ++prntCtr;
                if (prntCtr==10) {
		   printf(".");
		   fflush(stdout);
                   prntCtr=0;
                }
                }
		if (oldBc > bestCost+pNiCost) {
                   if (verbose) {
		   printf("\n\t\tImproved [BEST SCORE %d] ",bestCost+pNiCost);
		   /* printf(" @ iteration %d of level I", numTbrTrees); */
		   fflush(stdout);
                   }
		   oldBc = bestCost+pNiCost;
		}
	}
        if (verbose) {
	printf("\n");
        }

	oldBc = bestCost+pNiCost;
        if (verbose) {
	printf("\tRandomization level II: ");
	fflush(stdout);
        }
        prntCtr=0;
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		TBReckDay(arr, numTbrTrees, TH);
                if (verbose) {
                ++prntCtr;
                if (prntCtr==10) {
		   printf(".");
		   fflush(stdout);
                   prntCtr=0;
                }
                }
		if (oldBc > bestCost+pNiCost) {
                   if (verbose) {
		   printf("\n\t\tImproved [BEST SCORE %d] ",bestCost+pNiCost);
		   /* printf(" @ iteration %d of level II", numTbrTrees); */
		   fflush(stdout);
                   }
		   oldBc = bestCost+pNiCost;
		}
	}
        if (verbose) {
	printf("\n");
        }

	oldBc = bestCost+pNiCost;
        if (verbose) {
	printf("\tRandomization level III: ");
	fflush(stdout);
        }
        prntCtr=0;
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		nQ4RandOpt(arr, numTbrTrees, TH);
		TBRmaxMin(TH);
                freeRandTaxaQ(TH);
                if (verbose) {
                ++prntCtr;
                if (prntCtr==10) {
		   printf(".");
		   fflush(stdout);
                   prntCtr=0;
                }
                }
		if (oldBc > bestCost+pNiCost) {
                   if (verbose) {
		   printf("\n\t\tImproved [BEST SCORE %d] ",bestCost+pNiCost);
		   /* printf(" @ iteration %d of level III", numTbrTrees); */
		   fflush(stdout);
                   }
		   oldBc = bestCost+pNiCost;
		}
	}
        if (verbose) {
	printf("\n");
	fflush(stdout);
        }

    for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
	free(arr[numTbrTrees]);
    }
    free(arr);
}

void freeRandTaxaQ(THREADED) {
  struct taxaQ* temp_node;
  while (taxaQHead!=NULL) {
    temp_node = taxaQHead->next;
    free(taxaQHead);
    taxaQHead = temp_node;
  }

  taxaQHead = NULL;
  taxaQTail = NULL;

}

void TBReckDay(int **arr, int row, THREADED) {
	int taxaOver = 3, ctr;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

	sourceTree = genIstSubTask(TH);
	ctr = 0;
	taxaOver += 1; 
	mMTPWL = getEckDayMinTree(sourceTree, arr[row][ctr], TH);
	ctr = ctr + 1;
	while (taxaOver < matrix.num_taxa) {
		childTree =                                                   \
         genNextEckDayTask(sourceTree, mMTPWL.taxa, mMTPWL.pos, TH); 
		freeTree(sourceTree);
		sourceTree = childTree;
		taxaOver += 1; 
		mMTPWL = getEckDayMinTree(sourceTree, arr[row][ctr], TH);
		ctr = ctr + 1;
	}
  	childTree =                                                           \
         genNextEckDayTask(sourceTree, mMTPWL.taxa, mMTPWL.pos, TH); 
	if (childTree->costOfSubTask < bestCost) {
		  bestCost = childTree->costOfSubTask;
	}

	correctEckDayIndices(childTree, TH);
	tbr(childTree, TH);

	freeTree(childTree);
	freeTree(sourceTree);
}

void correctEckDayIndices (struct subTask *sourceSubTask, THREADED) {
  int arrIndx, indx, indxOfMinusOne;
                                        /* for all nodes starting from RHS */
  for(arrIndx=sourceSubTask->sizeOfSubTaskArray-1;arrIndx>=0;arrIndx=arrIndx-1){
                                        /* if a leaf node, stack it */
    if (sourceSubTask->node[arrIndx].id > 0) {
      sourceSubTask->node[arrIndx].rChildIndx = LEAF;
      sourceSubTask->node[arrIndx].lChildIndx = LEAF;
      push(arrIndx);
    }
                                        /* if node is internal node, pop 2 
                                           nodes from stack, & push this 
                                           node on stack. 2 poped nodes are 
                                           left & right children resp. */
    else if (sourceSubTask->node[arrIndx].id < 0) {
      if (sourceSubTask->node[arrIndx].id == -1) {
        indxOfMinusOne = arrIndx;
      }

      indx = pop();
      sourceSubTask->node[arrIndx].lChildIndx = indx;
      sourceSubTask->node[indx].parentIndx = arrIndx;

      indx = pop();
      sourceSubTask->node[arrIndx].rChildIndx = indx;
      sourceSubTask->node[indx].parentIndx = arrIndx;

      push(arrIndx);
    }
  }
                                        /* roots of rhs & lhs subtrees */
  sourceSubTask->root.lChildIndx = pop();
  sourceSubTask->root.rChildIndx = pop();

  indx = sourceSubTask->node[sourceSubTask->root.rChildIndx].id;

  sourceSubTask->node[sourceSubTask->root.rChildIndx].id = -1;
  sourceSubTask->node[indxOfMinusOne].id = indx;
  sourceSubTask->rhtSubTreeRootIndx = sourceSubTask->root.rChildIndx;

  sourceSubTask->node[0].parentIndx = ROOT;
  sourceSubTask->node[sourceSubTask->rhtSubTreeRootIndx].parentIndx = ROOT;
}

void nQ4RandOpt(int **arr, int row, THREADED) {
  struct taxaQ* temp_node;
  int ctr;
  for (ctr=0; ctr < matrix.num_taxa-3; ctr++) {
    if (taxaQHead==NULL) {
        taxaQHead = (struct taxaQ*)malloc(sizeof(struct taxaQ));
        taxaQHead->next = NULL;
      	taxaQHead->taxa = arr[row][ctr];
        taxaQTail = taxaQHead;
    }
    else {
        temp_node = (struct taxaQ*)malloc(sizeof(struct taxaQ));
        taxaQTail->next = temp_node;
      	temp_node->taxa = arr[row][ctr];
        temp_node->next = NULL;
        taxaQTail = temp_node;
    }
  }
}

void TBRmaxMin(THREADED) {
	int taxaOver = 3;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

	sourceTree = genIstSubTask(TH);
	mMTPWL = maxMini(sourceTree, TH);
	taxaOver += 1;
	while (taxaOver < matrix.num_taxa) {
		childTree =                                                  \
            genNextTask(sourceTree, mMTPWL.taxa, mMTPWL.pos, TH);
		freeTree(sourceTree);
		sourceTree = childTree;
		mMTPWL = maxMini(sourceTree, TH);
		taxaOver += 1;
	}

	childTree=genNextTask(sourceTree,mMTPWL.taxa,mMTPWL.pos,TH);
	freeTree(sourceTree);
	tbr(childTree, TH);
	if (childTree->costOfSubTask < bestCost) {
	  bestCost = childTree->costOfSubTask;
	}
	freeTree(childTree);
}

void TBRSmartGreedy(THREADED) {
	int taxaOver = 3;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

	sourceTree = genIstSubTask(TH);
	mMTPWL = smartGreedy(sourceTree, TH); 
	taxaOver += 1; 

	while (taxaOver < matrix.num_taxa) { 
		childTree =                                                  \
               genNextTask(sourceTree, mMTPWL.taxa, mMTPWL.pos, TH); 
		freeTree(sourceTree);
		sourceTree = childTree;
		mMTPWL = smartGreedy(sourceTree, TH);
		taxaOver += 1; 
	}

	childTree=genNextTask(sourceTree,mMTPWL.taxa,mMTPWL.pos,TH);
	freeTree(sourceTree);
	tbr(childTree, TH);
	if (childTree->costOfSubTask < bestCost) {
	  bestCost = childTree->costOfSubTask;
	}
	freeTree(childTree);
}

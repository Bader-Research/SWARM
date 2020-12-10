/******************************************************************************
* eckDay.c                                                                    *
* Decides initial upper bound using Eck & Dayhoff greedy algorithm            *
******************************************************************************/
#include "header.h"

struct maxMiniTaxaPosWLength getEckDayMinTree(struct subTask *sourceTree,    \
                                              int taxaNum, THREADED) {
	struct subTask *childTree;
	struct maxMiniTaxaPosWLength mMTPWL; 
	int currentLength = 0;
	int prevLength = 0;
	int pos;
	int currMinPos; 
	int currMinTaxa; 
	int arraySize = sourceTree->sizeOfSubTaskArray; 
	 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
	printf("\n\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        printf("~~~~~~~~~~~~\n");
	printf("\t\t\t\t CHECKING MINIMUM TREE LENGTH FOR ALL POSITIONS FOR");
        printf("TAXA: %d\n",taxaNum);
	printf("\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        printf("~~~~~~~~~\n");
#endif
                                        /* generate 1st case */
	childTree = genNextEckDayTask(sourceTree, taxaNum, 1, TH); 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
	printTree(childTree, 0, TH);
#endif
	prevLength = childTree->costOfSubTask;
	freeTree(childTree);

	currMinPos = 1;
	currMinTaxa = taxaNum; 

	for (pos=2; pos < arraySize; pos++) { 
	 childTree = genNextEckDayTask(sourceTree, taxaNum, pos, TH); 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
	 printTree(childTree, pos, TH);
#endif
	 currentLength = childTree->costOfSubTask; 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
	 printf("\n\t\t->->->->->->RETAIN MINIMUM LENGTH FROM CURRENT");
         printf(" & PREVIOUS LENGTHS<-<-<-<-<-<-\n");
         printf("\n\t\t\tCurr {MinLen:%3d @ Pos:%3d",currentLength, pos);
         printf(" w/ Taxa:%3d} vs. Prev {MinLen:%3d",taxaNum, prevLength);
         printf("@ Pos:%3d w/ Taxa:%3d}\n\n",currMinPos, currMinTaxa);
#endif
                                        /* if length < prev length, keep it */
	  if (currentLength < prevLength) {
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
	   printf("\t\tCurr Lgt (%d) < ",currentLength);
	   printf("Prev Lgt (%d)=> ",prevLength);
#endif
	   prevLength = currentLength; 
	   currMinPos = pos; 
	   currMinTaxa = taxaNum; 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
	   printf("Update Curr Lgt to: %d ",currentLength);
	   printf(" keeping Taxa: %d, Pos: %d\n\n",currMinTaxa,currMinPos);
#endif
	  } 
	  freeTree(childTree);
	}
		 
                                        /* return position & tree */
	mMTPWL.taxa = currMinTaxa;
	mMTPWL.pos = currMinPos;
	mMTPWL.length = prevLength;
		 
	return (mMTPWL); 
}

void computeEckDayCost (struct subTask *sourceSubTask, THREADED) {

  int rChildIndx, lChildIndx, arrIndx;
  struct subTaskArrayNode *leftChild, *rightChild, *parent;

                                        /* Loop to root of subtree */
  for (arrIndx=sourceSubTask->sizeOfSubTaskArray-1; arrIndx >= 1; arrIndx--) {
                                        /* For all parent nodes, */
	if (sourceSubTask->node[arrIndx].id < 0) {
                                        /* Gather children information */
		lChildIndx = sourceSubTask->node[arrIndx].lChildIndx;
		leftChild = &sourceSubTask->node[lChildIndx];

		rChildIndx = sourceSubTask->node[arrIndx].rChildIndx;
		rightChild = &sourceSubTask->node[rChildIndx];

		parent = &sourceSubTask->node[arrIndx];
		performFitchOp(leftChild , parent, rightChild, TH);
	}
  }
                                        /* Store cost of left-right trees 
                                           in the roots of subtrees, 
                                           compute cost of root node now */
  sourceSubTask->root.rChildIndx = 1;
  sourceSubTask->root.lChildIndx = 0;
  performFitchOp(&sourceSubTask->node[0] ,                                   \
                 &sourceSubTask->root, &sourceSubTask->node[1], TH);
  sourceSubTask->costOfSubTask = sourceSubTask->root.costUntilThisNode;
}

void eckDay(THREADED) {
	int taxaOver = 3, ctr, stackPtr, 
            lChildIndx, rChildIndx, *locationStack, keepLoc;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

	sourceTree = genIstSubTask(TH);
                                        /* apply maxmini first */
	taxaOver += 1; 
	mMTPWL = getEckDayMinTree(sourceTree, taxaOver, TH);
	while (taxaOver < matrix.num_taxa) {
#ifdef ECK_AND_DAY_VERBOSE
		printf("\t\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~\n");
		printf("\t\t\t\t\t\tADDING NEXT TAXA ORDER\n");
		printf("\t\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~");
#endif
		childTree = genNextEckDayTask(sourceTree, mMTPWL.taxa, mMTPWL.pos, TH); 
		freeTree(sourceTree);
		sourceTree = childTree;
#ifdef ECK_AND_DAY_VERBOSE
		printTree(sourceTree,mMTPWL.pos, TH);
#endif
		taxaOver += 1; 
		mMTPWL = getEckDayMinTree(sourceTree, taxaOver, TH);
	}
  	childTree = genNextEckDayTask(sourceTree, mMTPWL.taxa, mMTPWL.pos, TH); 
	if (childTree->costOfSubTask < bestCost) {
                                        /* Following lines of code just make 
                                           sure that root of RHS sub tree has 
                                           an internal ID of -1 */
		stackPtr = -1;
		locationStack = (int *)calloc(childTree->sizeOfSubTaskArray,  \
                                        sizeof(int));
		for (ctr=childTree->sizeOfSubTaskArray-1; ctr >= 0; ctr--) {
			if (childTree->node[ctr].id > 0) {
				stackPtr = stackPtr + 1;
				locationStack[stackPtr] = ctr;
			}
                                        /* if id < 0, cp sv & pop twice,
                                           compute cost & stack neg node */
			else if (childTree->node[ctr].id < 0) {
				if (childTree->node[ctr].id == -1) {
					keepLoc = ctr;
				}
				lChildIndx = locationStack[stackPtr];
				stackPtr = stackPtr - 1;
	
				rChildIndx = locationStack[stackPtr];
				stackPtr = stackPtr - 1;
	
				stackPtr = stackPtr + 1;
				locationStack[stackPtr] = ctr;
			}
		  }
		  lChildIndx = locationStack[stackPtr];
		  stackPtr = stackPtr - 1;

		  rChildIndx = locationStack[stackPtr];
		  stackPtr = stackPtr - 1;
	
		  if (childTree->node[rChildIndx].id != -1) {
			ctr = childTree->node[rChildIndx].id;
			childTree->node[rChildIndx].id =                     \
                                                  childTree->node[keepLoc].id;
			childTree->node[keepLoc].id = ctr;
		  } 
		  free(locationStack);
                                        /* by ED algorithm, set best cost */
		  bestCost = childTree->costOfSubTask;
		  flushBestCostStack(TH);
		  freeTree(sourceTree);
		  pushBestCostNode(childTree);
	}
	else {
		  freeTree(childTree);
		  freeTree(sourceTree);
	}
}

struct subTask *genNextEckDayTask(struct subTask *currTask, int taxa, int pos, THREADED) {

   int i;

   struct subTask *nextTask = (struct subTask *)calloc(1, sizeof(struct subTask));
                                        /* get highest internal node */
   nextTask->rhtSubTreeRootIndx = 1;
   nextTask->highestInternalNode = currTask->highestInternalNode-1;
                                       /* get next taxa to be added */
   nextTask->lastTaxaAdded = taxa;
                                       /* get # elements in the preordered 
                                          array of parent task */
   nextTask->sizeOfSubTaskArray = currTask->sizeOfSubTaskArray+2;
   nextTask->node=                                                            \
               (struct subTaskArrayNode *)calloc(nextTask->sizeOfSubTaskArray,\
                sizeof(struct subTaskArrayNode));
                                        /* cp elem 0-pos frm parent to child*/
   for (i = 0; i < pos; i++) {
     nextTask->node[i].id = currTask->node[i].id; 
     nextTask->node[i].costUntilThisNode = currTask->node[i].costUntilThisNode; 
     nextTask->node[i].costAtThisNode = currTask->node[i].costAtThisNode; 
     nextTask->node[i].stateVector =                                          \
                         (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
     memcpy(nextTask->node[i].stateVector, currTask->node[i].stateVector,     \
            matrix.num_pars_inf_sites*sizeof(int));
   }

                                        /* 2 elems frm "pos" get new values 
                                           (internal node, taxa) pair */
   nextTask->node[pos].id = nextTask->highestInternalNode;
   nextTask->node[pos].costUntilThisNode = 0;
   nextTask->node[pos].costAtThisNode = 0;
   nextTask->node[pos].stateVector = (int *)calloc(matrix.num_pars_inf_sites, \
                                     sizeof(int));

   nextTask->node[pos+1].id = nextTask->lastTaxaAdded;
   nextTask->node[pos+1].costUntilThisNode = 0;
   nextTask->node[pos+1].costAtThisNode = 0;
   nextTask->node[pos+1].stateVector=(int *)calloc(matrix.num_pars_inf_sites, \
                                     sizeof(int));
   memcpy(nextTask->node[pos+1].stateVector,                                  \
          matrix.reord_sites_enc[nextTask->lastTaxaAdded-1],                  \
          matrix.num_pars_inf_sites*sizeof(int));
                                       /* cp elem "pos" in parent to "pos+2" 
                                          in child */
   for (i=pos; i < currTask->sizeOfSubTaskArray; i++) {
     nextTask->node[i+2].id = currTask->node[i].id; 
     nextTask->node[i+2].costUntilThisNode =                                  \
                                     currTask->node[i].costUntilThisNode; 
     nextTask->node[i+2].costAtThisNode = currTask->node[i].costAtThisNode; 
     nextTask->node[i+2].stateVector=(int *)calloc(matrix.num_pars_inf_sites, \
                                      sizeof(int));
     memcpy(nextTask->node[i+2].stateVector, currTask->node[i].stateVector,   \
                                      matrix.num_pars_inf_sites*sizeof(int));
   }

   nextTask->root.stateVector = (int *)calloc(matrix.num_pars_inf_sites,      \
                                      sizeof(int));
                                        /* set parent child relationship */
   setParentChildIndices(nextTask);
                                        /* compute cost of the new child */
   computeEckDayCost(nextTask, TH);
   return (nextTask);
}

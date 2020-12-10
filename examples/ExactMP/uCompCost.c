/******************************************************************************
* uCompCost.c                                                                 *
* Computes cost in unoptimized way. "u" prefix indicates "unoptimized"        *
* This was the initial version, used when taxa is <= 8                        *
******************************************************************************/
#include "header.h"

void computeCost (struct subTask *sourceSubTask, int prevRhtSubTreeRootIndx, int pos, THREADED) {

  int parentIndx, rChildIndx, lChildIndx, arrIndx, endIndx, strtIndx;
  struct subTaskArrayNode *leftChild, *rightChild, *parent;

                                           /* Start from where the node was 
                                              inserted in parent tree */
  strtIndx = pos;
                                           /* End at the root of the subtree 
                                              (if RHS: -1, if LHS: 0) */
  endIndx=(pos > prevRhtSubTreeRootIndx) ?                                 \
           sourceSubTask->rhtSubTreeRootIndx : 0;

                                           /* Loop until root of the subtree 
                                              is encountered */
  arrIndx = strtIndx;
  while (arrIndx >= endIndx) {
                                           /* For all internal nodes which are 
                                              also parent nodes, */
	if (sourceSubTask->node[arrIndx].id < 0) {
                                           /* Gather children information */
		lChildIndx = sourceSubTask->node[arrIndx].lChildIndx;
		leftChild = &sourceSubTask->node[lChildIndx];

		rChildIndx = sourceSubTask->node[arrIndx].rChildIndx;
		rightChild = &sourceSubTask->node[rChildIndx];

		parent = &sourceSubTask->node[arrIndx];
                                           /* Send left child, parent, right 
                                              child for Fitch's operation */
		performFitchOp(leftChild , parent, rightChild, TH);
                                           /* If didnt reach the root of 
                                              left/right subtree, jump to the
                                              parent of current node */
		if(arrIndx > 0 && arrIndx != sourceSubTask->rhtSubTreeRootIndx){
		  parentIndx = sourceSubTask->node[arrIndx].parentIndx;
		  arrIndx = parentIndx;
		}
                                           /* If reached the root of 
                                              left/right subtree, break */
		else {
		  break;
		}
	}
  }
                                           /* Both subtrees (left & right) are 
                                              computed and cost kept
                                              in the roots of subtrees, compute
                                              cost of root node now */
  sourceSubTask->root.rChildIndx = sourceSubTask->rhtSubTreeRootIndx;
  sourceSubTask->root.lChildIndx = 0;
  performFitchOp(&sourceSubTask->node[0],&sourceSubTask->root,               \
                 &sourceSubTask->node[sourceSubTask->rhtSubTreeRootIndx],TH);
  sourceSubTask->costOfSubTask = sourceSubTask->root.costUntilThisNode;
}

void performFitchOp(struct subTaskArrayNode *lChild,                         \
                    struct subTaskArrayNode *parent,                         \
                    struct subTaskArrayNode *rChild, THREADED) {
   int site = 0;
   int length = 0;

#ifdef FITCH_OPERATION_VERBOSE
   /*
   printf("\n\t\t\t\t<><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
   printf("\t\t\t\t\t          FITCH OPERATION         \n");
   printf("\t\t\t\t<><><><><><><><><><><><><><><><><><><><><><><><><><>");
   printf("\n\t\t\t\t\tLeft child: %d, Parent: %d, Right Child: %d\n",      \
                                  lChild->id,parent->id,rChild->id);
   */
   printf("\t\t\t\t\t  Fitch Op {L:P:R}: {%4d:%4d:%4d}\n",      \
                                  lChild->id,parent->id,rChild->id);
#endif


   for (site=0; site < matrix.num_pars_inf_sites; site++) {
      if ((lChild->stateVector[site] & rChild->stateVector[site])==0) {
        parent->stateVector[site] =                                          \
                        lChild->stateVector[site] | rChild->stateVector[site];
	length = length + 1;
      }
      else {
        parent->stateVector[site] =                                          \
                        lChild->stateVector[site] & rChild->stateVector[site];
      }
   }
   parent->costAtThisNode = length;
   parent->costUntilThisNode = lChild->costUntilThisNode + rChild->costUntilThisNode + length;
}

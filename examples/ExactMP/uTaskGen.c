/******************************************************************************
* uTaskGen.c                                                                  *
* Generates task in unoptimized way. "u" prefix indicates "unoptimized"       *
* This was the initial version, used when taxa is <= 8                        *
******************************************************************************/
#include "header.h"

struct subTask *genNextTask(struct subTask *currTask, int taxa, int pos, THREADED) {

   int i;

   struct subTask *nextTask = (struct subTask *)calloc(1, sizeof(struct subTask));

                                  /* set chld task's pos of int node "-1" ,
                                     this is rhs subtree root
                                     if new node inserted after "-1" 
                                     then pos of "-1" in child is 
                                     previous + 2, else it stays same */
   if (pos < currTask->rhtSubTreeRootIndx) {
      nextTask->rhtSubTreeRootIndx = currTask->rhtSubTreeRootIndx+2;
   }
   else {
      nextTask->rhtSubTreeRootIndx = currTask->rhtSubTreeRootIndx;
   }

                                  /* get highest internal node of the 
                                     child task */
   nextTask->highestInternalNode = currTask->highestInternalNode-1;
                                  /* get next taxa to be added */
   nextTask->lastTaxaAdded = taxa;
                                  /* get # elements in the preordered 
                                     array of parent task */
   nextTask->sizeOfSubTaskArray = currTask->sizeOfSubTaskArray+2;

   nextTask->node =                                                    \
                   (struct subTaskArrayNode *)
                   calloc(nextTask->sizeOfSubTaskArray,                \
                   sizeof(struct subTaskArrayNode));

                                   /* copy all elements from 0 to pos 
                                      from parent to child task */
   for (i = 0; i < pos; i++) {
     nextTask->node[i].id = currTask->node[i].id; 
     nextTask->node[i].costUntilThisNode = currTask->node[i].costUntilThisNode; 
     nextTask->node[i].costAtThisNode = currTask->node[i].costAtThisNode; 
     nextTask->node[i].stateVector =                                    \
                   (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
     memcpy(nextTask->node[i].stateVector,                              \
          currTask->node[i].stateVector, matrix.num_pars_inf_sites*sizeof(int));
   }

                                   /* 2 elems from "pos" get new values 
                                      (internal node, taxa) pair */
   nextTask->node[pos].id = nextTask->highestInternalNode;
   nextTask->node[pos].costUntilThisNode = 0;
   nextTask->node[pos].costAtThisNode = 0;
   nextTask->node[pos].stateVector =                                   \
                        (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));

   nextTask->node[pos+1].id = nextTask->lastTaxaAdded;
   nextTask->node[pos+1].costUntilThisNode = 0;
   nextTask->node[pos+1].costAtThisNode = 0;
   nextTask->node[pos+1].stateVector =                                  \
                       (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
   memcpy(nextTask->node[pos+1].stateVector,                           \
          matrix.reord_sites_enc[nextTask->lastTaxaAdded-1],           \
          matrix.num_pars_inf_sites*sizeof(int));
                                    /* copy all elements from "pos" 
                                       in parent to "pos+2" in child */
   for (i=pos; i < currTask->sizeOfSubTaskArray; i++) {
     nextTask->node[i+2].id = currTask->node[i].id; 
     nextTask->node[i+2].costUntilThisNode=currTask->node[i].costUntilThisNode; 
     nextTask->node[i+2].costAtThisNode = currTask->node[i].costAtThisNode; 
     nextTask->node[i+2].stateVector = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
     memcpy(nextTask->node[i+2].stateVector,                            \
         currTask->node[i].stateVector, matrix.num_pars_inf_sites*sizeof(int));
   }

   nextTask->root.stateVector=                                          \
                       (int *)calloc(matrix.num_pars_inf_sites,sizeof(int));
                                    /* set par child relation for child tree*/
   setParentChildIndices(nextTask);
                                    /* compute cost of the new child */
   computeCost(nextTask, currTask->rhtSubTreeRootIndx, pos, TH);
   return (nextTask);

}

void setParentChildIndices (struct subTask *sourceSubTask) {
  int arrIndx;
  int indx;

                                    /* for all nodes in array starting frmRHS */
  for (arrIndx = sourceSubTask->sizeOfSubTaskArray-1; arrIndx >= 0;       \
                                 arrIndx = arrIndx-1) {
                                    /* if node is leaf node, push it to stack */
    if (sourceSubTask->node[arrIndx].id > 0) {
      sourceSubTask->node[arrIndx].rChildIndx = LEAF;
      sourceSubTask->node[arrIndx].lChildIndx = LEAF;
      push(arrIndx);
    }
                                    /* if node is internal node, pop 2 nodes 
                                       from stack, & push this 
                                       node on stack. 2 poped nodes are 
                                       left & right children resp. */
    else if (sourceSubTask->node[arrIndx].id < 0) {
      indx = pop();
      sourceSubTask->node[arrIndx].lChildIndx = indx;
      sourceSubTask->node[indx].parentIndx = arrIndx;

      indx = pop();
      sourceSubTask->node[arrIndx].rChildIndx = indx;
      sourceSubTask->node[indx].parentIndx = arrIndx;
      push(arrIndx);
    }
  }
                                     /* Parent of roots of rhs & lhs 
                                        subtrees is main root */
  sourceSubTask->node[0].parentIndx = ROOT;
  sourceSubTask->node[sourceSubTask->rhtSubTreeRootIndx].parentIndx = ROOT;
  
                                     /* Children of main root are roots 
                                        of rhs & lhs subtrees */
  sourceSubTask->root.lChildIndx = pop();
  sourceSubTask->root.rChildIndx = pop();

}

struct subTask *genIstSubTask(THREADED) {
  struct subTask *currSubTask;
  struct subTaskArrayNode leftChild, rightChild;

  currSubTask = (struct subTask *)calloc(1, sizeof(struct subTask));

  currSubTask->highestInternalNode = -1;
  currSubTask->lastTaxaAdded = 3;
  currSubTask->sizeOfSubTaskArray = 4;

  currSubTask->node =                                           \
       (struct subTaskArrayNode *)calloc(4, sizeof(struct subTaskArrayNode ));

  currSubTask->node[0].id= 1;
  currSubTask->node[0].costUntilThisNode = 0;
  currSubTask->node[0].parentIndx = ROOT;
  currSubTask->node[0].lChildIndx = LEAF;
  currSubTask->node[0].rChildIndx = LEAF;
  currSubTask->node[0].stateVector =                            \
                   (int *)calloc(matrix.num_pars_inf_sites, sizeof (int));
  memcpy(currSubTask->node[0].stateVector,                      \
          matrix.reord_sites_enc[0],matrix.num_pars_inf_sites*sizeof(int));

  currSubTask->node[2].id = 2;
  currSubTask->node[2].costUntilThisNode = 0;
  currSubTask->node[2].parentIndx = 1;
  currSubTask->node[2].lChildIndx = LEAF;
  currSubTask->node[2].rChildIndx = LEAF;
  currSubTask->node[2].stateVector =                            \
                   (int *)calloc(matrix.num_pars_inf_sites, sizeof (int));
  memcpy(currSubTask->node[2].stateVector,matrix.reord_sites_enc[1], \
                   matrix.num_pars_inf_sites*sizeof(int));
  leftChild = currSubTask->node[2];

  currSubTask->node[3].id = 3;
  currSubTask->node[3].costUntilThisNode = 0;
  currSubTask->node[3].parentIndx = 1;
  currSubTask->node[3].lChildIndx = LEAF;
  currSubTask->node[3].rChildIndx = LEAF;
  currSubTask->node[3].stateVector = (int *)calloc(matrix.num_pars_inf_sites, \
                                     sizeof (int));
  memcpy(currSubTask->node[3].stateVector,matrix.reord_sites_enc[2],          \
                                     matrix.num_pars_inf_sites*sizeof(int));
  rightChild = currSubTask->node[3];

  currSubTask->node[1].id = -1;
  currSubTask->node[1].costUntilThisNode = 0;
  currSubTask->node[1].parentIndx = ROOT;
  currSubTask->node[1].lChildIndx = 2;
  currSubTask->node[1].rChildIndx = 3;
  currSubTask->node[1].stateVector = (int *)calloc(matrix.num_pars_inf_sites, \
                                     sizeof (int));
  performFitchOp(&leftChild, &currSubTask->node[1], &rightChild, TH);

  leftChild = currSubTask->node[0];
  rightChild = currSubTask->node[1];
  currSubTask->root.id= 0;
  currSubTask->root.lChildIndx = 0;
  currSubTask->root.rChildIndx = 1;
  currSubTask->root.costUntilThisNode = 0;
  currSubTask->root.stateVector = (int *)calloc(matrix.num_pars_inf_sites,    \
                                     sizeof (int));
  performFitchOp(&leftChild, &currSubTask->root, &rightChild, TH);

  currSubTask->costOfSubTask = currSubTask->root.costUntilThisNode;
  currSubTask->rhtSubTreeRootIndx = 1;

  return(currSubTask);
}

/******************************************************************************
* uDataStr.c                                                                  *
* Data structures for unoptimized code. "u" prefix indicates "unoptimized"    *
* This was the initial version, used when taxa is <= 8                        *
******************************************************************************/
#include "header.h"

void printBestSolutions(THREADED){
        struct bestCostStack *delNode;

	while (bestCostStackHead != NULL) {
                delNode = bestCostStackHead;
                bestCostStackHead=bestCostStackHead->down;
#ifdef BnB_VERBOSE
		printTree(delNode->bCNode, -1, TH);
#endif

		freeTree(delNode->bCNode);
		free(delNode);
	}
        bestCostStackHead = NULL;

}

void pushBestCostNode(struct subTask *addBCNode) {
        struct bestCostStack *newNode;
	++numBestSolutions;
	if (bestCostStackHead == NULL) {
		bestCostStackHead =                                        \
                  (struct bestCostStack *)malloc(sizeof(struct bestCostStack));
		bestCostStackHead->down = NULL;
		bestCostStackHead->bCNode= addBCNode;
	}
	else {
		newNode =                                                  \
                  (struct bestCostStack *)malloc(sizeof(struct bestCostStack));
		newNode->down = bestCostStackHead;
		newNode->bCNode = addBCNode;
		bestCostStackHead = newNode;
	}
}

void flushBestCostStack(THREADED) {
        struct bestCostStack *delNode;
	numBestSolutions=0;
	while (bestCostStackHead != NULL) {
                delNode = bestCostStackHead;
                bestCostStackHead=bestCostStackHead->down;
		freeTree(delNode->bCNode);
		free(delNode);
	}
        bestCostStackHead = NULL;
}

struct taskNode *popTask(void) {
	struct taskStack *temp;
	struct taskNode *retVal;
	if (taskSTop != NULL) {
		temp = taskSTop;
		taskSTop = taskSTop->down;
		retVal = temp->tNode;
		free(temp);
	}
	return (retVal);
}

void pushTask(struct taskNode *srcTNode) {
	struct taskStack *temp;
	if (taskSTop==NULL) {
		temp = (struct taskStack *)malloc(sizeof(struct taskStack));
		temp->down = NULL;
		taskSTop = temp;
	        temp->tNode=(struct taskNode *)malloc(sizeof(struct taskNode));
		temp->tNode->partialTask = srcTNode->partialTask;
		temp->tNode->taxaRemHead = srcTNode->taxaRemHead;
		temp->tNode->taxaRemTail = srcTNode->taxaRemTail;
		temp->tNode->taxaDoneHead = srcTNode->taxaDoneHead;
		temp->tNode->taxaDoneTail = srcTNode->taxaDoneTail;
	}
	else {
		temp = (struct taskStack *)malloc(sizeof(struct taskStack));
		temp->down = taskSTop;
		temp->tNode=(struct taskNode *)malloc(sizeof(struct taskNode));
		temp->tNode->partialTask = srcTNode->partialTask;
		temp->tNode->taxaRemHead = srcTNode->taxaRemHead;
		temp->tNode->taxaRemTail = srcTNode->taxaRemTail;
		temp->tNode->taxaDoneHead = srcTNode->taxaDoneHead;
		temp->tNode->taxaDoneTail = srcTNode->taxaDoneTail;
		taskSTop = temp;
	}
}

void freeStackNode(struct taskNode *srcNode, THREADED) {
   struct taxaQ *temp;
   freeTree(srcNode->partialTask);
   while (srcNode->taxaRemHead != NULL) {
     temp = srcNode->taxaRemHead;
     srcNode->taxaRemHead = srcNode->taxaRemHead->next;
     free(temp);
   }
   while (srcNode->taxaDoneHead != NULL) {
     temp = srcNode->taxaDoneHead;
     srcNode->taxaDoneHead = srcNode->taxaDoneHead->next;
     free(temp);
   }
   free(srcNode);
}

int dequeueTaxaQ(struct taskNode *srcNode) {
   struct taxaQ *temp;
   int retVal;

   retVal = -1;
   if (srcNode->taxaRemHead != NULL) {
     retVal = srcNode->taxaRemHead->taxa;
     temp = srcNode->taxaRemHead;
     srcNode->taxaRemHead = srcNode->taxaRemHead->next;
     free(temp);
   }
   return (retVal);
}

struct taskNode *copyTaxaQ(struct taskNode *srcNode) {
   struct taskNode *destNode;
   struct taxaQ *temp_node, *tQNode;
   destNode = (struct taskNode *)calloc(1, sizeof(struct taskNode));
   destNode->taxaRemHead = NULL;
   destNode->taxaRemTail = NULL;
   for (tQNode = srcNode->taxaRemHead; tQNode != NULL; tQNode=tQNode->next) {
    if (destNode->taxaRemHead==NULL) {
       destNode->taxaRemHead = (struct taxaQ*)malloc(sizeof(struct taxaQ));
       destNode->taxaRemHead->next = NULL;
       destNode->taxaRemHead->taxa = tQNode->taxa;
       destNode->taxaRemTail = destNode->taxaRemHead;
    }
    else {
        temp_node = (struct taxaQ*)malloc(sizeof(struct taxaQ));
        destNode->taxaRemTail->next = temp_node;
        temp_node->taxa = tQNode->taxa;
        temp_node->next = NULL;
        destNode->taxaRemTail = temp_node;
    }
   }

   destNode->taxaDoneHead = NULL;
   destNode->taxaDoneTail = NULL;
   for (tQNode = srcNode->taxaDoneHead; tQNode != NULL; tQNode=tQNode->next) {
    if (destNode->taxaDoneHead==NULL) {
       destNode->taxaDoneHead = (struct taxaQ*)malloc(sizeof(struct taxaQ));
       destNode->taxaDoneHead->next = NULL;
       destNode->taxaDoneHead->taxa = tQNode->taxa;
       destNode->taxaDoneTail = destNode->taxaDoneHead;
    }
    else {
        temp_node = (struct taxaQ*)malloc(sizeof(struct taxaQ));
	assert(temp_node);
        destNode->taxaDoneTail->next = temp_node;
        temp_node->taxa = tQNode->taxa;
        temp_node->next = NULL;
        destNode->taxaDoneTail = temp_node;
    }
   }
   return (destNode);
}

void enqueueTaxa(void)  {
  struct taxaQ* temp_node;
  int i = 0;
  for (i=4; i <= matrix.num_taxa; i++) {
    if (taxaQHead==NULL) {
        taxaQHead = (struct taxaQ*)malloc(sizeof(struct taxaQ));
        taxaQHead->next = NULL;
        taxaQHead->taxa = i;
        taxaQTail = taxaQHead;
    }
    else {
        temp_node = (struct taxaQ*)malloc(sizeof(struct taxaQ));
        taxaQTail->next = temp_node;
        temp_node->taxa = i;
        temp_node->next = NULL;
        taxaQTail = temp_node;
    }
  }
}

int pop(void) {
	struct integerStack *temp;
	int retVal;
	if (stackTop != NULL) {
		temp = stackTop;
		stackTop = stackTop->down;
		retVal = temp->id;
		free(temp);
	}
	return (retVal);
}

void push(int i) {
	struct integerStack *temp;
	if (stackTop==NULL) {
		temp=(struct integerStack *)malloc(sizeof(struct integerStack));
		temp->down = NULL;
		stackTop = temp;
		temp->id = i;
	}
	else {
		temp=(struct integerStack *)malloc(sizeof(struct integerStack));
		temp->down = stackTop;
		temp->id = i;
		stackTop = temp;
	}
}

void enqueueAddOrdTaxa(int taxa) {
  struct taxaQ *temp_node;
    if (tAddOrdQHead==NULL) {
      tAddOrdQHead = (struct taxaQ*)malloc(sizeof(struct taxaQ));
      tAddOrdQHead->next = NULL;
      tAddOrdQHead->taxa = taxa;
      tAddOrdQTail = tAddOrdQHead;
    }
    else {
      temp_node = (struct taxaQ*)malloc(sizeof(struct taxaQ));
      tAddOrdQTail->next = temp_node;
      temp_node->taxa = taxa;
      temp_node->next = NULL;
      tAddOrdQTail = temp_node;
    }
}

void deleteTaxa(int taxa) {
  struct taxaQ *preDelNode = NULL;
  struct taxaQ* delNode = NULL;
  if (taxaQHead!=NULL) {
     if (taxaQHead->taxa == taxa) {
       delNode = taxaQHead->next;
       free(taxaQHead);
       taxaQHead = delNode;
     }
     else {
       for (preDelNode=taxaQHead;preDelNode!=NULL;preDelNode=preDelNode->next){
         if (preDelNode->next->taxa==taxa) {
           delNode = preDelNode->next;
           preDelNode->next = delNode->next;
           free(delNode);
	   break;
         }
       }
     }
  }
}

void freeTree(struct subTask *sourceTree) {
   int i;
   for (i=0; i < sourceTree->sizeOfSubTaskArray; i++) {
      free(sourceTree->node[i].stateVector);
   }
   free(sourceTree->node);
   free(sourceTree->root.stateVector);
   free(sourceTree);
}

void clearMemory(void) {
   int col_v1, row_v1;
   for (col_v1=0; col_v1 < matrix.num_taxa; col_v1++) {
      //free(matrix.taxon_table[col_v1]);   // GIVES A SEG FAULT
   }
   //free(matrix.taxon_table);              // GIVES A SEG FAULT
   for (col_v1=0; col_v1 < matrix.num_taxa; col_v1++) {
      free(matrix.taxons[col_v1]);
   }
   free(matrix.taxons);
   for (col_v1=0; col_v1 < matrix.num_taxa; col_v1++) {
      free(matrix.state_encoding[col_v1]);
   }
   free(matrix.state_encoding);
   for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
      free(matrix.unique_states[col_v1]);
   }
   free(matrix.unique_states);
   free(matrix.num_st_at_site);

   for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
      free(matrix.unique_states_reps[col_v1]);
   }
   free(matrix.unique_states_reps);

   for (row_v1=0; row_v1 < matrix.num_taxa; row_v1++) {
      free(matrix.reord_sites_enc[row_v1]);
   }
   free(matrix.reord_sites_enc);
}

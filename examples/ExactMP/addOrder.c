/******************************************************************************
* addOrder.c                                                                  *
* Decides the initial addition order of taxa using max_mini method            *
******************************************************************************/
#include "header.h"

void taxaAddOrder(THREADED) {
	int taxaOver = 3;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;
	struct taxaQ *t;

	sourceTree = genIstSubTask(TH);
                                        /* enQ all rem taxa */
	enqueueTaxa();
                                        /* apply maxmini, & get 1st add order*/
	mMTPWL = maxMini(sourceTree, TH); 
	enqueueAddOrdTaxa(mMTPWL.taxa); 
	taxaOver += 1; 
                                        /* for all remaining taxa */
	while (taxaOver < matrix.num_taxa) { 
#ifdef TAXA_ADD_ORDER_VERBOSE
		printf("\t\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~\n");
		printf("\t\t\t\t\t\tADDING NEXT TAXA ORDER\n");
		printf("\t\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~");
#endif
                                        /* gen next task w\ prev taxa, pos */
		childTree = genNextTask(sourceTree, mMTPWL.taxa, mMTPWL.pos, TH); 
		freeTree(sourceTree);
		sourceTree = childTree;
#ifdef TAXA_ADD_ORDER_VERBOSE
		printTree(sourceTree,mMTPWL.pos, TH);
#endif
		mMTPWL = maxMini(sourceTree, TH);
		enqueueAddOrdTaxa(mMTPWL.taxa); 
		taxaOver += 1; 
	}
                                        /* if cost < best cost then update */
	if (mMTPWL.length < bestCost) {
	  bestCost = mMTPWL.length;
	  flushBestCostStack(TH);
  	  childTree = genNextTask(sourceTree, mMTPWL.taxa, mMTPWL.pos, TH);
	  freeTree(sourceTree);
	  pushBestCostNode(childTree);
	}
	else {
	  freeTree(sourceTree);
	}

        if (verbose) { 
	printf("\t\t\t[BEST SCORE %d]\n",bestCost+pNiCost);
        }
	/*fprintf(fp,"\t\t\t[BEST SCORE %d]\n",bestCost+pNiCost);*/

        if (verbose) {
        printf("\tTaxa addition order: 1 2 3 ");
        }
        /*fprintf(fp,"\tTaxa addition order: 1 2 3 ");*/
	taxaQueue[0] = 1;
	taxaQueue[1] = 2;
	taxaQueue[2] = 3;
        for (taxaOver=3,t=tAddOrdQHead;t!=NULL;t=t->next,taxaOver=taxaOver+1) {
           taxaQueue[taxaOver] = t->taxa;
           if (verbose) {
           printf("%d ",t->taxa);
           }
           /*fprintf(fp,"%d ",t->taxa);*/
        }
        if (verbose) {
        printf("\n");
        }
        /*fprintf(fp,"\n");*/
}

struct maxMiniTaxaPosWLength maxMini(struct subTask *sourceTree, THREADED) {
	struct maxMiniTaxaPosWLength prevmMTPWL, currmMTPWL; 
	struct taxaQ *tQNode;

	int taxaNum;

	prevmMTPWL.length = 0;
	prevmMTPWL.taxa = 0;
	prevmMTPWL.pos = 0;
	currmMTPWL.length = 0; 
	currmMTPWL.taxa = 0; 
	currmMTPWL.pos = 0; 

#ifdef MAX_MINI_VERBOSE
	printf("\n\t\t******************************************************");
        printf("*******************************************\n");
	printf("\n\t\t\t\t\t\t\tMAXMINI FUNCTION\n");
	printf("\n\t\t******************************************************");
        printf("*******************************************\n");
#endif
                                        /* for all remaining taxa */
	for (tQNode = taxaQHead; tQNode != NULL; tQNode=tQNode->next) { 
		taxaNum = tQNode->taxa; 
                                        /* get tree with minimum length */
#ifdef MAXMINI
		currmMTPWL = getMinTree(sourceTree, taxaNum, TH); 
#endif
#ifdef MAX_MINI_VERBOSE
		printf("\n\t\t->->->->->->RETAIN MAXIMUM LENGTH OF ALL ");
                printf("MINIMUM LENGTHS<-<-<-<-<-<-\n");
                printf("\t\tCurr {MinLen:%3d ",currmMTPWL.length);
                printf("@ Pos:%3d ",currmMTPWL.pos);
                printf(" w/ Taxa:%3d} ",currmMTPWL.taxa);
                printf(" vs. Prev  {MinLen:%3d ",prevmMTPWL.length);
                printf(" @ Pos:%3d ",prevmMTPWL.pos);
                printf("w/ Taxa:%3d}\n",prevmMTPWL.taxa);
#endif
                                        /* length > prev length, keep it */
		if (currmMTPWL.length > prevmMTPWL.length) { 
#ifdef MAX_MINI_VERBOSE
			printf("\t\tCurr Min Lgt (%d) > Prev Min Lgt (%d)=> ",\
                               currmMTPWL.length,prevmMTPWL.length);
#endif
			prevmMTPWL.pos = currmMTPWL.pos; 
			prevmMTPWL.taxa = currmMTPWL.taxa; 
			prevmMTPWL.length = currmMTPWL.length; 
#ifdef MAX_MINI_VERBOSE
                        printf("Update Curr Max Len to ");
                        printf("{MaxLen:%3d @ Pos:%3d w/ Taxa:%3d}\n",       \
                         prevmMTPWL.length,prevmMTPWL.pos,prevmMTPWL.taxa);
#endif
		}
	}
                                        /* remove taxa added */
	deleteTaxa(prevmMTPWL.taxa);
	 
	return (prevmMTPWL); 
} 
 
struct maxMiniTaxaPosWLength getMinTree(struct subTask *sourceTree,          \
                                        int taxaNum, THREADED) {
	struct subTask *childTree; 
	int currentLength = 0; 
	int prevLength = 0; 
	 
	struct maxMiniTaxaPosWLength mMTPWL; 
	 
	int pos;
	int currMinPos; 
	int currMinTaxa; 
	int arraySize = sourceTree->sizeOfSubTaskArray; 
	 
#ifdef GET_MIN_TREE_VERBOSE
	printf("\n\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        printf("~~~~~~~~~~\n");
	printf("\t\t\t\t CHECKING MINIMUM TREE LENGTH FOR ALL POSITIONS FOR ");
        printf("TAXA: %d\n",taxaNum);
	printf("\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        printf("~~~~~~~~\n");
#endif
                                        /* gen 1st case, add taxa in pos0 */
	childTree = genNextTask(sourceTree, taxaNum, 0, TH); 
#ifdef GET_MIN_TREE_VERBOSE
	printTree(childTree, 0, TH);
#endif
                                        /* set initial tree length */
	prevLength = childTree->costOfSubTask;
	freeTree(childTree);

	currMinPos = 0;
	currMinTaxa = taxaNum; 

                                        /* for all remaining positions */
	for (pos=1; pos < arraySize; pos++) { 
		if (sourceTree->node[pos].id != -1) { 
			childTree = genNextTask(sourceTree,taxaNum,pos, TH); 
#ifdef GET_MIN_TREE_VERBOSE
	                printTree(childTree, pos, TH);
#endif
			currentLength = childTree->costOfSubTask; 
#ifdef GET_MIN_TREE_VERBOSE
			printf("\n\t\t->->->->->->RETAIN MINIMUM LENGTH FROM ");
                        printf("CURRENT & PREVIOUS LENGTHS<-<-<-<-<-<-\n");
                        printf("\n\t\t\tCurr {MinLen:%3d ",currentLength);
                        printf(" @ Pos:%3d w/",pos);
                        printf(" Taxa:%3d} vs.",taxaNum);
                        printf(" Prev {MinLen:%3d @",prevLength);
                        printf(" Pos:%3d w/",currMinPos);
                        printf("Taxa:%3d}\n\n",currMinTaxa);
#endif
                                        /* if length < prev length, keep it */
			if (currentLength < prevLength) {
#ifdef GET_MIN_TREE_VERBOSE
				printf("\t\tCurr Lgt (%d) < Prev Lgt (%d)=> ",\
                                currentLength,prevLength);
#endif
				prevLength = currentLength; 
				currMinPos = pos; 
				currMinTaxa = taxaNum; 
#ifdef GET_MIN_TREE_VERBOSE
			       printf("Update Curr Lgt to: %d ",currentLength);
                               printf(" & keeping Taxa: %d", currMinTaxa);
                               printf(" Pos: %d\n\n",currMinPos);
#endif
			} 
			freeTree(childTree);
		}
	}
		 
                                        /* return position & tree */
	mMTPWL.taxa = currMinTaxa;
	mMTPWL.pos = currMinPos;
	mMTPWL.length = prevLength;
		 
	return (mMTPWL); 
}

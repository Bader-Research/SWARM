/******************************************************************************
* aggGreedy.c                                                                 *
* Computes initial upper bound using Improved Eck & Dayhoff Greedy method     *
******************************************************************************/
#include "header.h"

void getIstBestCost(THREADED) {
	int taxaOver = 3;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

#ifdef GET_IST_BESTCOST_VERBOSE
	printf("\n\t\t\t\t\t$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\t\t\t\t\t    IN GET IST BEST COST FUNCTION\n");
	printf("\t\t\t\t\t$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
#endif

	sourceTree = genIstSubTask(TH);
#ifdef GET_IST_BESTCOST_VERBOSE
	printTree(sourceTree,-1, TH);
#endif

	enqueueTaxa();
	mMTPWL = smartGreedy(sourceTree, TH); 
	taxaOver += 1; 

	while (taxaOver < matrix.num_taxa) {
#ifdef GET_IST_BESTCOST_VERBOSE
		printf("\n\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		printf("\t\t\t\t\tADDING NEXT TAXA IN GET IST BEST COST\n");
		printf("\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
#endif
		childTree = genNextTask(sourceTree, mMTPWL.taxa, mMTPWL.pos, TH); 
		freeTree(sourceTree);
		sourceTree = childTree;
#ifdef GET_IST_BESTCOST_VERBOSE
		printTree(sourceTree,mMTPWL.pos, TH);
#endif
		mMTPWL = smartGreedy(sourceTree, TH);
		taxaOver += 1; 
	}
	childTree=genNextTask(sourceTree,mMTPWL.taxa,mMTPWL.pos,TH);
	freeTree(sourceTree);
	bestCost = mMTPWL.length;
	pushBestCostNode(childTree);
}

struct maxMiniTaxaPosWLength smartGreedy(struct subTask *sourceTree,THREADED) { 
	struct maxMiniTaxaPosWLength prevmMTPWL, currmMTPWL; 
	struct taxaQ *tQNode;

	int taxaNum;

	prevmMTPWL.length = 1000000000;
	prevmMTPWL.taxa = 0;
	prevmMTPWL.pos = 0;
	currmMTPWL.length = 0; 
	currmMTPWL.taxa = 0; 
	currmMTPWL.pos = 0; 

#ifdef SMART_GREEDY_VERBOSE
	printf("\n\t*******************************************************");
        printf("******************************************\n");
        fflush(stdout);
	printf("\n\t\t\t\t\t\tSMART GREEDY FUNCTION\n");
        fflush(stdout);
	printf("\n\t*******************************************************");
        printf("******************************************\n");
        fflush(stdout);
#endif

	for (tQNode = taxaQHead; tQNode != NULL; tQNode=tQNode->next) { 
		taxaNum = tQNode->taxa; 
		currmMTPWL = getMinTree(sourceTree, taxaNum, TH); 
#ifdef SMART_GREEDY_VERBOSE
		printf("\n\t\t->->->->->->RETAIN MINIMUM LENGTH OF ALL ");
                printf("MINIMUM LENGTHS<-<-<-<-<-<-\n");
        	fflush(stdout);
		printf("\t\tRecd {MinLen:%3d @ ",currmMTPWL.length);
		printf("Pos:%3d w/ Taxa:%3d} ",currmMTPWL.pos, taxaNum);
                printf(" vs. Prev {MinLen:%3d ",prevmMTPWL.length);
                printf(" @ Pos:%3d w/",prevmMTPWL.pos);
                printf(" Taxa:%3d}\n\n",prevmMTPWL.taxa);
        	fflush(stdout);
#endif
		if (currmMTPWL.length < prevmMTPWL.length) { 
#ifdef SMART_GREEDY_VERBOSE
			printf("\tCurr Min Lgt (%d) < Prev Min Lgt (%d)=> ", \
                                currmMTPWL.length,prevmMTPWL.length);
        		fflush(stdout);
#endif
			prevmMTPWL.pos = currmMTPWL.pos; 
			prevmMTPWL.taxa = currmMTPWL.taxa; 
			prevmMTPWL.length = currmMTPWL.length; 
#ifdef SMART_GREEDY_VERBOSE
			printf("Updating: {CurrMinLen:%3d ",prevmMTPWL.length);
			printf(" @ Pos:%3d w/",prevmMTPWL.pos);
			printf(" Taxon: %3d}\n\n",prevmMTPWL.taxa);
        		fflush(stdout);
#endif
		}
	}
	deleteTaxa(prevmMTPWL.taxa); 
	return (prevmMTPWL); 
} 

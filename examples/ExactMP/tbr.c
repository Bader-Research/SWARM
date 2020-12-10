/******************************************************************************
* tbr.c                                                                       *
* TBRs and input tree                                                         *
******************************************************************************/

#include "header.h"

struct tbrNode {                          /* Represents branch of source tree */
   int pos;                               /* pos of start vertex of branch */
   int size;                              /* size of start vertex of branch 
                                             (# children including itself) */
   int lcs;                               /* left child size of start vertex of 
                                             branch */
   int rcs;                               /* size of right child of start 
                                             vertex of branch */
};                                        /* (size of a node = 
                                             # children + 1 for itself) */

struct edgeStack {                        /* stack holding branches of source 
                                             tree for TBRing */
   int frm;
   int to;
   int lr;
   struct edgeStack *down;
};

                                          /* top of stack holds TBR edges */
struct edgeStack *eStop=NULL, *eSend=NULL;
                                          /* global arrays */
int *locStack;                            /* holds loc of node in src array */
int **svStack;                            /* 2D array holding state vector */
int stPtr;                                /* stack pointer */
int *lChildSv;                            /* state vector of left child */
int *rChildSv;                            /* state vector of right child */
int *parentSv;                            /* state vector of parent */

struct edgeStack popEdge(THREADED);
void pushEdge(int , int , int , THREADED);

int tbrCost(int *, int, THREADED);

void tbr(struct subTask *src, THREADED) {
   struct tbrNode *tbrNd;                  /* stores int nodes of source */
   struct edgeStack edge;

   int lci;                                /* left child index */
   int rci;                                /* right child index */
   int pi;                                 /* parent's index */
   int ctr=0, ctr2=0;                      /* temporary counter variables */
   int sArrSize, bArrSize;                 /* src & branch arr sizes for TBR */

   int *tbrSrc;                            /* TBR src arr containing only IDs */
   int *tbrBra;                            /* TBR branch array (that is to be 
                                              swapped) containing only IDs */
   int *tbrMain;                           /* TBR tree with branch swapped */
   int saCtr, baCtr;                       /* counters, TBR src & branch arr */
   int preCost = 1000053612, cost;         /* cost of TBRed tree, 
                                              initialize to very large value */
  
                                           /* initialize space for pointers */
   tbrMain = (int*)calloc(src->sizeOfSubTaskArray, sizeof(int));
      assert(tbrMain);
   tbrBest = (int*)calloc(src->sizeOfSubTaskArray, sizeof(int));
      assert(tbrBest);
   lChildSv = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
      assert(lChildSv);
   rChildSv = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
      assert(rChildSv);
   parentSv = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
      assert(parentSv);
   locStack = (int *)calloc(src->sizeOfSubTaskArray, sizeof(int));
      assert(locStack);
   svStack = (int **)calloc(src->sizeOfSubTaskArray, sizeof(int *));
      assert(svStack);
   for (ctr=0; ctr < src->sizeOfSubTaskArray; ctr++) {
	svStack[ctr] = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
   	   assert(svStack[ctr]);
   }

#ifdef TBR_VERBOSE
   printf("\t\t\t!@#$%Received tree to be TBR-ed!@#$%\n");
   printTree(src, -1, TH);
#endif

                                            /* allocate space for tbr array 
                                               equal to # internal nodes 
					       'cause, this stores information 
                                               of internal nodes only */
   tbrNd = (struct tbrNode *)calloc(abs(src->highestInternalNode),           \
           sizeof(struct tbrNode));
   assert(tbrNd);
                                            /* for all nodes (int & ext) 
                                               in the source array */
					    /* excluding root of LHS subtree 
                                               get the node starting from RHS */
   for (ctr = src->sizeOfSubTaskArray-1; ctr > 0; ctr--) {
#ifdef TBR_VERY_VERBOSE
	   printf("Current node id in the array: %d ",src->node[ctr].id);
#endif
					    /* if node id > 0, identify which 
                                               child it is of parent (left */
				   	    /* or right) and update the size 
                                               of that child in parent to be */
				 	    /* 1 (lcs or rcs) */
	   if (src->node[ctr].id > 0) {
#ifdef TBR_VERY_VERBOSE
		   printf(">0\n");
#endif
	      pi = src->node[ctr].parentIndx;
	      lci = src->node[pi].lChildIndx;
	      rci = src->node[pi].rChildIndx;
	      if (src->node[lci].id == src->node[ctr].id) {
		      tbrNd[abs(src->node[pi].id)-1].lcs = 1;
	      }
	      else if (src->node[rci].id == src->node[ctr].id) {
		      tbrNd[abs(src->node[pi].id)-1].rcs = 1;
	      }
	   }
					    /* if node id < 0, check if it is 
                                               not root of RHS sub tree */
	   else if (src->node[ctr].id < 0 && src->node[ctr].id != -1) {
#ifdef TBR_VERY_VERBOSE
		   printf("<0 && != -1\n");
#endif
		                            /* compute its size (# children  
                                               including itself) */
	      tbrNd[abs(src->node[ctr].id)-1].pos = ctr;
	      tbrNd[abs(src->node[ctr].id)-1].size=                          \
   tbrNd[abs(src->node[ctr].id)-1].lcs+tbrNd[abs(src->node[ctr].id)-1].rcs + 1;
					    /* identify which child it is of 
					       its parent (leftor right) and 
                                               update the size of that child 
                                               in the parent to be
					       1 (lcs or rcs) */
	      pi = src->node[ctr].parentIndx;
	      lci = src->node[pi].lChildIndx;
	      rci = src->node[pi].rChildIndx;
	      if (src->node[lci].id == src->node[ctr].id) {
		      tbrNd[abs(src->node[pi].id)-1].lcs =                  \
                                         tbrNd[abs(src->node[ctr].id)-1].size;
	      }
	      else if (src->node[rci].id == src->node[ctr].id) {
		      tbrNd[abs(src->node[pi].id)-1].rcs =                  \
                                         tbrNd[abs(src->node[ctr].id)-1].size;
	      }

	                                     /* stack left edge the right */
	      lci = src->node[ctr].lChildIndx;
	      rci = src->node[ctr].rChildIndx;
	      pushEdge (src->node[ctr].id, src->node[lci].id, 0, TH);
#ifdef TBR_VERY_VERBOSE
	      printf("pushing left edge: %d, %d\n",                          \
                                         src->node[ctr].id, src->node[lci].id);
#endif
	      pushEdge (src->node[ctr].id, src->node[rci].id, 1, TH);
#ifdef TBR_VERY_VERBOSE
	      printf("pushing right edge: %d, %d\n",src->node[ctr].id,       \
                                                    src->node[rci].id);
#endif
	   }
                                             /* if node is root of RHSsubtree */
	   else if (src->node[ctr].id == -1) {
#ifdef TBR_VERY_VERBOSE
		   printf(" == -1\n");
#endif
		                             /* compute its size (# children 
                                                including itself) */
	      tbrNd[abs(src->node[ctr].id)-1].pos = ctr;
	      tbrNd[abs(src->node[ctr].id)-1].size=tbrNd[abs(src->node[ctr].id)-1].lcs+tbrNd[abs(src->node[ctr].id)-1].rcs + 1;
	                                     /* push left edge followed by 
                                                right edge on stack */
	      lci = src->node[ctr].lChildIndx;
	      rci = src->node[ctr].rChildIndx;
	      pushEdge (src->node[ctr].id, src->node[lci].id, 0, TH);
#ifdef TBR_VERY_VERBOSE
	      printf("pushing left edge: %d, %d\n",                          \
                                        src->node[ctr].id, src->node[lci].id);
#endif
	      pushEdge (src->node[ctr].id, src->node[rci].id, 1, TH);
#ifdef TBR_VERY_VERBOSE
	      printf("pushing right edge: %d, %d\n",src->node[ctr].id,       \
                                        src->node[rci].id);
#endif
	   }
   }

                                              /* after the whole array is over,
                                                 cover root of LHS subtree */
#ifdef TBR_VERY_VERBOSE
   printf("current node id in the array: %d \n",src->node[ctr].id);
#endif
   ctr = 0;
   if (src->node[ctr].id < 0) {
		                              /* compute its size (# children 
                                                 including itself) */
     tbrNd[abs(src->node[ctr].id)-1].pos = ctr;
     tbrNd[abs(src->node[ctr].id)-1].size=                                   \
   tbrNd[abs(src->node[ctr].id)-1].lcs+tbrNd[abs(src->node[ctr].id)-1].rcs + 1;
 
	                                      /* push left edge followed by 
                                                 right edge on stack */
     lci = src->node[ctr].lChildIndx;
     rci = src->node[ctr].rChildIndx;
     pushEdge (src->node[ctr].id, src->node[lci].id, 0, TH);
#ifdef TBR_VERY_VERBOSE
     printf("pushing left edge: %d, %d\n",src->node[ctr].id,src->node[lci].id);
#endif
     pushEdge (src->node[ctr].id, src->node[rci].id, 1, TH);
#ifdef TBR_VERY_VERBOSE
     printf("pushing right edge: %d, %d\n",src->node[ctr].id,src->node[rci].id);
#endif
   }

#ifdef TBR_VERBOSE
   printf("\n\n");
   for (ctr=0; ctr < abs(src->highestInternalNode); ctr++) {
     printf("INODE: %d --> INODEPOS: %d, LCS: %d, RCS: %d, INODESIZE: %d\n", \
      -ctr-1, tbrNd[ctr].pos, tbrNd[ctr].lcs, tbrNd[ctr].rcs, tbrNd[ctr].size);
   }
   printf("\n\n");
#endif

                                          /* having identified all edges that 
                                             are to be TBR-ed */
					  /* for each edge, cut the branch and 
                                             place in all feasible valid */
					  /* locations in the source tree and 
                                             compute the cost */
   edge = popEdge(TH);
                                          /* check if edges are not over */
   while (edge.frm !=0 && edge.to != 0) {
	sArrSize = 0;
	bArrSize = 0;
	saCtr = 0;
	baCtr = 0;
	                                  /* if it is right edge */
        if (edge.lr) {
		                          /* source array size (array excluding
                                             branch array) equals 0 to */
					  /* start of the branch, again from 
                                             end of branch to end of array */
		sArrSize = tbrNd[abs(edge.frm)-1].pos +                     \
                           tbrNd[abs(edge.frm)-1].lcs +                     \
                           src->sizeOfSubTaskArray -                    \
                           (tbrNd[abs(edge.frm)-1].pos +                    \
                           tbrNd[abs(edge.frm)-1].lcs +                    \
                           tbrNd[abs(edge.frm)-1].rcs+1);

		                          /* branch array size is equal to 
                                             size of branch which is the size */
					  /* of the right child of the 
                                             originating node of the branch */
		bArrSize = 1 + tbrNd[abs(edge.frm)-1].rcs;

		                          /* allocate memory for source and 
                                             branch arrays */
		tbrSrc = (int *)calloc(sArrSize, sizeof(int));
		tbrBra = (int *)calloc(bArrSize, sizeof(int));

		                          /* copy the respective portions of 
                                             source into TBR arrays */
#ifdef TBR_VERIFY
		printf("\n\n\t\t\tRight edge (%d ---> %d)",edge.frm, edge.to);
		printf("\nTBR Source array (size: %d): ",sArrSize);
#endif
		if (tbrNd[abs(edge.frm)-1].pos!=0) {
			for (ctr=0; ctr<tbrNd[abs(edge.frm)-1].pos; ctr++) {
				tbrSrc[saCtr] = src->node[ctr].id;
				saCtr = saCtr + 1;
#ifdef TBR_VERIFY
				printf("%d ",src->node[ctr].id);
#endif
			}
		}
		for (ctr=tbrNd[abs(edge.frm)-1].pos+1;                        \
           ctr<tbrNd[abs(edge.frm)-1].pos+tbrNd[abs(edge.frm)-1].lcs+1; ctr++){
			tbrSrc[saCtr] = src->node[ctr].id;
			saCtr = saCtr + 1;
#ifdef TBR_VERIFY
			printf("%d ",src->node[ctr].id);
#endif
		}
		for (ctr=tbrNd[abs(edge.frm)-1].pos + \
                     tbrNd[abs(edge.frm)-1].lcs + \
                     tbrNd[abs(edge.frm)-1].rcs+1;\
                     ctr < src->sizeOfSubTaskArray; ctr++) {
			tbrSrc[saCtr] = src->node[ctr].id;
			saCtr = saCtr + 1;
#ifdef TBR_VERIFY
			printf("%d ",src->node[ctr].id);
#endif
		}

#ifdef TBR_VERIFY
		printf("\nTBR Branch array (size: %d): ",bArrSize);
#endif
		ctr=tbrNd[abs(edge.frm)-1].pos;
		tbrBra[baCtr] = src->node[ctr].id;
		baCtr = baCtr + 1;
#ifdef TBR_VERIFY
		printf("%d ",src->node[ctr].id);
#endif
		for (ctr=tbrNd[abs(edge.frm)-1].pos+                         \
                     tbrNd[abs(edge.frm)-1].lcs+1;                           \
                     ctr < tbrNd[abs(edge.frm)-1].pos+                       \
                     tbrNd[abs(edge.frm)-1].lcs+1+tbrNd[abs(edge.frm)-1].rcs;\
                     ctr++) {
			tbrBra[baCtr] = src->node[ctr].id;
			baCtr = baCtr + 1;
#ifdef TBR_VERIFY
			printf("%d ",src->node[ctr].id);
#endif
		}
	}
	                                 /* if it is left edge, repeat process 
                                            followed for right edge */
	else {
#ifdef TBR_VERIFY
	printf("\n\n\t\t\tLeft edge (%d ----> %d)\n",edge.frm, edge.to);
#endif

	bArrSize = tbrNd[abs(edge.frm)-1].lcs+1;
	sArrSize = src->sizeOfSubTaskArray -                               \
                   (tbrNd[abs(edge.frm)-1].pos+                            \
                   tbrNd[abs(edge.frm)-1].lcs+1) +                         \
                   tbrNd[abs(edge.frm)-1].pos;                             \

		tbrSrc = (int *)calloc(sArrSize, sizeof(int));
		tbrBra = (int *)calloc(bArrSize, sizeof(int));

#ifdef TBR_VERIFY
		printf("TBR Branch array (size: %d): ",bArrSize);
#endif
		for (ctr=tbrNd[abs(edge.frm)-1].pos;                 \
                     ctr <= tbrNd[abs(edge.frm)-1].pos+              \
                     tbrNd[abs(edge.frm)-1].lcs; ctr++) {            \
			tbrBra[baCtr] = src->node[ctr].id;
			baCtr = baCtr + 1;
#ifdef TBR_VERIFY
			printf("%d ",src->node[ctr].id);
#endif
		}

#ifdef TBR_VERIFY
		printf("\nTBR Source array (size: %d): ",sArrSize);
#endif
		if (tbrNd[abs(edge.frm)-1].pos==0) {
			for (ctr=tbrNd[abs(edge.frm)-1].pos+            \
                             tbrNd[abs(edge.frm)-1].lcs+1;              \
                             ctr < src->sizeOfSubTaskArray; ctr++) {
				tbrSrc[saCtr] = src->node[ctr].id;
				saCtr = saCtr + 1;
#ifdef TBR_VERIFY
				printf("%d ",src->node[ctr].id);
#endif
			}
		}
		else if (tbrNd[abs(edge.frm)-1].pos!=0) {
			for (ctr=0; ctr < tbrNd[abs(edge.frm)-1].pos; ctr++) {
				tbrSrc[saCtr] = src->node[ctr].id;
				saCtr = saCtr + 1;
#ifdef TBR_VERIFY
				printf("%d ",src->node[ctr].id);
#endif
			}
			for (ctr=tbrNd[abs(edge.frm)-1].pos+               \
                             tbrNd[abs(edge.frm)-1].lcs+1;                 \
                             ctr<src->sizeOfSubTaskArray;ctr++) {
				tbrSrc[saCtr] = src->node[ctr].id;
				saCtr = saCtr + 1;
#ifdef TBR_VERIFY
				printf("%d ",src->node[ctr].id);
#endif
			}
		}
	}

#ifdef TBR_VERBOSE
	printf("\n\nMAIN SOURCE: ");
	for (ctr=0; ctr < src->sizeOfSubTaskArray; ctr++) {
		printf(" %d ",src->node[ctr].id);
	}
	printf("\nBRANCH: ");
	for (ctr=0; ctr < bArrSize; ctr++) {
		printf("%d ",tbrBra[ctr]);
	}
	printf("\nSOURCE: ");
	for (ctr=0; ctr < sArrSize; ctr++) {
		printf("%d ",tbrSrc[ctr]);
	}
	printf("\n");
#endif
	         
	                                /* TBR branch sorted out, create main 
                                           TBR array by branch swapping */
	for (ctr=0; ctr < sArrSize; ctr++) {
	  memcpy(tbrMain, tbrSrc, ctr*sizeof(int));
	  memcpy(tbrMain+ctr, tbrBra, bArrSize*sizeof(int));
	  memcpy(tbrMain+ctr+bArrSize, tbrSrc+ctr, (sArrSize-ctr)*sizeof(int));

		                        /* compute cost of the branch-swapped 
                                           TBR array */
	  cost = tbrCost(tbrMain, src->sizeOfSubTaskArray, TH);
	  if (cost < preCost) {
	   memcpy(tbrBest,tbrMain,src->sizeOfSubTaskArray*sizeof(int));
	   preCost = cost;
	  }
#ifdef TBR_VERBOSE
	  printf("\nTBR-ed ARRAY #%d with cost %d: ",ctr+1, cost);
	  for (ctr2=0; ctr2 < src->sizeOfSubTaskArray; ctr2++) {
	   printf("%d ",tbrMain[ctr2]);
  	  }
	  printf("\n");
#endif
	}

   	free(tbrSrc);
   	free(tbrBra);
		                         /* keep checking for all branches */
   	edge = popEdge(TH);
   }

   free(tbrNd);
   if (preCost < bestCost) {
	bestCost = preCost;
   }

   free(lChildSv);
   free(rChildSv);
   free(parentSv);
   free(locStack);
   for (ctr=0; ctr < src->sizeOfSubTaskArray; ctr++) {
	free(svStack[ctr]);
   }
   free(svStack);
   free(tbrMain);
   free(tbrBest);
}

                                        /* compute cost of TBR-ed array */
					/* this arr contains only node IDs 
                                           unlike previous forms where */
					/* array contains all information about
                                           the nodes (children) etc */
					/* No stack is implemented, its just 
                                           simulated, so frequent */
					/* system call can be avoided */
int tbrCost(int *tbrMain, int sizeOfArray, THREADED) {
	int ctr, site, lChildIndx, rChildIndx, length;

	stPtr = -1;
	length = 0;

	                                /* starting from RHS of the array, 
                                           go until beginning */
	for (ctr=sizeOfArray-1; ctr >=0; ctr--) {
		                        /* if id > 0, copy state vectors and 
                                           push to stack */
		if (tbrMain[ctr] > 0) {
			stPtr = stPtr + 1;
			locStack[stPtr] = ctr;
			memcpy(svStack[stPtr],                         \
                               matrix.reord_sites_enc[tbrMain[ctr]-1], \
                               matrix.num_pars_inf_sites*sizeof(int));
		}
		                        /* if id < 0, copy state vectors and 
                                           pop twice, compute cost and */
					/* push the negative node onto stack */
		else if (tbrMain[ctr] < 0) {
			lChildIndx = locStack[stPtr];
			memcpy(lChildSv,                                 \
                               svStack[stPtr],                           \
                               matrix.num_pars_inf_sites*sizeof(int));
			stPtr = stPtr - 1;

			rChildIndx = locStack[stPtr];
			memcpy(rChildSv, svStack[stPtr],                 \
                               matrix.num_pars_inf_sites*sizeof(int));
			stPtr = stPtr - 1;

			stPtr = stPtr + 1;
			locStack[stPtr] = ctr;
			for (site=0;site < matrix.num_pars_inf_sites; site++) {
                                         /* If parent needs a union of child 
                                            states, incr length */
                                         /* A intersection B is NULL, there is 
                                            no common state */
			   if ((lChildSv[site] & rChildSv[site])==0) {
			        parentSv[site]=lChildSv[site] | rChildSv[site];
				length = length + 1;
			   }
                                         /* Else, update parent's state to 
                                            be a intersection */
                                         /* A intersection B is not NULL, 
                                            there is common state */
			   else {
			        parentSv[site]=lChildSv[site] & rChildSv[site];
			   }
			}
			memcpy(svStack[stPtr], parentSv,                     \
                               matrix.num_pars_inf_sites*sizeof(int));
		}
	}

	                                 /* compute cost for root nodes of LHS 
                                            & RHS subtrees */
	lChildIndx = locStack[stPtr];
	memcpy(lChildSv, svStack[stPtr],                                     \
               matrix.num_pars_inf_sites*sizeof(int));
	stPtr = stPtr - 1;

	rChildIndx = locStack[stPtr];
	memcpy(rChildSv,svStack[stPtr],matrix.num_pars_inf_sites*sizeof(int));
	stPtr = stPtr - 1;

	stPtr = stPtr + 1;
	locStack[stPtr] = ctr;
	for (site=0; site < matrix.num_pars_inf_sites; site++) {
                                         /* If parent needs a union of child 
                                            states, incr length */
                                         /* A intersection B is NULL, there is 
                                            no common state */
	   if ((lChildSv[site] & rChildSv[site])==0) {
	        parentSv[site] = lChildSv[site] | rChildSv[site];
		length = length + 1;
	   }
                                         /* Else, update parent's state to be 
                                            a intersection */
                                         /* A intersection B is not NULL, 
                                            there is common state */
	   else {
	        parentSv[site] = lChildSv[site] & rChildSv[site];
	   }
	}

	return (length);
}

struct edgeStack popEdge(THREADED) {
	struct edgeStack *temp;
	struct edgeStack retVal;
	retVal.frm = 0;
	retVal.to = 0;
	retVal.lr = 0;
	if (eStop != NULL) {
		temp = eStop;
		eStop = eStop->down;
		retVal.frm = temp->frm;
		retVal.to = temp->to;
		retVal.lr = temp->lr;
		free(temp);
	}

	return (retVal);
}

void pushEdge(int frm, int to, int lr, THREADED) {
	struct edgeStack *temp;
	if (eStop==NULL) {
		temp = (struct edgeStack *)malloc(sizeof(struct edgeStack));
		temp->down = NULL;
		eStop = temp;
		temp->frm = frm;
		temp->to = to;
		temp->lr = lr;
	}
	else {
		temp = (struct edgeStack *)malloc(sizeof(struct edgeStack));
		temp->down = eStop;
		temp->frm = frm;
		temp->to = to;
		temp->lr = lr;
		eStop = temp;
	}
}

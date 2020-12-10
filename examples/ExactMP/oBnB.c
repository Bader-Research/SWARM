#include "header.h"

void bnb(THREADED) {
  int ****workSpace, *plbValues, ***fitchArray, ***mpTree, **wosac;
  int level=0, rsize, lsize, cctr, nctr, taxa, cost, indx, prevTID=0, currTID=0;
  int done, rem, statesOver, statesRem, site, plbLen, plbIndx, mask,
      bpos, maskRem, maskOver;
  int sp=-1, rci=0, lci=0, si=0, length=0,
  	  lc, rc, pa, flag, var, nflg, lim, r, c, ctr, val;
  int lhsSize, rhsSize, row, col, chldCol, pos, mpRows;
  int lhsOffSet, skipSize, *stack;

  skipSize = 5 + matrix.num_pars_inf_sites,
  lhsOffSet = 4;
      
  stack = (int *)calloc(matrix.num_taxa, sizeof(int));
  assert(stack);

  plbValues = (int *)calloc(matrix.num_taxa-6, sizeof(int));
  assert(plbValues);
  for (plbIndx=7; plbIndx < matrix.num_taxa; plbIndx++) {
    for (plbLen=0, site=0; site < matrix.num_pars_inf_sites; site++) {
     for (statesOver=0, done=0; done < plbIndx; done++) {
         statesOver = statesOver |                                             \
                                matrix.reord_sites_enc[taxaQueue[done]-1][site];
     }
     for (statesRem=0, rem=plbIndx; rem < matrix.num_taxa; rem++) {
         statesRem = statesRem |                                               \
                                 matrix.reord_sites_enc[taxaQueue[rem]-1][site];
     }

     mask = 1;
     for (bpos=0; bpos < 32; bpos++) {
         mask = mask << bpos;
         maskRem = statesRem & mask;
         maskOver = statesOver & mask;
         if (maskRem != 0 && maskOver == 0)  {
            ++plbLen;
         }
     }

      plbValues[plbIndx-7] = plbLen;
    }
  }

  fitchArray = (int ***)calloc(16, sizeof(int **));
  assert(fitchArray);
  for (r=0; r < 16; r++) {
     fitchArray[r] = (int **)calloc(16, sizeof(int *));
     assert(fitchArray[r]);
  }
  for (r=0; r < 16; r++) {
     for (c=0; c < 16; c++) {
        fitchArray[r][c] = (int *)calloc(2, sizeof(int));
        assert(fitchArray[r][c]);
     }
  }
  for (r=0; r < 16; r++) {
     for (c=0; c < 16; c++) {
         val = r & c;
         if (val==0) {
                 fitchArray[c][r][0] = r | c;
                 fitchArray[c][r][1] = 1;
         }
         else {
                 fitchArray[c][r][0] = val;
                 fitchArray[c][r][1] = 0;
         }
     }
  }

                                        // allocate memory for the work space
					// this includes all the levels of the 
					// mp tree starting from row with 7 
                                        // taxa. However, other than row 0 in 
                                        // the work space, each row contains
					// only (2*level-5) trees can be 
					// allocated (done to keep memory
					// under control). Row 0 contains all
                                        // (2*7-5)!! = 945 trees
  lim = matrix.num_taxa-6;
  workSpace = (int ****)calloc(lim, sizeof(int ***));
  assert(workSpace);
  workSpace[0] = (int ***)calloc(945, sizeof(int **));
  assert(workSpace[0]);
  for (c=0; c < 945; c++) {
      workSpace[0][c] = (int **)calloc(2, sizeof(int *));
      workSpace[0][c][0] = (int *)calloc(12*skipSize+lhsOffSet, sizeof(int));
      workSpace[0][c][1] = (int *)calloc(12*skipSize+1, sizeof(int));
  }
  for (r=1; r < lim; r++) {
      workSpace[r] = (int ***)calloc((2*(r+7)-5), sizeof(int **));
      for (c=0; c < 2*(r+7)-5; c++) {
          workSpace[r][c] = (int **)calloc(2, sizeof(int *));
          workSpace[r][c][0] = (int *)calloc((2*(r+7)-1)*skipSize+lhsOffSet,   \
			                                           sizeof(int));
          workSpace[r][c][1] = (int *)calloc((2*(r+7)-1)*skipSize+1,           \
			                                           sizeof(int));
      }
  }

  lim = matrix.num_taxa-6;
  wosac = (int **)calloc(lim, sizeof(int *));
  for (r=0; r < lim; r++) {
      wosac[r] = (int *)calloc(2, sizeof(int));
  }

  wosac[0][1] = 945;
  for (r=1; r < lim; r++) {
      wosac[r][1] = 2*(r+7)-5;
  }

  mpRows = 5;
  mpTree = (int ***)calloc(mpRows+1, sizeof(int **));
 
                                        // This piece of code just generates
					// all trees at kth level, where 
					// k=mpRows. This algorithm is very 
					// inefficient though, but it is used 
					// only once in the program...
					// It recursively finds all trees in
					// next level to the current level
					// until the requested level is reached
  for (row=1; row <= mpRows; row++) {
   mpTree[row] = (int **)calloc(numNodesInRow(row, TH), sizeof(int *));
   for (col=0; col < numNodesInRow(row, TH); col++) {
       mpTree[row][col] = (int *)calloc(2*(row+1), sizeof(int ));
   }
  }

  //mpTree[1][0] = (int *)calloc(4, sizeof(int ));
  mpTree[1][0][0] =  1;
  mpTree[1][0][1] = -1;
  mpTree[1][0][2] =  2;
  mpTree[1][0][3] =  3;

  for (row=1; row < mpRows; row++) {
   for (chldCol=0, col=0; col < numNodesInRow(row, TH); col++) {
      for (pos=0; pos < 2*(row+1); pos++) {
         if (mpTree[row][col][pos] != -1) {
            for (ctr=0; ctr < pos; ctr++) {
               mpTree[row+1][chldCol][ctr] = mpTree[row][col][ctr];
	    }
            mpTree[row+1][chldCol][pos] = -row-1;
            mpTree[row+1][chldCol][pos+1] = taxaQueue[row+2];
            for (ctr=pos; ctr < 2*(row+1); ctr++) {
               mpTree[row+1][chldCol][ctr+2] = mpTree[row][col][ctr];
	    }
	    chldCol += 1;
	 }
      }
   } 
  }

					// For all nodes in the current level
					// do the following-
					// (a) set important fields (in offset)
					// (b) copy states of leaf nodes
					// (c) set parent-child relationships
					// (d) compute cost of each of 945 trees
  for (col=0; col < numNodesInRow(mpRows, TH); col++) {
     lhsSize = 0;
     for (pos = 0; pos < 2*(mpRows+1); pos++) {
         if (mpTree[5][col][pos] != -1) {
            ++lhsSize;
         } else {
            chldCol = pos;
            break;
	 }
     }
     rhsSize = 0;
     for (pos = chldCol; pos < 2*(mpRows+1); pos++) {
         ++rhsSize;
     }
                                        // SET IMPORTANT FIELDS OF LHS TREE
                                        // Size of LHS tree
     workSpace[0][col][0][0] = (lhsSize+1)*skipSize + lhsOffSet;  
                                        // Size of RHS tree
     workSpace[0][col][0][1] = rhsSize*skipSize + 1;       
                                        // Highest Internal Node
     workSpace[0][col][0][2] = -5;
                                        // Last taxa added
     workSpace[0][col][0][3] = taxaQueue[5];
                                        // dummy
     workSpace[0][col][0][4] = ROOT;                       
                                        // Next position in taxa Queue
     workSpace[0][col][0][5] = 7;
                                        // cost of task
     workSpace[0][col][0][8] = 0;                          
                                        // parent of lhs root is ROOT
     workSpace[0][col][0][lhsOffSet+skipSize+1] = 4;          
                                        // parent of rhs root is ROOT
     workSpace[0][col][1][2] = 4;                          

                                        // SET IMPORTANT FIELDS OF RHS TREE
                                        // Size of RHS tree
     workSpace[0][col][1][0] = rhsSize*skipSize + 1;       

                                        // Convert 1D array into 2D LHS/RHS
					// array with states loaded
					// This for loop covers LHS tree and 
					// breaks at RHS tree
     for (pos = 0; pos < 2*(mpRows+1); pos++) {
                                        // Get position where RHS tree starts
         if (mpTree[5][col][pos] == -1) {
		 chldCol = pos;
		 break;
	 }

                                        // copy taxon ID
	 workSpace[0][col][0][(pos+1)*skipSize+lhsOffSet] = mpTree[5][col][pos];
                                        // if leaf node, copy state vector
	 if (mpTree[5][col][pos] > 0) {
            memcpy(&workSpace[0][col][0][(pos+1)*skipSize+lhsOffSet+5],        \
			    matrix.reord_sites_enc[mpTree[5][col][pos]-1],     \
			    matrix.num_pars_inf_sites*sizeof(int));
	    workSpace[0][col][0][((pos+1)*skipSize)+lhsOffSet+4] = 0;
	 }
     }

					// This for loop covers RHS tree
     for (row = chldCol; row < 2*(mpRows+1); row++) {
                                        // copy taxon ID
	 workSpace[0][col][1][(row-chldCol)*skipSize+1] = mpTree[5][col][row];
                                        // if leaf node, copy state vector
	 if (mpTree[5][col][row] > 0) {
            memcpy(&workSpace[0][col][1][(row-chldCol)*skipSize+1+5],          \
			    matrix.reord_sites_enc[mpTree[5][col][row]-1],     \
			    matrix.num_pars_inf_sites*sizeof(int));
	    workSpace[0][col][1][(row-chldCol)*skipSize+1+4] = 0;
	 }
     }
  }

                                        // set LHS tree parent child relation
  for (col=0; col < numNodesInRow(mpRows, TH); col++) {
      row = workSpace[0][col][0][0];
      for (sp=-1,ctr=row-skipSize;ctr>=lhsOffSet+skipSize;ctr=ctr-skipSize) {
       if (workSpace[0][col][0][ctr] > 0) {
         stack[++sp] = ctr;
       } else {
         workSpace[0][col][0][stack[sp]+1] = ctr;
	 workSpace[0][col][0][ctr+2] = stack[sp];

	 sp = sp-1;
	 workSpace[0][col][0][stack[sp]+1] = ctr;
	 workSpace[0][col][0][ctr+3] = stack[sp];

	 stack[sp] = ctr;
       }
      }

                                        // set RHS tree parent child relation
      row = workSpace[0][col][1][0];
      for (sp=-1, ctr=row-skipSize; ctr > 0; ctr = ctr-skipSize) {
       if (workSpace[0][col][1][ctr] > 0) {
         stack[++sp] = ctr;
       } else { 
         workSpace[0][col][1][stack[sp]+1] = ctr;
	 workSpace[0][col][1][ctr+2] = stack[sp];

	 sp = sp-1;
	 workSpace[0][col][1][stack[sp]+1] = ctr;
	 workSpace[0][col][1][ctr+3] = stack[sp];

	 stack[sp] = ctr;
       }
      }
      
                                        // compute LHS cost
  for (si=workSpace[0][col][0][0]-skipSize; si > 4; si = si - skipSize) {
        if (workSpace[0][col][0][si] < 0) {
          lci = workSpace[0][col][0][si+2];
          rci = workSpace[0][col][0][si+3];
          for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {
             lc = workSpace[0][col][0][lci+site+5];
	     rc = workSpace[0][col][0][rci+5+site];
             pa = lc | rc;
             var = lc & rc;
             flag = !var;
             nflg = ~flag;
             nflg = nflg + 1;
             pa = pa & nflg;
             pa = pa | var;
             length = length+flag;
             workSpace[0][col][0][si+5+site] = pa;
          }
          workSpace[0][col][0][si+4] = workSpace[0][col][0][lci+4]+workSpace[0][col][0][rci+4]+length;
        }
  }

                                        // compute RHS cost
  for (si=workSpace[0][col][1][0]-skipSize; si >= 1; si = si - skipSize) {
        if (workSpace[0][col][1][si] < 0) {
          lci = workSpace[0][col][1][si+2];
          rci = workSpace[0][col][1][si+3];
          for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {
             lc = workSpace[0][col][1][lci+site+5];
	     rc = workSpace[0][col][1][rci+5+site];
             pa = lc | rc;
             var = lc & rc;
             flag = !var;
             nflg = ~flag;
             nflg = nflg + 1;
             pa = pa & nflg;
             pa = pa | var;
             length = length+flag;
             workSpace[0][col][1][si+5+site] = pa;
          }
          workSpace[0][col][1][si+4] = workSpace[0][col][1][lci+4]+workSpace[0][col][1][rci+4]+length;
        }
  }

                                        // compute ROOT cost
  for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {
     lc = workSpace[0][col][0][lhsOffSet+skipSize+5+site];
     rc = workSpace[0][col][1][6+site];
     pa = lc | rc;
     var = lc & rc;
     flag = !var;
     nflg = ~flag;
     nflg = nflg + 1;
     pa = pa & nflg;
     pa = pa | var;
     length = length+flag;
     workSpace[0][col][0][9+site] = pa;
  }
  workSpace[0][col][0][8] = workSpace[0][col][0][lhsOffSet+skipSize+4]+workSpace[0][col][1][5]+length;
 }

 for (row=1; row < mpRows; row++) {
   for (col=0; col < numNodesInRow(row, TH); col++) {
       free(mpTree[row][col]);
   }
   free(mpTree[row]);
 }
 
 for (col=0; col < numNodesInRow(mpRows, TH); col++) {
       free(mpTree[mpRows][col]);
 }
 free(mpTree[row]);
 free(mpTree);

  /* Actual computation starts here */
   while (1) {
     if (MYTHREAD==0) assignJobs(TH);
     if (request[MYTHREAD]==0) {
       cctr = taskID[MYTHREAD];
       level = 0;
       nctr = 0;
       if (verbose) {
         printf("T%d: Decomposing task %d at level %d, bestCost: %d\n",MYTHREAD, cctr, level, bestCost+pNiCost);
         fflush(stdout);
       }
       lsize = workSpace[0][cctr][0][0];
       rsize = workSpace[0][cctr][0][1];
       taxa = taxaQueue[workSpace[0][cctr][0][5]];
       cost = workSpace[0][cctr][0][8] + plbValues[0];
       loadDistribution[MYTHREAD] += 1;
       if (cost <= bestCost) {
         for (indx=lhsOffSet+skipSize; indx<=lsize-skipSize; indx+=skipSize) {
           genTaskLeft(workSpace[level][cctr],workSpace[level+1][nctr],taxa,indx,stack,TH);
           nctr = nctr + 1;
         }
         for (indx=1+skipSize; indx <= rsize-skipSize; indx += skipSize) {
           genTaskRight(workSpace[level][cctr],workSpace[level+1][nctr],taxa,indx,stack,TH);
           nctr = nctr + 1;
         }
         ++level;
         for (ctr=1; ctr <= matrix.num_taxa-7; ctr++) {
           wosac[ctr][0]=0;
         }
         currTID = 0;
         prevTID = 0;
         while (wosac[1][0] <= wosac[1][1]) {
          prevTID = currTID;
          currTID = wosac[1][0];
          if (prevTID==wosac[1][1] && currTID==0) break;
#ifdef debug
          printf("\nStatus at last level: ");
          fflush(stdout);
          for (ctr=0; ctr < matrix.num_taxa-7; ctr++) {
                printf("wosac[%d][0]: %d   ",ctr, wosac[ctr][0]);
                fflush(stdout);
          }
#endif
          if (level==matrix.num_taxa-7) {
#ifdef debug
             printf("\nwosac[%d][0]:",level);
             fflush(stdout);
#endif
             for (wosac[level][0]=0; wosac[level][0] < wosac[level][1]; wosac[level][0]++) {
#ifdef debug
                printf("  %d",wosac[level][0]);
                fflush(stdout);
#endif
                if (workSpace[level][wosac[level][0]][0][8] < bestCost) {
                   pthread_mutex_lock(&(lock));
                   bestCost = workSpace[level][wosac[level][0]][0][8];
	           printf("\t\t\t\t\t\tT%d: [BEST SCORE: %d]\n",MYTHREAD, bestCost+pNiCost);
                   fflush(stdout);
                                                   // Copy Tree to Solution Queue
                   lsize = workSpace[level][wosac[level][0]][0][0];
                   rsize = workSpace[level][wosac[level][0]][0][1];
                   memcpy(solQ[0][0],workSpace[level][wosac[level][0]][0],lsize*sizeof(int));
                   memcpy(solQ[0][1],workSpace[level][wosac[level][0]][1],rsize*sizeof(int));
	           solQCtr=0;
                   pthread_mutex_unlock(&(lock));
	        }
                else if (workSpace[level][wosac[level][0]][0][8] == bestCost) {
                   pthread_mutex_lock(&(lock));
                   if (solQCtr < keepTrees-1) {
                                                   // Copy Tree to Solution Queue
	             ++solQCtr;
                     lsize = workSpace[level][wosac[level][0]][0][0];
                     rsize = workSpace[level][wosac[level][0]][0][1];
                     memcpy(solQ[solQCtr][0],workSpace[level][wosac[level][0]][0],lsize*sizeof(int));
                     memcpy(solQ[solQCtr][1],workSpace[level][wosac[level][0]][1],rsize*sizeof(int));
                   }
                   pthread_mutex_unlock(&(lock));
	        }
             }
             loadDistribution[MYTHREAD] += wosac[level][1];
#ifdef debug
             printf("\n");
             fflush(stdout);
#endif
          } else {
             nctr = 0;
             cctr = wosac[level][0];
             lsize = workSpace[level][cctr][0][0];
             rsize = workSpace[level][cctr][0][1];
             taxa = taxaQueue[workSpace[level][cctr][0][5]];
             cost = workSpace[level][cctr][0][8] + plbValues[level];
             if (cost <= bestCost) {
                for (indx=lhsOffSet+skipSize; indx<=lsize-skipSize; indx+=skipSize) {
                   genTaskLeft(workSpace[level][cctr],workSpace[level+1][nctr],taxa,indx,stack,TH);
                   nctr = nctr + 1;
                }
                for (indx=1+skipSize; indx <= rsize-skipSize; indx += skipSize) {
                   genTaskRight(workSpace[level][cctr],workSpace[level+1][nctr],taxa,indx,stack,TH);
                   nctr = nctr + 1;
                }
	        wosac[level+1][0] = 0;
	        ++wosac[level][0];
	        ++level;
                loadDistribution[MYTHREAD] += nctr;
	     } else {
	        ++wosac[level][0];
                loadDistribution[MYTHREAD] += 1;
	     }
          }
          for (cctr=level; cctr >=1; cctr--) {
              if (wosac[cctr][0]==wosac[cctr][1]) {
		 --level;
	      } else {
                 level = cctr;
		 break;
	      }
          }
          if (level==0) break;
#ifdef debug
          printf("\n1-l: Wrapping up to level %d, taskID: %d",level, taskID[MYTHREAD]);
          fflush(stdout);
#endif
          if (MYTHREAD==0) assignJobs(TH);
         }
       }
       request[MYTHREAD] = 1;
       if (MYTHREAD==0) assignJobs(TH);
     }
     if (MYTHREAD==0) {
        done = 1;
        if (nextID==945) {
          for (ctr=0; ctr < THREADS; ctr++) {
            if (request[ctr]==0) done=0;
          }
        }
        if (done==0) break;
     } else {
        if (nextID == 945 && request[MYTHREAD]==1) {
#ifdef debug
          printf("T%d breaking\n");
          fflush(stdout);
#endif
          break;
        }
     }
   }
}

void assignJobs(THREADED) {
 int i;
 for (i=0; i < THREADS; i++) {
   if (nextID < 945) {
     if (request[i] == 1) {
       taskID[i] = nextID;
       request[i] = 0;
#ifdef debug
       printf("assigning task %d to thread %d\n",nextID, i);
       fflush(stdout);
#endif
       ++nextID;
     }
   } 
 }
}

int numNodesInRow(int row, THREADED) {
  if (row == 1) return 1;
  else return (2*row-1)*numNodesInRow(row-1, TH);
}

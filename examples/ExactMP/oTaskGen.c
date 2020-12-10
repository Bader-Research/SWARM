#include "header.h"
#define CPU

/*  THIS BLOCK OF LOGIC PROVIDES A FAST METHOD TO COMPUTE AND-OR
 * LC -> LEFT CHILD
 * RC -> RIGHT CHILD
 * PA -> PARENT
 * VAR -> VARIABLE
 * FLAG -> FLAG
 * LENGTH -> LENGTH
 *
 *   pa = lc | rc;
 *   var = lc & rc;
 *   flag = !var;
 *   nflg = ~flag;
 *   nflg = nflg + 1;
 *   pa = pa & nflg;
 *   pa = pa | var;
 *   length = length+flag;
 */

/* THIS BLOCK OF LOGIC PROVIDES A FAST METHOD TO COMPUTE AND-OR
 * New logic for packing 8 sites together, & computing states in one go
 * LC -> LEFT CHILD
 * RC -> RIGHT CHILD
 * PA -> PARENT
 * VAR -> VARIABLE
 * FLAG -> FLAG
 * LENGTH -> LENGTH
 * MASK -> 0001 0001 0001 0001 0001 0001 0001 0001
 * 
 *   1. and = lc & rc
 *   2. or  = lc | rc
 *   3. s1  = and >> 1
 *   4. s2  = s1 >> 1
 *   5. s3  = s2 >> 1
 *   6. s4  = s3 >> 1
 *   7. ones = s1 | s2 
 *   8. ones = ones | s3
 *   9. ones = ones | s4
 *  10. bits = bits & mask
 *  11. bits = ~bits
 *  12. length = 0;
 *      while (bits) {
 *        length++;
 *        bits &= (bits - 1);
 *      }
 *  13. localMask = ones << 1
 *  14. localMask = localMask << 1
 *  15. localMask = localMask << 1
 *  16. localMask = localMask << 1
 *  17. result = or & localMask
 *  18. result = result | and 
 *
 */

/* 
 * OPTIMIZED LOGIC FOR IF THEN ELSE BLOCKS (not applicable always)
 * 
 * TAXON = nextTask[0][ctr]
 * TAXON = TAXON >> 31
 * FLAG = !TAXON
 * 
 * IF TAXON > 0 THEN ABOVE LOGIC GIVES FLAG = 1
 * ELSE IF TAXON < 0 THEN ABOVE LOGIC GIVES FLAG = 0
 *
 */

/*        RHS TREE                       LHS TREE
 * 0 size				0 size of lhs tree
 * 1 id					1 size of rhs tree
 * 2 parent				2 highest internal node
 * 3 left child				3 last taxon added
 * 4 right child			4 ROOT
 * 5 cost				5 pos of taxa in taxaQ to be added next
 * 6 sv0				6 left child
 * 7 sv1				7 right child
 * 8 sv2				8 ROOT COST 
 * 9 sv3				9 sv0
 * 10 sv4				10 
 * 11 sv5				11 
 * 12 sv6				12 
 * 13 sv7				13 
 * 14 sv8				14 
 * 15 sv9				15 
 * 16 id				16 
 * 17 parent				17 
 * 18 left child			18 sv9
 * 19 right child			19 id
 * 20 cost				20 parent
 * 21 sv0				21 left child
 * 30 sv9				22 right child
 * 31 id				23 cost
 * 32 parent                            34 id
 */

int **genIstOptSubTask(int **currSubTask, THREADED) {
  int indx=0, length, skipSize, lhsOffSet, rhsOffSet,
      pa, lc, rc, site, flag, var, lci, rci, nflg;

  skipSize = 5 + matrix.num_pars_inf_sites,
  lhsOffSet = 4;
  rhsOffSet = 1+matrix.num_taxa;

  currSubTask = (int **)calloc(2, sizeof(int *));
  assert(currSubTask);
  currSubTask[0] = (int *)calloc(lhsOffSet+2*skipSize, sizeof(int));
  currSubTask[1] = (int *)calloc(1+2*skipSize+skipSize, sizeof(int));

                                      // Size of LHS including offset
  currSubTask[0][0] = 2*skipSize+lhsOffSet;       
                                      // Size of RHS including 1
  currSubTask[0][1] = 2*skipSize+1+skipSize; 
                                      // Highest Internal node
  currSubTask[0][2] = -1;                    
                                      // Last taxa added
  currSubTask[0][3] = 3;                     
                                      // Not used, dummy
  currSubTask[0][4] = ROOT;                  
                                      // next position in taxa Q to be added
  currSubTask[0][5] = 3;                     
                                      // cost of whole tree
  currSubTask[0][8] = 0;                     
                                      // parent of lhs root is ROOT
  currSubTask[0][lhsOffSet+skipSize+1] = 4;  
                                      // parent of rhs root is ROOT
  currSubTask[1][2] = 4;                     
                                      // skip states of ROOT & move to 1st taxa
				      // ->> 8 places to reach Root->Cost +
				      // ->> matrix.num_pars_inf_sites places
				      //     to skip Root's States +
				      // ->> 1 place to reach to ID of the taxa
  indx = 9+matrix.num_pars_inf_sites;
                                      // Set the ID to be 1
  currSubTask[0][indx] = 1;
                                      // Set cost until this node = 0
  indx = indx+4;
  currSubTask[0][indx] = 0;
  ++indx;
                                      // Starting from State 0, copy all states
  memcpy(&currSubTask[0][indx],matrix.reord_sites_enc[0],                     \
		  matrix.num_pars_inf_sites*sizeof(int));

                                      // Left Tree over, now move to Right tree
                                      // Set size of right tree
  indx = 0;
  currSubTask[1][indx] = 2*skipSize+skipSize+1;
                                      // Set ID of root of right tree
  ++indx;
  currSubTask[1][indx] = -1;
                                      // Rite left-child-location of this root
  indx = indx + 2;
  currSubTask[1][indx] = 1+skipSize;
                                      // Rite right-child-location of this root
  ++indx;
  currSubTask[1][indx] = 1+2*skipSize;
                                      // Set ID of 1st leaf taxon of RHS tree
  indx = indx+2+matrix.num_pars_inf_sites;
  currSubTask[1][indx] = 2;
                                      // Set parent of this leaf node to be the
				      // internal node (-1) at position 1
  ++indx;
  currSubTask[1][indx] = 1;
                                      // Set cost of this leaf node = 0
  currSubTask[1][5+skipSize] = 0;
                                      // Starting frm State 0, copy all states
  indx = indx+4;
  memcpy(&currSubTask[1][indx],matrix.reord_sites_enc[1],                     \
		  matrix.num_pars_inf_sites*sizeof(int));

                                      // Set ID of 2nd leaf taxon of RHS tree
  indx = indx+matrix.num_pars_inf_sites;
  currSubTask[1][indx] = 3;
                                      // Set parent of this leaf node to be the
				      // internal node (-1) at position 1
  ++indx;
  currSubTask[1][indx] = 1;
                                      // Set cost of this leaf node = 0
  indx = indx+3;
  currSubTask[1][indx] = 0;
                                      // Starting from State 0, copy all states
  ++indx;
  memcpy(&currSubTask[1][indx],matrix.reord_sites_enc[2],                      \
		  matrix.num_pars_inf_sites*sizeof(int));

                                      // Compute the cost of RHS tree using
				      // Fitch Operation
  lci = 1+skipSize;
  rci = 1+2*skipSize;
  for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {
           lc = currSubTask[1][lci+site+5];
	   rc = currSubTask[1][rci+5+site];
           pa = lc | rc;
           var = lc & rc;
           flag = !var;
           nflg = ~flag;
           nflg = nflg + 1;
           pa = pa & nflg;
           pa = pa | var;
           length = length+flag;
				      // Store res in corres. parent state place
           currSubTask[1][6+site] = pa;
  }
				      // Final result stored in -1's cost place
  currSubTask[1][5] = length;

                                      // Compute the cost of full tree, after
				      // having computed RHS tree cost using 
				      // Fitch operation
  for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {
     lc = currSubTask[0][lhsOffSet+skipSize+5+site];
     rc = currSubTask[1][6+site];
     pa = lc | rc;
     var = lc & rc;
     flag = !var;
     nflg = ~flag;
     nflg = nflg + 1;
     pa = pa & nflg;
     pa = pa | var;
     length = length+flag;
     currSubTask[0][9+site] = pa;
  }

  currSubTask[0][8] = currSubTask[0][lhsOffSet+skipSize+4] + currSubTask[1][5] \
		                                           + length;
  return(currSubTask);
}

void genTaskLeft(int **currTask,int **nextTask,int taxa,int pos, int *stack, THREADED) {

   int ctsl, ctsr, ntsl, ntsr, skipSize, lhsOffSet, rhsOffSet,
       ctr=0, sp=-1, hin=0,
       rci=0, lci=0, si=0, length=0, 
       site=0, lc, rc, pa, flag, var, nflg;

  skipSize = 5 + matrix.num_pars_inf_sites,
  lhsOffSet = 4;
  rhsOffSet = 1+matrix.num_taxa;

				      // Get size of LHS tree
   ctsl = currTask[0][0];
				      // Get size of RHS tree
   ctsr = currTask[0][1];
				      // Get highest internal node
   hin = currTask[0][2]-1;

				      // Compute length of new LHS tree
      ntsl = ctsl + 2*skipSize;
				      // Compute length of new RHS tree
      ntsr = ctsr;

				      // Write size of LHS tree in location 0
      nextTask[0][0] = ntsl;
				      // Write size of RHS tree in location 1
      nextTask[0][1] = ntsr;
				      // Write Highest Internal Node in loc 2
      nextTask[0][2] = hin;
				      // Write taxa added in location 3
      nextTask[0][3] = taxa;
				      // Write location of next taxa in 
				      // taxaQueue to be added
      nextTask[0][5] = currTask[0][5]+1;
				      // Set parent of LHS root to be ROOT
      nextTask[0][lhsOffSet+skipSize+1] = 4;

				      // Write size of RHS tree in location 0
				      // of RHS tree
      nextTask[1][0] = ntsr;
				      // Set parent of RHS root to be ROOT
      nextTask[1][2] = 4;

                                      // Copy all elements from 0 to pos 
				      // from parent to child task
      memcpy(&nextTask[0][lhsOffSet+skipSize],                                 \
		      &currTask[0][lhsOffSet+skipSize],                        \
		      (pos-lhsOffSet-skipSize)*sizeof(int));

                                      // 2 elems from "pos" get new values 
				      // (internal node, taxa) pair
      nextTask[0][pos] = hin;
      nextTask[0][pos+skipSize] = taxa;
      nextTask[0][pos+skipSize+4] = 0;

				      // Copy states of leaf node newly added
      memcpy(&nextTask[0][pos+skipSize+5],                                     \
		      matrix.reord_sites_enc[taxa-1],                          \
		      matrix.num_pars_inf_sites*sizeof(int));

                                      // For RHS tree, copy all elements from 
				      // 0 to pos from parent to child task
      memcpy(&nextTask[0][pos+2*skipSize], &currTask[0][pos],                  \
		      (ctsl-pos)*sizeof(int));

                                      // since this is for left tree, copy all 
				      // elements from 0 to pos from parent
      memcpy(&nextTask[1][1], &currTask[1][1], (ntsr-1)*sizeof(int));

                                      // set parent child relationship for the 
				      // child tree
      for (sp=-1,ctr=ntsl-skipSize;ctr>=lhsOffSet+skipSize;ctr=ctr-skipSize) {
       if (nextTask[0][ctr] > 0) {
         stack[++sp] = ctr;
       } else {
         nextTask[0][stack[sp]+1] = ctr;
	 nextTask[0][ctr+2] = stack[sp];

	 sp = sp-1;
	 nextTask[0][stack[sp]+1] = ctr;
	 nextTask[0][ctr+3] = stack[sp];

	 stack[sp] = ctr;
       }
      }

                                      // compute cost of left half of the newly
				      // created tree using optimized logic
  si = pos;
  do {
        lci = nextTask[0][si+2];
        rci = nextTask[0][si+3];
        for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {

#ifdef MEMONL
           nextTask[0][si+5+site] = fitchArray[nextTask[0][lci+site+5]][nextTask[0][rci+5+site]][0];
           length = length + fitchArray[nextTask[0][lci+site+5]][nextTask[0][rci+5+site]][1];
#endif

#ifdef MEMVAR
           lc = nextTask[0][lci+site+5];
	   rc = nextTask[0][rci+5+site];
           nextTask[0][si+5+site] = fitchArray[lc][rc][0];
           length = length + fitchArray[lc][rc][1];
#endif

#ifdef CPU
           lc = nextTask[0][lci+site+5];
	   rc = nextTask[0][rci+5+site];
	   pa = lc | rc;
	   var = lc & rc;
	   flag = !var;
           nflg = ~flag;
           nflg = nflg + 1;
	   pa = pa & nflg;
	   pa = pa | var;
	   length = length + flag;
           nextTask[0][si+5+site] = pa;
#endif
        }

                                       /* New Check added */
        /*
        if (nextTask[0][lci+4]+nextTask[0][rci+4]+length > bestCost) {
           nextTask[0][8] = 1000053612;
           goto SKIPALLLEFT;
        }*/

        nextTask[0][si+4] = nextTask[0][lci+4]+nextTask[0][rci+4]+length;
        si = nextTask[0][si+1];
  } while (si!=4);

  for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {

#ifdef MEMONL
     nextTask[0][9+site] = fitchArray[nextTask[0][lhsOffSet+skipSize+5+site]][nextTask[1][6+site]][0];
     length = length + fitchArray[nextTask[0][lhsOffSet+skipSize+5+site]][nextTask[1][6+site]][1];
#endif

#ifdef MEMVAR
     lc = nextTask[0][lhsOffSet+skipSize+5+site];
     rc = nextTask[1][6+site];
     nextTask[0][9+site] = fitchArray[lc][rc][0];
     length = length + fitchArray[lc][rc][1];
#endif
	 
#ifdef CPU
     lc = nextTask[0][lhsOffSet+skipSize+5+site];
     rc = nextTask[1][6+site];
     pa = lc | rc;
     var = lc & rc;
     flag = !var;
     nflg = ~flag;
     nflg = nflg + 1;
     pa = pa & nflg;
     pa = pa | var;
     length = length+flag;
     nextTask[0][9+site] = pa;
#endif
  }

  nextTask[0][8] = nextTask[0][lhsOffSet+skipSize+4]+nextTask[1][5]+length;

/*
SKIPALLLEFT:
  lc = 0;
*/

}

void genTaskRight(int **currTask,int **nextTask,int taxa,int pos, int *stack, THREADED) {

   int ctsl, ctsr, ntsl, ntsr, skipSize, lhsOffSet, rhsOffSet,
       ctr=0, sp=-1, hin=0,
       rci=0, lci=0, si=0, length=0, 
       site=0, lc, rc, pa, flag, var, nflg;

   skipSize = 5 + matrix.num_pars_inf_sites,
   lhsOffSet = 4;
   rhsOffSet = 1+matrix.num_taxa;

				      // Get size of LHS tree
   ctsl = currTask[0][0];
				      // Get size of RHS tree
   ctsr = currTask[0][1];
				      // Get highest internal node
   hin = currTask[0][2]-1;

      ntsl = ctsl;
      ntsr = ctsr + 2*skipSize;

      nextTask[0][0] = ntsl;
      nextTask[0][1] = ntsr;
      nextTask[0][2] = hin;
      nextTask[0][3] = taxa;
      nextTask[0][5] = currTask[0][5]+1;
      nextTask[0][lhsOffSet+skipSize+1] = 4;

      nextTask[1][0] = ntsr;
      nextTask[1][2] = 4;

                                      // Copy all elements of LHS tree 
				      // from parent
      memcpy(&nextTask[0][lhsOffSet+skipSize],                                 \
		      &currTask[0][lhsOffSet+skipSize],                        \
		      (ntsl-lhsOffSet-skipSize)*sizeof(int));

                                      // copy all elements from 0 to pos 
				      // from parent to RHS tree
      memcpy(&nextTask[1][1], &currTask[1][1], (pos-1)*sizeof(int));

                                      // 2 elems from "pos" get new values 
				      // (internal node, taxa) pair
      nextTask[1][pos] = hin;
      nextTask[1][pos+skipSize] = taxa;
      nextTask[1][pos+skipSize+4] = 0;

                                      // copy states of newly added leaf node
      memcpy(&nextTask[1][pos+skipSize+5], matrix.reord_sites_enc[taxa-1],     \
		                    matrix.num_pars_inf_sites*sizeof(int));

                                      // copy all remaining elements
      memcpy(&nextTask[1][pos+2*skipSize], &currTask[1][pos],                  \
		      (ntsr-pos-2*skipSize)*sizeof(int));

                                      // set parent child relationship for the 
				      // child tree
      for (sp=-1, ctr=ntsr-skipSize; ctr > 0; ctr = ctr-skipSize) {
       if (nextTask[1][ctr] > 0) {
         stack[++sp] = ctr;
       } else { 
         nextTask[1][stack[sp]+1] = ctr;
	 nextTask[1][ctr+2] = stack[sp];

	 sp = sp-1;
	 nextTask[1][stack[sp]+1] = ctr;
	 nextTask[1][ctr+3] = stack[sp];

	 stack[sp] = ctr;
       }
      }

                                      // compute cost of right half of the newly
				      // created tree using optimized logic
  si = pos;
  do {
        lci = nextTask[1][si+2];
        rci = nextTask[1][si+3];
        for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {

#ifdef MEMONL
           nextTask[1][si+5+site] = fitchArray[nextTask[1][lci+site+5]][nextTask[1][rci+5+site]][0];
           length = length + fitchArray[nextTask[1][lci+site+5]][nextTask[1][rci+5+site]][1];
#endif

#ifdef MEMVAR
           lc = nextTask[1][lci+site+5];
	   rc = nextTask[1][rci+5+site];
           nextTask[1][si+5+site] = fitchArray[lc][rc][0];
           length = length + fitchArray[lc][rc][1];
#endif

#ifdef CPU
           lc = nextTask[1][lci+site+5];
	   rc = nextTask[1][rci+5+site];
	   pa = lc | rc;
	   var = lc & rc;
	   flag = !var;
           nflg = ~flag;
           nflg = nflg + 1;
	   pa = pa & nflg;
	   pa = pa | var;
	   length = length + flag;
           nextTask[1][si+5+site] = pa;
#endif
        }

                                       /* New Check added */
        /*
        if (nextTask[1][lci+4]+nextTask[1][rci+4]+length > bestCost) {
           nextTask[0][8] = 1000053612;
           goto SKIPALLRIGHT;
        }
        */

        nextTask[1][si+4] = nextTask[1][lci+4]+nextTask[1][rci+4]+length;
        si = nextTask[1][si+1];
  } while (si!=4);

  for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {

#ifdef MEMONL
     nextTask[0][9+site] = fitchArray[nextTask[0][lhsOffSet+skipSize+5+site]][nextTask[1][6+site]][0];
     length = length + fitchArray[nextTask[0][lhsOffSet+skipSize+5+site]][nextTask[1][6+site]][1];
#endif

#ifdef MEMVAR
     lc = nextTask[0][lhsOffSet+skipSize+5+site];
     rc = nextTask[1][6+site];
     nextTask[0][9+site] = fitchArray[lc][rc][0];
     length = length + fitchArray[lc][rc][1];
#endif

#ifdef CPU
     lc = nextTask[0][lhsOffSet+skipSize+5+site];
     rc = nextTask[1][6+site];
     pa = lc | rc;
     var = lc & rc;
     flag = !var;
     nflg = ~flag;
     nflg = nflg + 1;
     pa = pa & nflg;
     pa = pa | var;
     length = length+flag;
     nextTask[0][9+site] = pa;
#endif
  }

  nextTask[0][8] = nextTask[0][lhsOffSet+skipSize+4]+nextTask[1][5]+length;

/*
SKIPALLRIGHT:
  lc = 0;
*/

}

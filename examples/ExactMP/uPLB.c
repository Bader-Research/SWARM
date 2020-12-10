/******************************************************************************
* uPLB.c                                                                      *
* Computes PLB  in unoptimized way. "u" prefix indicates "unoptimized"        *
* This was the initial version, used when taxa is <= 8                        *
******************************************************************************/

#include "header.h"

int purdomsLowerBnd(struct taskNode *tNode) {

  int site, added, length, tDoneIndx, tRemIndx, bitsInInt, remStatesORed, maskVal, maskPos; 
  struct taxaQ *taxaRemNode, *taxaDoneNode; 

  length = 0;
  bitsInInt = sizeof(int)*8;

  for (site = 0; site < matrix.num_pars_inf_sites; site++) { 
  	remStatesORed = 0;
	for (taxaRemNode=tNode->taxaRemHead; taxaRemNode != NULL; taxaRemNode=taxaRemNode->next) { 
		tRemIndx = taxaRemNode->taxa-1; 
  		remStatesORed = matrix.reord_sites_enc[tRemIndx][site] | remStatesORed;
	}
	maskVal = 1;
        for (maskPos=0; maskPos < bitsInInt-1; maskPos=maskPos+1) {
		added = 0;
                if (remStatesORed & maskVal) {
			for (taxaDoneNode=tNode->taxaDoneHead; taxaDoneNode != NULL;taxaDoneNode=taxaDoneNode->next) { 
				tDoneIndx = taxaDoneNode->taxa-1; 
	        		if (remStatesORed & matrix.reord_sites_enc[tDoneIndx][site]) {
					added = 1; 
             				break; 
	        		} 
			} 
             		if (!added) { 
	        		++length; 
             		} 
			remStatesORed = maskStateInPos(remStatesORed, maskPos);
               	}
               	maskVal = maskVal << 1;
        }
  }
  return (length);
}

int maskStateInPos(int stateVector, int bitPos) {
	int mask;

	mask = 0;
	mask = 1 << bitPos;
	mask = ~mask;
	stateVector = stateVector & mask;

	return (stateVector);
}

void showBits(int stateVector) {
	int bitsInInt, pos, mask=0, bit;

	bitsInInt = sizeof(int)*8;
	mask = 1 << (bitsInInt-2);
	for (pos=bitsInInt-2; pos >= 0; pos=pos-1) {
		bit = (stateVector & mask) ? 1: 0;
		printf("%d",bit);
		mask = mask >> 1;
	}
	printf("\n");
}

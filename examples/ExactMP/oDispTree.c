/******************************************************************************
* oDispTree.c                                                                 *
* Only for printing the tree, does no computation                             *
******************************************************************************/
#include "header.h"

#define LP 1000000
#define RP 2000000
#define COMA 3000000

void displayTree(int **cst, int lhsOffSet, int skipSize, THREADED) {
  int ctr=0, tmp=0;

  printf("\n\t\t\t\t\t---------------------------------------\n");
  fflush(stdout);
  printf("\t\t\t\t\t     TREE IN PREORDERED ARRAY FORM (T%d)\n", MYTHREAD);
  fflush(stdout);
  printf("\t\t\t\t\t---------------------------------------\n");
  fflush(stdout);
  printf("\n\t\tHIGHEST INTERNAL NODE: %d\t",cst[0][2]);
  fflush(stdout);
  printf("LAST TAXON ADDED: %d\t",cst[0][3]);
  fflush(stdout);
  printf("COST OF TASK: %d\n\n",cst[0][8]);
  fflush(stdout);

  printf("\tLCHILD\tNODE\tRCHILD\tPARENT\t");
  fflush(stdout);
#ifdef PRINT_SV
  printf("  COST UNTIL HERE\t{STATE VECTOR}\n");
  fflush(stdout);
#endif
#ifndef PRINT_SV
  printf("  COST UNTIL HERE\n");
  fflush(stdout);
#endif

  printf("ROOT NODE\n");
  fflush(stdout);
  printf("\t%3d\t root\t%3d\t root\t\t%3d\t\t",cst[0][lhsOffSet+skipSize],     \
		                                    cst[1][1], cst[0][8]);
  fflush(stdout);
#ifndef PRINT_SV
  printf("\n");
#endif

#ifdef PRINT_SV
  printf("{");
  fflush(stdout);
  for (ctr=9; ctr < 9+matrix.num_pars_inf_sites; ctr++) {
         printf("%2d ",cst[0][ctr]);
         fflush(stdout);
  }
  printf("}\n");
  fflush(stdout);
#endif

  printf("LEFT TREE\n");
  fflush(stdout);
  if (cst[0][lhsOffSet+skipSize] < 0) {
    printf("\t%3d\t%3d\t%3d\t root\t\t%3d\t\t",                                \
		   cst[0][cst[0][lhsOffSet+skipSize+2]],                       \
		   cst[0][lhsOffSet+skipSize],                                 \
		   cst[0][cst[0][lhsOffSet+skipSize+3]],                       \
		   cst[0][lhsOffSet+skipSize+4]);
    fflush(stdout);
  } else {
    printf("\t null\t%3d\t null\t root\t\t%3d\t\t", cst[0][lhsOffSet+skipSize],\
		   cst[0][lhsOffSet+skipSize+4]);
    fflush(stdout);
  }

#ifndef PRINT_SV
  printf("\n");
  fflush(stdout);
#endif

#ifdef PRINT_SV
  printf("{");
  fflush(stdout);
  for (ctr=lhsOffSet+skipSize+5; ctr < lhsOffSet+2*skipSize; ctr++) {
         printf("%2d ",cst[0][ctr]);
         fflush(stdout);
  }
  printf("}\n");
  fflush(stdout);
#endif

  for (ctr=lhsOffSet+2*skipSize; ctr < cst[0][0]; ctr=ctr+skipSize) {
     if (cst[0][ctr] > 0) {
           printf("\t null\t%3d\t null\t%3d\t\t%3d\t\t", cst[0][ctr],          \
			   cst[0][cst[0][ctr+1]], cst[0][ctr+4]);
           fflush(stdout);
     } else {
           printf("\t%3d\t%3d\t%3d\t%3d\t\t%3d\t\t",cst[0][cst[0][ctr+2]],     \
	   cst[0][ctr], cst[0][cst[0][ctr+3]], cst[0][cst[0][ctr+1]],          \
	   cst[0][ctr+4]);
           fflush(stdout);
     }
#ifndef PRINT_SV
  printf("\n");
  fflush(stdout);
#endif

#ifdef PRINT_SV
     printf("{");
     fflush(stdout);
     for (tmp=ctr+5; tmp < ctr+5+matrix.num_pars_inf_sites; tmp++) {
         printf("%2d ",cst[0][tmp]);
         fflush(stdout);
     }
     printf("}\n");
     fflush(stdout);
#endif
  }

  printf("RIGHT TREE\n");
  fflush(stdout);
  if (cst[1][1] < 0) {
    printf("\t%3d\t%3d\t%3d\t root\t\t%3d\t\t",cst[1][cst[1][3]],              \
		   cst[1][1], cst[1][cst[1][4]], cst[1][5]);
    fflush(stdout);
  } else {
    printf("\t null\t%3d\t null\t root\t\t%3d\t\t", cst[1][1], cst[1][5]);
    fflush(stdout);
  }
#ifndef PRINT_SV
  printf("\n");
  fflush(stdout);
#endif

#ifdef PRINT_SV
  printf("{");
  fflush(stdout);
  for (ctr=6; ctr < 6+matrix.num_pars_inf_sites; ctr++) {
         printf("%2d ",cst[1][ctr]);
         fflush(stdout);
  }
  printf("}\n");
  fflush(stdout);
#endif

  for (ctr=1+skipSize; ctr < cst[1][0]; ctr=ctr+skipSize) {
     if (cst[1][ctr] > 0) {
           printf("\t null\t%3d\t null\t%3d\t\t%3d\t\t", cst[1][ctr],          \
			   cst[1][cst[1][ctr+1]], cst[1][ctr+4]);
           fflush(stdout);
     } else {
           printf("\t%3d\t%3d\t%3d\t%3d\t\t%3d\t\t",cst[1][cst[1][ctr+2]],     \
	   cst[1][ctr], cst[1][cst[1][ctr+3]], cst[1][cst[1][ctr+1]],          \
	   cst[1][ctr+4]);
           fflush(stdout);
     }
#ifndef PRINT_SV
  printf("\n");
  fflush(stdout);
#endif

#ifdef PRINT_SV
     printf("{");
     fflush(stdout);
     for (tmp=ctr+5; tmp < ctr+5+matrix.num_pars_inf_sites; tmp++) {
         printf("%2d ",cst[1][tmp]);
         fflush(stdout);
     }
     printf("}\n");
     fflush(stdout);
#endif
  }
}

void displayTreeInNexusFormat(int **cst, int lhsOffSet, int skipSize, THREADED){
  int ctr=0, tmp=0;
  int **childParentArray, sizeOfArray;

  int currSize, loop, pos, *newArray, *currArray;

                                         /* First sort parent nodes out and 
                                            write their child nodes in an array 
                                            Array is indexed by abs(parent #),
                                            [0] is left & [1] is right child
                                         */
  sizeOfArray = abs(cst[0][2])+1;

  childParentArray = (int **)calloc(sizeOfArray, sizeof(int *));
  assert(childParentArray);
  for (ctr=0; ctr < sizeOfArray; ctr++) {
    childParentArray[ctr] = (int *)calloc(2, sizeof(int));
    assert(childParentArray[ctr]);
  }

  /* printf("ROOT NODE\n"); */
  childParentArray[0][0] = cst[0][lhsOffSet+skipSize];
  childParentArray[0][1] = cst[1][1];

  /* printf("LEFT TREE\n"); */
  if (cst[0][lhsOffSet+skipSize] < 0) {
    tmp = abs(cst[0][lhsOffSet+skipSize]);
    childParentArray[tmp][0] = cst[0][cst[0][lhsOffSet+skipSize+2]];
    childParentArray[tmp][1] = cst[0][cst[0][lhsOffSet+skipSize+3]];
  } 

  for (ctr=lhsOffSet+2*skipSize; ctr < cst[0][0]; ctr=ctr+skipSize) {
     if (cst[0][ctr] < 0) {
       tmp = abs(cst[0][ctr]);
       childParentArray[tmp][0] = cst[0][cst[0][ctr+2]];
       childParentArray[tmp][1] = cst[0][cst[0][ctr+3]];
     }
  }

  /* printf("RIGHT TREE\n"); */
  if (cst[1][1] < 0) {
       tmp = abs(cst[1][1]);
       childParentArray[tmp][0] = cst[1][cst[1][3]];
       childParentArray[tmp][1] = cst[1][cst[1][4]];
  }

  for (ctr=1+skipSize; ctr < cst[1][0]; ctr=ctr+skipSize) {
     if (cst[1][ctr] < 0) {
       tmp = abs(cst[1][ctr]);
       childParentArray[tmp][0] = cst[1][cst[1][ctr+2]];
       childParentArray[tmp][1] = cst[1][cst[1][ctr+3]];
     }
  }

  /* 
  for (ctr=0; ctr < sizeOfArray; ctr++) {
    printf("%d\t%d\t%d\n",childParentArray[ctr][0],-ctr, childParentArray[ctr][1]);
  } 
  */
 
                                              /* Now run the algorithm to expand
                                                 all internal nodes, level-by-level
                                              */
  loop = 1;
  currSize = 5;
  currArray = (int *)calloc(currSize, sizeof(int));
  currArray[0] = LP;
  currArray[1] = childParentArray[0][0];
  currArray[2] = COMA;
  currArray[3] = childParentArray[0][1];
  currArray[4] = RP;
  while (loop) {
    pos = 0;
    loop = 0;
    while (pos < currSize) {
      if (currArray[pos] < 0) {
        currSize += 4;
        newArray = (int *)calloc(currSize, sizeof(int));
        memcpy(&newArray[0], &currArray[0], pos*sizeof(int));
        newArray[pos] = LP;
        newArray[pos+1] = childParentArray[abs(currArray[pos])][0];
        newArray[pos+2] = COMA;
        newArray[pos+3] = childParentArray[abs(currArray[pos])][1];
        newArray[pos+4] = RP;
        if (newArray[pos+1] < 0) loop = 1;
        if (newArray[pos+3] < 0) loop = 1;
        pos += 5;
        memcpy(&newArray[pos], &currArray[pos-4], (currSize-pos)*sizeof(int));
        free(currArray);
        currArray = newArray;
      }
      ++pos;
    } /* pos < currSize */
  } /* loop */ 

  for (tmp=0; tmp < currSize; tmp++) {
    if (currArray[tmp] == LP) {
       printf("(");
       fflush(stdout);
    } else if (currArray[tmp] == RP) {
       printf(")");
       fflush(stdout);
    } else if (currArray[tmp] == COMA) {
       printf(",");
       fflush(stdout);
    } else {
       printf("%s",matrix.taxons[currArray[tmp]-1]);
       fflush(stdout);
    }
  }
  free(currArray);
  for (ctr=0; ctr < sizeOfArray; ctr++) {
    free(childParentArray[ctr]);
  }
  free(childParentArray);
}

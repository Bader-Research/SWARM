#include "header.h"

void nj(THREADED) {                      /* Output file */
                                         /* Getting the distance matrix from the sequence matrix */
   int ctr, i, j, k;

   NJ.taxon_table = (double **)calloc(matrix.num_taxa, sizeof(double *));
   assert(NJ.taxon_table);
   for (ctr=0; ctr < matrix.num_taxa; ctr++) {
	   NJ.taxon_table[ctr] = (double *)calloc(matrix.num_taxa, sizeof(double));
	   assert(NJ.taxon_table[ctr]);
   }

                                         /* THIS SECTION CALCULATES SIMILARITY MATRIX */
   for(i=0; i < matrix.num_taxa; i++) {
      for(j=i+1; j < matrix.num_taxa; j++) {
	 NJ.taxon_table[i][j]=0.0;
         for(k=0; k < matrix.num_sites; k++) {                                     
		                        /* Exclude sites with missing entries */
            if (matrix.taxon_table[i][k]!='?' && matrix.taxon_table[j][k]!='?' && matrix.taxon_table[j][k]!=' ' && matrix.taxon_table[j][k]!=' ' && matrix.taxon_table[i][k]!='-' && matrix.taxon_table[j][k]!='-') {
		                        /* This being a symmetrical matrix, calculate the upper triangle alone */
               if ((matrix.taxon_table[i][k]!=matrix.taxon_table[j][k])) {
                   NJ.taxon_table[i][j] = NJ.taxon_table[i][j]+1.0;
	       }
	    }
	 }
      }
    }

                                        /* Fill the lower triangular matrix */
    for(i=0;i<matrix.num_taxa;i++) {
       for(j=i+1;j<matrix.num_taxa;j++) {
          NJ.taxon_table[j][i] = NJ.taxon_table[i][j];
       }
    }

                                         /* THIS SECTION PRINTS OUT SIMILARITY MATRIX */
#ifdef NJ_VERBOSE
    printf("\nSIMILARITY MATRIX FOR NJ ALGORITHM...\n");
    for (i=0; i< matrix.num_taxa; i++) {
        for (j=0; j < matrix.num_taxa; j++) {
            printf("%f ",NJ.taxon_table[i][j]);
	}
        printf("\n");
    }
#endif

                                         /* CALL NJ TREE WITH SK MODIFICATION */
   neighborJoin(TH);

   for (ctr=0; ctr < matrix.num_taxa; ctr++) {
	   free(NJ.taxon_table[ctr]);
   }
   free(NJ.taxon_table);
}


/********************************************* Neighbor-Joining tree ****************************************************/

void neighborJoin (THREADED) {

    double *skParam;                          /* S&K parameter */
    int *addToNJTree;                         /* Array holding taxa to be joined to NJ tree */
    char **branch;                            /* #i entry in this list hold the subtree of which i is a part */
    int maxConstSize;                         /* Max size of tree list */
    int ***njMatrix;
    char *tmpbranch;                          /* temp variable and final tree (listing of taxas) respectively */

    int i, j, k;                              /* Scratch variables */
    int flag;                                 /* Flags */
    double NJtreescore;                       /* Holds the NJ score */
    int N;
    double Dij;                               /* Studier and Keppler(S&K) method parameters */

    double Sij, minSij;                       /* S&K Paramrters */
    int mini, minj;                           /* #taxas that are closest negihbours currently */

    NJUniqStates = 1;                         /* Initializing to one unique state */
    MPScoreFromNJ = 0;                        /* MP Score= length initialised to zero */

    NJ.uniq_chars = (char *) malloc (MAXSTATES*sizeof(char));
    assert(NJ.uniq_chars);

                                              /* The first unique character is first character in sequence matrix */
    NJ.uniq_chars[0]= (char)(matrix.taxon_table[0][0]);   

    skParam = (double *)calloc(matrix.num_taxa, sizeof(double));
    assert(skParam);

    addToNJTree = (int *)calloc(matrix.num_taxa, sizeof(int));
    assert(addToNJTree);

    branch = (char **)calloc(matrix.num_taxa, sizeof(char *));
    assert(branch);

    maxConstSize = ((int)ceil(log10((double)matrix.num_taxa))+4)*matrix.num_taxa + 1;

	                                                   /* Allocate memory to branch, & initialize it with taxa # */
    for (i = 0; i < matrix.num_taxa ; i++) {
        branch[i] = (char *)calloc(maxConstSize+10, sizeof(char));
	assert(branch[i]);
        sprintf (branch[i],"%d",i+1);                                      
    }

    tmpbranch = (char *)calloc(maxConstSize, sizeof(char));
    assert(tmpbranch);

    NJtreescore = 0;                                       /* Initializing NJ score to 0 */

    for ( i = 0; i < matrix.num_taxa; i++ ) {              /* Initially, all taxa are to be considered in NJ Algorithm */
        addToNJTree[i] = TRUE;
    }

                                                           /* Matrix to hold the states for all taxa */
    njMatrix = ( int ***)calloc(matrix.num_taxa, sizeof(int **));           
    assert(njMatrix);
    for (i=0; i < matrix.num_taxa; i++) {
        njMatrix[i] = (int **)calloc(matrix.num_sites, sizeof(int *));
        assert(njMatrix[i]);
        for (j=0; j<matrix.num_sites; j++) {
	                                                   /* Initially given an arbitary amount of memory as */ 
		                                           /* the number of unique states is not known */
            njMatrix[i][j] = (int *)calloc(MAXSTATES,sizeof(int));            
	    assert(njMatrix[i][j]);
            for (k=0; k<MAXSTATES; k++) {
                njMatrix[i][j][k] = 0;                     /* Initialising all states to 000..00 */
	    }
	}
    }

                                                           /* THIS SECTION ASSIGNS STATE VALUES */
   NJUniqStates = 1;                                       /* Initialising to single state */
   for (i=0; i<matrix.num_taxa; i++) {
                                                           /* Find # of uniques states, & characters representing them */
	                                                   /* in entire sequence matrix */
       for (j=0; j<matrix.num_sites; j++) {
           for (k=0; k<NJUniqStates; k++) {
               if (matrix.taxon_table[i][j]==NJ.uniq_chars[k]) {
		                                           /* Set 1 for the corresponding state if exists */
                   njMatrix[i][j][k] = 1;
                   break;
	       }
	   }
           if (k==NJUniqStates) {                          /* New state found*/
              njMatrix[i][j][NJUniqStates] = 1;            /* Set 1 for state */
                                                           /* Add new state to the list of unique states */
              NJ.uniq_chars[NJUniqStates] = (char)matrix.taxon_table[i][j];
              NJUniqStates = NJUniqStates + 1;             /* Increment the number of unique states by one */
           }
       }
   }

#ifdef NJ_VERBOSE
   printf("\nNumber of unique states = %d : ",NJUniqStates);
   for (i=0; i<NJUniqStates; i++) {
       printf("%c",NJ.uniq_chars[i]);
   }
#endif

                                                           /* Calculating parameters for NJ and assigning neighbors */
	                                                   /* Recursively find nearest neighbors, & join until only */ 
                                                           /* two effective nodes are left to join */
   for (N=matrix.num_taxa; N > 2; N--) {
        for (i = 0; i < matrix.num_taxa; i++ ) {
            if (addToNJTree[i]) {                          /* Only for taxa that are yet to be joined to tree */
                skParam[i] = 0.0;                          /* Calculate S&K Parameter */
                for ( j = 0; j < matrix.num_taxa; j++ ) {
                    if ((addToNJTree[j]) && (i != j)) {
                        if ( i < j )
                            skParam[i] += NJ.taxon_table[i][j];
                        else
                            skParam[i] += NJ.taxon_table[j][i];
                    }
                    skParam[i] = skParam[i]/((double)N - 2.0);
		}
            }
        }

                                                           /* Find nearest neighbors of all available taxa */
        minSij = DBL_INF;
        mini = -1;
        minj = -1;
        for (i = 0; i < matrix.num_taxa; i++) {
 	    if (addToNJTree[i]) {
                for (j = i + 1; j < matrix.num_taxa; j++ ) {
                    if (addToNJTree[j]) {
                        Sij = NJ.taxon_table[i][j] - ( skParam[i] + skParam[j] );
                        if ( Sij < minSij ) {
                            minSij = Sij;
                            mini = i;
                            minj = j;
                        }
                    }
                }
            }
        }

                                                           /* Calculating length of subtree formed until now */
	                                                   /* (using Fitch's Parsimony) */
      for (j=0; j< matrix.num_sites; j++) {
          flag=0;
          for (k=0; k<NJUniqStates; k++) {
              flag = flag + (njMatrix[mini][j][k])*(njMatrix[minj][j][k]);
	  }
          if (flag==0) {
	                                                   /* If no states is common, union of states of daughter */
		                                           /* taxa is assigned to parent and length incremented */
             MPScoreFromNJ = MPScoreFromNJ + 1;    
             for (k=0; k<NJUniqStates; k++) {
                 njMatrix[mini][j][k] = njMatrix[mini][j][k]+njMatrix[minj][j][k];
	     }
	  }
          else {
	                                                   /* Otherwise intersection of states */
              njMatrix[mini][j][k] = njMatrix[mini][j][k]*njMatrix[minj][j][k];       
	  }
      }
      strcpy (tmpbranch, branch[mini]);
      strcpy (tmpbranch, ",");
      strcpy (tmpbranch, branch[mini]);
                                                           /* Copy newly added taxa to tree list */
      sprintf (branch[mini], "(%s)", tmpbranch);   

                                                           /* Re-estimating distances after joining taxa to subtree */
#ifdef NJ_VERBOSE
      printf("\nRe-estimating distances after joining taxa to subtree...\n");
#endif
      addToNJTree[minj] = FALSE;
      Dij = NJ.taxon_table[mini][minj];
      for (j = 0; j < mini; j++ ) {
            if (addToNJTree[j]) {
                if (minj < j)
                    NJ.taxon_table[j][mini] = 0.5 * ( NJ.taxon_table[j][mini] + NJ.taxon_table[minj][j] - Dij );
                else
                    NJ.taxon_table[j][mini] = 0.5 * ( NJ.taxon_table[j][mini] + NJ.taxon_table[j][minj] - Dij );
            }
      }
      for (j = mini+1; j < matrix.num_taxa; j++) {
            if (addToNJTree[j]) {
                if (minj < j)
                    NJ.taxon_table[mini][j] = 0.5 * ( NJ.taxon_table[mini][j] + NJ.taxon_table[minj][j] - Dij );
                else
                    NJ.taxon_table[mini][j] = 0.5 * ( NJ.taxon_table[mini][j] + NJ.taxon_table[j][minj] - Dij );
            }
      }
      NJtreescore += Dij;
#ifdef NJ_VERBOSE
      printf("NJ score = %f\n",NJtreescore);
#endif
   }
    
#ifdef NJ_VERBOSE
   printf("NJ score = %f\n",NJtreescore);
   printf("Adding last taxa..\n");
#endif

   mini = minj = -1;
   for ( i = 0; i < matrix.num_taxa; i++ ) {
       if (addToNJTree[i]) {
          if ( mini < 0 )                                  /* Adding last two taxa */
                mini = i;
          else
                minj = i;
       }
   }

   Dij = NJ.taxon_table[mini][minj];
   NJtreescore+= Dij;
                                                           /* adding length due to last taxa alone */
   for (j=0; j< matrix.num_sites; j++) {
       flag=0;
       for (k=0; k<NJUniqStates; k++) {
           flag = flag + (njMatrix[mini][j][k])*(njMatrix[minj][j][k]);
       }
       if (flag==0) {
          MPScoreFromNJ = MPScoreFromNJ +1;               /* If no states is common, union of states of daughter */
	                                                  /* taxa is assigned to parent and length incremented */
          for (k=0; k<NJUniqStates; k++) {
              njMatrix[mini][j][k] = njMatrix[mini][j][k]+njMatrix[minj][j][k];
	  }
       }
       else {
	                                                  /* Otherwise intersection of states */
         njMatrix[mini][j][k] = njMatrix[mini][j][k]*njMatrix[minj][j][k];       
       }
   }

   printf("NJ Score: %.2f, and corresponding MP Score = %d\n",NJtreescore, MPScoreFromNJ);

   free(skParam);
   free(addToNJTree);
   for (i = 0; i < matrix.num_taxa ; i++) {
	free(branch[i]);
   }
   free(branch);
   free(tmpbranch);

   for (i=0; i < matrix.num_taxa; i++) {
        for (j=0; j<matrix.num_sites; j++) {
	    free(njMatrix[i][j]);
	}
   }
   for (i=0; i < matrix.num_taxa; i++) {
        free(njMatrix[i]);
   }
   free(njMatrix);
   free(NJ.uniq_chars);

}

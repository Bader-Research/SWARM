/******************************************************************************
* preProcess.c                                                                *
* Does initial preprocessing: site encoding, reordering, removing PNI sites   *
******************************************************************************/
#include "header.h"

void preProcess(char *fname, THREADED) {
   int ptr;
   FILE *fp;
   char **seq_matrix;                  /* Holds input sequence matrix */
   char **taxa_name;                   /* Holds names of taxa */
   int **enc_st_tbl;                   /* Holds state encoding of the matrix */
   char **unq_st_tbl;                  /* Holds unique states in the site */
   int dup_st;                         /* Flag: = 1 if state prev occured */
   char *site_st_tbl;                  /* Holds state of current site */
   int *num_states;                    /* Indx 0 holds # states in site 0 */
   int row_v1, col_v1, ctr2;           /* Scratch pad */
   char ch;

   dna = 1;

   fp = fopen(fname,"r");   
   if (fp==NULL) {
     perror(fname);
     exit(-1);
   }
                                       /* Get # of taxa & sites */
   ptr=fscanf(fp,"%d%d\n",&matrix.num_taxa,&matrix.num_sites);

   seq_matrix = (char **)calloc(matrix.num_taxa, sizeof(char *));
   assert(seq_matrix);
   for (col_v1=0; col_v1 < matrix.num_taxa; col_v1++) {
      seq_matrix[col_v1] = (char *)calloc(matrix.num_sites, sizeof(char));
      assert(seq_matrix[col_v1]);
   }
   taxa_name = (char **)calloc(matrix.num_taxa, sizeof(char *));
   assert(taxa_name);
   for (col_v1=0; col_v1 < matrix.num_taxa; col_v1++) {
      taxa_name[col_v1] = (char *)calloc(TAXA_NAME_LENGTH, sizeof(char));
      assert(taxa_name[col_v1]);
   }

                                       /* Read input */
   for (row_v1=0; row_v1 < matrix.num_taxa; row_v1++) {
      fscanf(fp,"%s",taxa_name[row_v1]);       
      for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
         ctr2 = fscanf(fp,"%c",&ch);
         while (ctr2 && (ch==' ')) {         /* Ignore all blanks */
            ctr2 = fscanf(fp,"%c",&ch);
         }
         assert(ctr2 > 0);
         seq_matrix[row_v1][col_v1] = ch;
	 if (ch != 'A' && ch != 'C' && ch != 'T' && ch != 'G') {
	    /*
            printf("%c ",ch);
	    fflush(stdout);
            */
            dna=0;
	 } 
      }
   }

   if (dna) {
    if (verbose) {
    printf("\nInput matrix contains DNA sequence\n");
    fflush(stdout);
    }
   } else {
    if (verbose) {
    printf("Input matrix is non-DNA sequence\n");
    fflush(stdout);
    }
   }

   site_st_tbl = (char *)calloc(matrix.num_taxa,sizeof(char));
   num_states  = (int *)calloc(matrix.num_sites,sizeof(int));

                                       /* memory for state encodings table */
   enc_st_tbl = (int **)calloc(matrix.num_taxa, sizeof(int *));
   assert(enc_st_tbl);
   for (col_v1=0; col_v1 < matrix.num_taxa; col_v1++) {
      enc_st_tbl[col_v1] = (int *)calloc(matrix.num_sites, sizeof(int));
      assert(enc_st_tbl[col_v1]);
   }
                                       /* memory for unique states @ each site*/
   unq_st_tbl = (char **)calloc(matrix.num_sites, sizeof(char *));
   assert(unq_st_tbl);

   for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
                                       /* Encode 1st st to 1 in state table */
      enc_st_tbl[0][col_v1] = 1;
                                       /* Init # of states per site to 1 */
      num_states[col_v1] = 1;
   }

                                       /* For all sites */
   for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
	                               /* 1st st of st tbl = 1st st in seq mat*/
       site_st_tbl[0] = seq_matrix[0][col_v1];
                                       /* For all taxon */
       for (row_v1=1; row_v1 < matrix.num_taxa; row_v1++) {
	                               /* If state occured, mark as duplicate*/
	  dup_st=0;
	  for (ctr2=0; ctr2 < num_states[col_v1]; ctr2++) {
	      if (seq_matrix[row_v1][col_v1]==site_st_tbl[ctr2]) {
		 dup_st=1;
		 break;
	      }
	  }
	                               /* If state didnt occur before */
	  if (!dup_st) {
	                               /* Append to st list for the site */
	     site_st_tbl[num_states[col_v1]] = seq_matrix[row_v1][col_v1];
	                               /* Incr # of states for this site */
	     num_states[col_v1] += 1;
	  }
       }
	                               /* memory for unique states @ the site */
       unq_st_tbl[col_v1] = (char *)calloc(num_states[col_v1],sizeof(char));
       assert(unq_st_tbl[col_v1]);

       for (ctr2=0; ctr2 < num_states[col_v1]; ctr2++) {
          unq_st_tbl[col_v1][ctr2] = site_st_tbl[ctr2];
       }

                                       /* For all taxa in the matrix */
       for (row_v1=0; row_v1 < matrix.num_taxa; row_v1++) {
                                       /* For all states in the current site */
           for (ctr2=0; ctr2 < num_states[col_v1]; ctr2++) {
                                       /* State encoding = 1 left shfted by 
                                          its rank in state list of the site */
	       if (seq_matrix[row_v1][col_v1]==site_st_tbl[ctr2]) {
		  enc_st_tbl[row_v1][col_v1] = 1 << ctr2;
		  break;
	       }
	   }
       }
   }

#ifdef PREPROCESS_VERBOSE
   for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
      printf("Site %2d has %2d states: ",col_v1,num_states[col_v1]);
      for (row_v1=0; row_v1 < num_states[col_v1]; row_v1++) {
          printf("%c ",unq_st_tbl[col_v1][row_v1]); 
      }
      printf("\n");
   }

   printf("\n");
   for (row_v1=0; row_v1 < matrix.num_taxa; row_v1++) {
      printf("%s\t",taxa_name[row_v1]);
      for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
         printf("%4c",seq_matrix[row_v1][col_v1]);
      }
      printf("\n");
   }

   printf("\n\t");
   for (row_v1=0; row_v1 < matrix.num_taxa; row_v1++) {
      for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
         printf("%3d ",enc_st_tbl[row_v1][col_v1]);
      }
      printf("\n\t");
   }
#endif

   matrix.taxon_table = seq_matrix;
   matrix.taxons = taxa_name;
   matrix.state_encoding = enc_st_tbl;
   matrix.num_st_at_site = num_states;
   matrix.unique_states = unq_st_tbl;

   free(site_st_tbl);

   fclose(fp);

}

void reorderSites(THREADED) {
   int row_v1, col_v1, ctr2;
   int **unq_st_reps;
   int **rol_table;
   int rank=0;

   struct reorder_list {
	   int site;
	   int rank;
	   struct reorder_list *left;
	   struct reorder_list *right;
   };
   
   struct reorder_list *rol_head = NULL;
   struct reorder_list *rol_tail = NULL;
   struct reorder_list *curr_rol_node = NULL;
   struct reorder_list *temp_rol_node = NULL;

   struct pNonInfSites {
	   int site;
	   struct pNonInfSites *next;
   };

   struct pNonInfSites *tempPniNode = NULL;
   struct pNonInfSites *pNiHead = NULL;
   struct pNonInfSites *pNiTail  = NULL;

   int num_pi_sites=0;                   /* # of parsimony informative sites */

   pNiCost=0;
                                         /* memory to hold reps of states */
   unq_st_reps = (int **)calloc(matrix.num_sites,sizeof(int *));
                                         /* For all sites */
   for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
                                         /* memory for unique states @ site */
      unq_st_reps[col_v1] =                                                 \
                   (int *)calloc(matrix.num_st_at_site[col_v1],sizeof(int ));
                                         /* Compare each unique state @ site 
                                            with all taxon in the site */
      for (ctr2=0; ctr2 < matrix.num_st_at_site[col_v1]; ctr2++) {
         for (row_v1=0; row_v1 < matrix.num_taxa; row_v1++) {
                                         /* if repretition, incr count */
	    if (matrix.taxon_table[row_v1][col_v1] ==                        \
                                          matrix.unique_states[col_v1][ctr2]) {
	        unq_st_reps[col_v1][ctr2] += 1;
	    }
         }
      }
   }
                                             /* Store array globally */
   matrix.unique_states_reps = unq_st_reps;
                                             /* algorithm to reorder sites */
                                             /* For all sites */
   for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
      rank = 0;
                                             /* proceed iff site has > 1 st */
      if (matrix.num_st_at_site[col_v1] > 1) {
                                             /* for all states in this site */
	 for (ctr2=0; ctr2 < matrix.num_st_at_site[col_v1]; ctr2++) {
		                             /* if state has > 1 reps note it */
	    if (unq_st_reps[col_v1][ctr2] > 1) {
	       rank = rank+1;
	    }
         }
	                                     /* if # states repeating multiple 
                                                times is > 1, this is PI site */
	                                     /* enqueue it in descending order 
                                                in the reorder list */
	 if (rank > 1) {
            num_pi_sites += 1;               /* count # of PI sites */
		                             /* if start of the list */ 
            if (rol_head==NULL) {
		                             /* memory for head node */ 
               rol_head =                                                     \
                    (struct reorder_list *)malloc(sizeof(struct reorder_list));
	       rol_head->right = NULL;
	       rol_head->left = NULL;
	       rol_tail = rol_head;
	       rol_head->site = col_v1;
	       rol_head->rank = rank;
	    }
	    else {
	       if (rank >= rol_head->rank) {
                  temp_rol_node =                                            \
                    (struct reorder_list *)malloc(sizeof(struct reorder_list));
		  temp_rol_node->left = NULL;
		  temp_rol_node->right = rol_head;
		  rol_head->left = temp_rol_node;
		  temp_rol_node->site = col_v1;
		  temp_rol_node->rank = rank;
		  rol_head = temp_rol_node;
	       }
	       else if (rank <= rol_tail->rank) {
                  temp_rol_node =                                            \
                    (struct reorder_list *)malloc(sizeof(struct reorder_list));
		  temp_rol_node->left = rol_tail;
		  temp_rol_node->right = NULL;
		  temp_rol_node->site = col_v1;
		  temp_rol_node->rank = rank;
		  rol_tail->right = temp_rol_node;
		  rol_tail = temp_rol_node;
	       }
	                                     /* Intermediate cond: get interval
                                                 where it fits in the list */
	       else {
                  for (curr_rol_node=rol_head->right; curr_rol_node!=NULL;    \
                                          curr_rol_node=curr_rol_node->right) {
			                     /* if rank falls betn a, b 
                                                st. a >= rank >= b , & 
                                                (a,rank,b) are consecutive */
                     
                     if ((curr_rol_node->left->rank >= rank) &&               \
                                               (rank >= curr_rol_node->rank)) {
                        temp_rol_node =                                       \
                    (struct reorder_list *)malloc(sizeof(struct reorder_list));

                        temp_rol_node->right = curr_rol_node;
			temp_rol_node->left = curr_rol_node->left;

			curr_rol_node->left->right = temp_rol_node;
			curr_rol_node->left = temp_rol_node;

		        temp_rol_node->site = col_v1;
		        temp_rol_node->rank = rank;

			break;
                     }
                  }
               }
	    }

#ifdef REORDER_VERBOSE
         printf("\n# parsimony informative sites until now: %d",num_pi_sites);
         printf(" Current list {site,rank} given below\n");
         for (curr_rol_node=rol_head; curr_rol_node!=NULL;                   \
                                       curr_rol_node=curr_rol_node->right) {
              printf("{%d,%d} ",curr_rol_node->site,curr_rol_node->rank);
         }
	 printf("\n");
#endif
         }
	                                  /* if PI site */
	                                  /* add constant length to cost = 
                                             num of states in site - 1 */
	 else {
            if (pNiHead ==NULL) {
		                          /* create memory for the head node */ 
               pNiHead  =                                                     \
                   (struct pNonInfSites *)malloc(sizeof(struct pNonInfSites )); 
	       pNiHead->next = NULL;
	       pNiTail  = pNiHead;
	       pNiHead->site = col_v1;
	    }
	    else {
                  tempPniNode =                                               \
                    (struct pNonInfSites *)malloc(sizeof(struct pNonInfSites));
		  tempPniNode->next = NULL;
		  pNiTail->next = tempPniNode;
		  tempPniNode->site = col_v1;
		  pNiTail = tempPniNode;
            }
	 }
      }
   }

#ifdef REORDER_VERBOSE
   printf("\n");
   for (col_v1=0; col_v1 < matrix.num_sites; col_v1++) {
      for (ctr2=0; ctr2 < matrix.num_st_at_site[col_v1]; ctr2++) {
          printf("Site: %d, State: %c, Reps: %d\n",col_v1,                    \
                matrix.unique_states[col_v1][ctr2], unq_st_reps[col_v1][ctr2]);
      }
      printf("\n");
   }
#endif
                                                                     
                                         /* Update global variable with number 
                                            of parsimony informative sites */
   matrix.num_pars_inf_sites = num_pi_sites;
	                                 /* memory for reordered sites table */
   rol_table = (int **)calloc(matrix.num_taxa,sizeof(int *));
   for (ctr2=0; ctr2 < matrix.num_taxa; ctr2++) {
     rol_table[ctr2] = (int *)calloc(num_pi_sites,sizeof(int ));
   }

	                                  /* Retain parsimony informative 
                                             sites alone with state enc */
   temp_rol_node = rol_head;
   for (col_v1=0; col_v1 < num_pi_sites; col_v1++) {
      ctr2 = temp_rol_node->site;
      for (row_v1=0; row_v1 < matrix.num_taxa; row_v1++) {
         rol_table[row_v1][col_v1] = matrix.state_encoding[row_v1][ctr2];
      }
      temp_rol_node = temp_rol_node->right;
   }

                                        /* Update global variable with 
                                           pars inf sites and their encoding */
   matrix.reord_sites_enc = rol_table;

#ifdef REORDER_VERBOSE
   printf("\nReordered and encoded table <encoded ");
   printf("state(original site number)>\n");
   printf("-------------");
   for (col_v1=0; col_v1 < num_pi_sites*8; col_v1++) {
      printf("-");
   }
   printf("\n");
   printf("   SITE:\t");
   temp_rol_node = rol_head;
   for (col_v1=0; col_v1 < num_pi_sites; col_v1++) {
      ctr2 = temp_rol_node->site;
      printf("%2d\t",ctr2);
      temp_rol_node = temp_rol_node->right;
   }
   printf("\n");
   printf("-------------");
   for (col_v1=0; col_v1 < num_pi_sites*8; col_v1++) {
      printf("-");
   }
   printf("\n");
   for (row_v1=0; row_v1 < matrix.num_taxa; row_v1++) {
      printf("Taxa %2d:\t",row_v1+1);
      for (col_v1=0; col_v1 < num_pi_sites; col_v1++) {
         printf("%2d\t",matrix.reord_sites_enc[row_v1][col_v1]);
      }
      printf("\n");
   }

   printf("\nParsimony Non Informative Sites with one or more states: ");
   printf("{Site (Number of states in that site)}:\n");
   pNiCost = 0;
   for (tempPniNode=pNiHead; tempPniNode!=NULL; tempPniNode=tempPniNode->next){
      printf("{%d (%d)} ",tempPniNode->site,                                 \
                                 matrix.num_st_at_site[tempPniNode->site]);
      pNiCost = pNiCost + matrix.num_st_at_site[tempPniNode->site]-1;
   }
   printf("\n");
   printf("\nCost due to parsimony non informative sites: %d\n",pNiCost);
#endif

	                                  /* Free reorder list */
   for (curr_rol_node=rol_head; curr_rol_node!=NULL;                         \
                                        curr_rol_node=curr_rol_node->right) {
      free(curr_rol_node);
   }

   pNiCost = 0;
   numPNISites = 0;
   for (tempPniNode=pNiHead; tempPniNode!=NULL; tempPniNode=tempPniNode->next){
      pNiCost = pNiCost + matrix.num_st_at_site[tempPniNode->site]-1;
      numPNISites = numPNISites + 1;
   }

   for (tempPniNode=pNiHead; tempPniNode!=NULL; tempPniNode=tempPniNode->next){
      free(tempPniNode);
   }

   numConstSites = matrix.num_sites - matrix.num_pars_inf_sites - numPNISites;

}

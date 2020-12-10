#ifndef _HEADER_H
#define _HEADER_H

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#include <pthread.h>
//#include <sys/time.h>
#include "swarm.h"

/*#define MEMVAR*/              /* Select for DNA */
#define CPU                     /* Select for non-DNA */

//#define PRINT_SV

#define keepTrees 100
double *loadDistribution;
pthread_mutex_t lock;

int dna,
    nextID,
    pNiCost,
    *taskID,
    ***solQ,
    verbose,
    *request,
    solQCtr,
    bestCost,
    *tbrBest,
    *taxaQueue,
    randOptLev,
    randOptFlag,
    numPNISites,
    NJUniqStates,                        
    numConstSites,
    MPScoreFromNJ,
    NUM_TBR_TREES,
    numBestSolutions;

int **initSol;

struct bestCostStack {
  struct subTask *bCNode;
  struct bestCostStack *down;
};

struct bestCostStack *bestCostStackHead;

struct integerStack {
   int id;
   struct integerStack *down;
};

struct integerStack *stackTop;

struct {
   int  num_taxa;                     /* Total number of input taxa (rows) */
   int  num_sites;                    /* Total number of input sites (cols) */
   char **taxons;                     /* Indx 0 => name of 0th taxon */
   char **taxon_table;                /* Indx r,c => rth taxa's cth character */
   int  **state_encoding;             /* Indx r,c => encoded taxon_table[r][c]*/
   int  **reord_sites_enc;            /* Reordered contains PI sites only */
   int  num_pars_inf_sites;           /* Total number of PI sites */
   int  *num_st_at_site;              /* Indx 0 => # unique states of site 0 */
   char **unique_states;              /* Holds unique states in a site */
   int  **unique_states_reps;         /* Indx 0 => # reps of indx 0 of */
                                      /* unique_states array */
} matrix;

struct subTaskArrayNode {             // (6+matrix.num_pars_inf_sites) integers
  int id;                             // 1 integer
  int parentIndx;                     // 1 integer
  int lChildIndx;                     // 1 integer
  int rChildIndx;                     // 1 integer
  int costAtThisNode;                 // 1 integer
  int costUntilThisNode;              // 1 integer
  int *stateVector;                   // matrix.num_pars_inf_sites integers
};

struct subTask {                      // (sizeOfSubTaskArray+1) X 
                                      //      (6+matrix.num_pars_inf_sites) + 4
  int sizeOfSubTaskArray;             // 1 integer
  int lastTaxaAdded;                  // 1 integer
  int highestInternalNode;            // 1 integer
  int costOfSubTask;                  // 1 integer
  int rhtSubTreeRootIndx;             // 0 integer NO LONGER PRESENT IN 1D ARR
  struct subTaskArrayNode root;       // 6 + matrix.num_pars_inf_sites
  struct subTaskArrayNode *node;      // sizeOfSubTaskArray X 
                                      //      (6 + matrix.num_pars_inf_sites)
};

struct taskNode {
  struct taxaQ *taxaRemHead;
  struct taxaQ *taxaRemTail;
  struct taxaQ *taxaDoneHead;
  struct taxaQ *taxaDoneTail;
  struct subTask *partialTask;
};

struct {
   char *uniq_chars;                  /* Char representing unique states */
   double **taxon_table;              /* Indx r,c => dist betn ith, jth taxa */
} NJ;

struct taskStack {
  struct taskNode *tNode;
  struct taskStack *down;
};

struct taxaQ {
  int taxa;
  struct taxaQ *next;
};

struct taskStack *taskSTop;

struct taxaQ *taxaQHead;
struct taxaQ *taxaQTail;

struct taxaQ *tAddOrdQHead;
struct taxaQ *tAddOrdQTail;

struct maxMiniTaxaPosWLength {
   int pos;
   int taxa;
   int length;
};

void assignJobs(THREADED);
void nj(THREADED);
void bnb(THREADED);
void eckDay(THREADED);
void TBRmaxMin(THREADED);
void bnb(THREADED);
void clearMemory(void);
void neighborJoin(THREADED);
void setPlbValues(THREADED);
void freeRandTaxaQ(THREADED);
void clearWorkSpace(THREADED);
void TBRSmartGreedy(THREADED);
void searchOptimize(THREADED);
void getKthLevel(int , THREADED);
int numNodesInRow(int , THREADED);
void allocMem4WorkSpace(THREADED);
void TBReckDay(int **, int , THREADED);
void enqueueAddOrdTaxa(int); 
void nQ4RandOpt(int **, int , THREADED);
int **genIstOptSubTask(int **, THREADED);
void displayTree(int **, int , int, THREADED);
void computeEckDayCost (struct subTask *, THREADED);
void genTaskLeft(int **, int **, int, int, int *, THREADED);
void genTaskRight(int **, int **, int, int, int *, THREADED);
void correctEckDayIndices (struct subTask *, THREADED);
void displayTreeInNexusFormat(int **, int , int , THREADED);
struct maxMiniTaxaPosWLength getEckDayMinTree(struct subTask *, int , THREADED);

int    pop(void);
void   push(int);
void   freeMemory(THREADED);
void   enqueueTaxa(void);
void   randOptimize(THREADED);
void   reorderSites(THREADED);
void   taxaAddOrder(THREADED);
void   showBits(int);
void   getIstBestCost(THREADED);
void   branchAndBound(THREADED);
void   deleteTaxa(int);
struct taskNode *popTask(void);
void   preProcess(char *, THREADED);
void   flushBestCostStack(THREADED);
struct subTask* dequeueTask(THREADED);
void   setFrameWork(char *, THREADED);
void   tbr(struct subTask *, THREADED);
struct subTask *genIstSubTask(THREADED);
void   checkUsage(int, char **, THREADED);
int    maskStateInPos(int , int);
void   freeTree(struct subTask *);
void   storeTree(struct subTask *, THREADED);
void   pushTask(struct taskNode *);
void   enqueueTask(struct subTask *, THREADED);
int    dequeueTaxaQ(struct taskNode *);
void   printTree(struct subTask *, int, THREADED);
void   freeStackNode(struct taskNode *, THREADED);
int    purdomsLowerBnd(struct taskNode *);
void   pushBestCostNode(struct subTask *);
void   subTaskPrint(struct subTask *, int, THREADED);
void   printTaxaStatus (struct taskNode *, THREADED);
struct taskNode *copyTaxaQ(struct taskNode *);
void   computeCost(struct subTask *, int, int, THREADED);
void   setParentChildIndices (struct subTask *);
void   printBestSolutions(THREADED);
struct maxMiniTaxaPosWLength maxMini(struct subTask *, THREADED);
struct maxMiniTaxaPosWLength smartGreedy(struct subTask *, THREADED);
struct maxMiniTaxaPosWLength getMinTree(struct subTask *, int , THREADED);
struct subTask *genNextTask(struct subTask *, int , int , THREADED);
struct subTask *genNextEckDayTask(struct subTask *, int , int , THREADED);
void performFitchOp(struct subTaskArrayNode *, struct subTaskArrayNode *,   \
                                          struct subTaskArrayNode *, THREADED);

FILE *fp;

#define MAXMINI
#define ROOT             -9999999
#define LEAF             -8888888
#define TRUE             1
#define FALSE            0
#define MAXSTATES        32
#define MAX_RAND         2147483647
#define DBL_INF          (double)MAX_RAND
#define TAXA_NAME_LENGTH 50

//#define TBR_VERY_VERBOSE
//#define TBR_VERIFY
//#define TBR_VERBOSE
//#define BnB_VERBOSE
//#define REORDER_VERBOSE
//#define MAX_MINI_VERBOSE
//#define BnB_VERY_VERBOSE
//#define BnB_VERY_VERY_VERBOSE
//#define GET_ALL_SOLUTIONS
#define GET_ONE_BEST_SOLUTION
//#define PREPROCESS_VERBOSE
//#define SMART_GREEDY_VERBOSE
//#define GET_MIN_TREE_VERBOSE
//#define TAXA_ADD_ORDER_VERBOSE
//#define FITCH_OPERATION_VERBOSE
//#define GET_IST_BESTCOST_VERBOSE

#endif

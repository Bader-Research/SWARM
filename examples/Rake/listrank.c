/**********************************************************************************
 NOTES:  This listranking code is modified to return the prefix value of the last element in the list after the List Ranking operation is done. This prefix value is the total no of leaves for the rake operation. It does not return the list head pointer.

Last changed July 30 2002

**********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "listrank.h"
#include "swarm_random.h"


typedef struct master_d
{
  LISTRANK_TYPE prefix;
  LDATA succ;
  LDATA index;
} master_t;


LDATA list_ranking(LDATA n, int k, list_t *List, THREADED)
{
  
  register LDATA i, j, s, target, current;
  register LISTRANK_TYPE val;
  LDATA block, start, finish, group, group_log, group_m1, times,
    shift, div, rem, v, *succ_ptr, succ, initial;
  master_t *Master;
  LDATA *Succ;
  master_t *s_master_ptr;
  list_t *s_list_ptr, *f_list_ptr;
  LDATA list_head_ptr=0;
  //added to listrank.c code to determine total no leaves in the tree.
  LDATA *totalnoleaves;

  SWARM_random_init(TH);
  SWARM_srandom(MYTHREAD+1,TH);

  s = k * THREADS;
  
  Master = (master_t *) SWARM_malloc((s+1)*sizeof(master_t),TH); 
  Succ = (LDATA *) SWARM_malloc(n*sizeof(LDATA),TH);
  totalnoleaves = (LDATA *) SWARM_malloc(sizeof(LDATA),TH);

  SWARM_Barrier();
  
  block = n/THREADS;
  start = block * MYTHREAD;
  finish = start + block;
  if (MYTHREAD == THREADS-1) finish = n;

  for (i=start ; i<finish ; i++) {
    Succ[i] = List[i].succ;

  
    if (List[i].succ < 0) {
      List[i].succ = -(s + 1);
    }
  }
  
  /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  /* IDENTIFY THE  MASTERS */
  
  s_master_ptr = Master + k * MYTHREAD;
  
  group = n/s;

	#ifdef HAVE_LOG2
  		group_log = (LDATA) floor(log2((double) group)+0.1);
	#else
  		group_log = (LDATA) floor((log((double) group)/log((double) 2))+0.1);
	#endif

  group_m1 = group - 1;
  
  times = (LDATA) 31/group_log;
  shift = 31 - times*group_log;
  div = k/times;
  rem = k - (times*div);
  start = start - group;
  current = -(k*MYTHREAD);
  
  for (i=0 ; i<div ; i++) {
    v = SWARM_random(TH);
    target = ((v = (v >> shift)) & group_m1) + (start += group);
    (++s_master_ptr) -> index = target;
    s_master_ptr -> succ = List[target].succ;
    List[target].succ = --current; 
    for (j=1;j<times;j++) {
      target = ((v >>= group_log) & group_m1) + (start += group);
      (++s_master_ptr) -> index = target;
      s_master_ptr -> succ = List[target].succ;
      List[target].succ = --current;
    }
  }
  
  v = SWARM_random(TH);
  v = (v >> shift);
  for (i=0;i<rem;i++) {
    target = ((v >>= group_log) & group_m1) + (start += group);
    (++s_master_ptr) -> index = target;
    s_master_ptr -> succ = List[target].succ;
    List[target].succ = --current;
  }
  
  SWARM_Barrier();

  on_one {
    if (List[list_head_ptr].succ < 0) {
      Master -> succ = initial = -List[list_head_ptr].succ;
      Master -> prefix = Master[initial].prefix = LISTRANK_IDENTITY;
    }
    else {
      initial = 0;
      Master -> prefix = LISTRANK_IDENTITY;
      Master -> index = list_head_ptr;
      Master -> succ = List[list_head_ptr].succ;
      List[list_head_ptr].succ = 0;
    } 
  }
  
  SWARM_Barrier();

  /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  /* TRAVERSE THE SUBLISTS */

  start = (k * MYTHREAD) + 1;
  finish = start + k;
  if ((MYTHREAD == 0) && (initial == 0))
    start = 0;
  
  for (i = start; i < finish; i++) {
    current = Master[i].index;
    val     = List[current].prefix;    
    current = Master[i].succ;
    if (current >= 0) { /* We aren't at the end of the list */
      target = List[current].succ;
      while (target >= 0) {
        val = List[current].prefix = (val LISTRANK_OPERATOR List[current].prefix);
	List[current].succ = - i;
	current            = target;
	target             = List[current].succ;
      }
      if (target > (- s - 1)) { /* We are at a new sublist */
        Master[- target].prefix = val;
	Master[i].succ =  -target;
      }
      else { /* We are at the end of the list */
        List[current].prefix = (val LISTRANK_OPERATOR List[current].prefix);
	List[current].succ = - i;
      }
    }
  } 
  
  SWARM_Barrier();

  /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */  
  /* RANK THE MASTERS: */

  on_one {
    succ = Master -> succ; 
    current = 0;
      
    for (i=0 ; i<s ; i++) {
      Master[succ].prefix =
	(Master[succ].prefix LISTRANK_OPERATOR Master[current].prefix);
      current = succ;
      succ = Master[current].succ;  
    } 
      
  }
  
  SWARM_Barrier();

  /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  /* COMPLETE THE RANK OF EVERY ELEMENT: */

  s_list_ptr = List + (block * MYTHREAD) -1;
  f_list_ptr = s_list_ptr + block;
  if (MYTHREAD == THREADS - 1) f_list_ptr = List + n - 1;
  succ_ptr = Succ + (block * MYTHREAD) - 1;
  
  while ((++s_list_ptr) <= f_list_ptr) {
    s_list_ptr->prefix = s_list_ptr->prefix LISTRANK_OPERATOR 
      Master[-s_list_ptr->succ].prefix;
    s_list_ptr -> succ =  *(++succ_ptr);
    //added to listrank.c code to determine total no leaves in the tree.
    if(s_list_ptr -> succ < 0)
      *totalnoleaves = s_list_ptr->prefix;
  } 

  SWARM_Barrier();
  SWARM_free(Master, TH);
  SWARM_free(Succ, TH);

  //return(SWARM_Bcast_l(list_head_ptr,TH));
  //added to listrank.c code to return total no leaves in the tree.
  return(*totalnoleaves);   
}





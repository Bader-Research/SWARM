/**********************************************************************************
 NOTES:   1. For listranking No of sublists (s) = k * THREADS. k=4 for FULL trees and k=8 for CAT and RAN trees. Hence make sure that the no of vertices >= sublists. Accordingly choose values of No of threads t and the no of vertices in the tree.More specifically,
t< no_of_vertices in the tree /4 for FULL trees or
t< no_of_vertices in the tree /8 for CAT/RAN trees 
          2. root is always in position 0 of array
	  3. each node has a sibling (except root)
	  4. a node has either no children or 2 children
	  5. long is used to deal with large trees (size >2GB).
	  6. For input levels larger than 27 typedef DATA as long in rake.c as well as typedef LDATA as long in listrank.h

Last changed Sep 11 2002

**********************************************************************************/

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include "swarm.h"
#include "swarm_random.h"
#include "listrank.h"

#define NIL -1
#define ALLONE 1
#define SERIAL 1
#define PARALLEL 1
#define ADD 0
#define MULT 1
#define TRUE 1
#define FALSE 0

typedef int DATA;

struct node{
  DATA a;			/* a,b are label values, c is leaf value */
  DATA b;
  DATA c;			
  int opp;                      /* opp == (+ / *) for internal nodes only */
  int LC;			/* LC is true if node is a left child */
  DATA index_left;
  DATA index_right;
  DATA index_parent;
  DATA index_me;
};

typedef struct node NODE;
typedef NODE * NODE_PTR;

/* serial solution for expression evaluation */
static DATA serial_express_eval(NODE* array,DATA vertices);

/* initializes non-root nodes for a FULL Binary tree */
static void internal(NODE *array, DATA positon, DATA vertices, DATA internal_nodes, 
	      int ran, THREADED);

/* creates a FULL, CAT or RANDOM binary tree*/
static void init_array(NODE *array, DATA vertices, int levels, int ran, int type, 
		THREADED);

/* evaluates expression like algorithm in JaJa's Book with data copy (step 2.3)*/
static DATA express_eval(NODE *array,int levels, DATA * leaves, NODE_PTR *A, NODE_PTR *B, 
		  THREADED);

/* evaluates expression with no data copy (step 2.3 is omitted). This is inefficient as compared to express_eval and hence is not used.*/
static DATA express_eval_nodatacopy(NODE *array,DATA leaves,int levels,NODE_PTR *A,THREADED);
 				
/* RAKE OPPS */
static void Rake_LC(DATA i, NODE_PTR *A, NODE *array);
static void Rake_RC(DATA i, NODE_PTR *A, NODE *array);
				
/* creates array with info about tree to pass to listranking */
static void crt_edg_arry(NODE *array, list_t *edg_arry, DATA vertices, THREADED);
				
/* used in testing, makes all the leaves hold a value of 1*/
static void make_c_1(NODE *array, DATA vertices,THREADED);
				
/* in tree creation->operation value.can be all 1-*, all 0 - +, or random*/
static int crt_opp(int ran, THREADED);


static void *Rake_main(THREADED)
{
  NODE *array;
  NODE_PTR *A, *B;
  double start_time1,start_time2,euler_time,lr_time;
  double rake_time,serial_time, leafarray_time,algo_time,speedup;
  list_t *edg_arry;
  DATA i,*leaves,vertices,serialresult,parallelresult,v=1,lrnoelems;
  int levels;
  FILE *fp;
				
  /* If you change this enum, change it in crt_opp too */
  enum RAN {PSEUDORAN, ALLMULT, ALLADD};
  int ran = ALLADD;
				
  /* FULL is a full binary tree. CAT is a caterpillar binary tree. 
     RANDOM is a random binary tree.
     If you change this enum change it in init_array too. */
  
  enum TYPE {FULL, CAT, RANDOM};
  int type;
    
  leaves=(DATA *)SWARM_malloc(sizeof(DATA),TH);
  fp=stdout;

  /* error reporting if argument is not specified */
  on_one_thread {
    if (THARGC != 2) {
      fprintf(stderr,"ERROR: Incorrect Number of arguments. First Argument : type of tree required (FULL:0 CAT:1 RANDOM:2), Second Argument : log2(|v|)\n");
      exit(-1);
    }
  }

  SWARM_Barrier();

  type      = atoi(THARGV[0]);
  levels    = atoi(THARGV[1]);
  vertices  = (v<<levels)-1;
  lrnoelems = (vertices-1)<<1;

  edg_arry = (list_t *)SWARM_malloc(sizeof(list_t)*lrnoelems,TH);
  array = (NODE *)SWARM_malloc(sizeof(NODE)*vertices, TH);

  SWARM_Barrier();

  init_array(array, vertices, levels, ran, type, TH);  

#if ALLONE
  make_c_1(array,vertices,TH);
#endif
 
#if SERIAL  
  on_one_thread {
    start_time1 = get_seconds();
    serialresult=serial_express_eval(array,vertices);
    serial_time=get_seconds()-start_time1;
  }
#endif

  SWARM_Barrier();

#if PARALLEL
  on_one_thread {
    start_time1 = start_time2= get_seconds();
  }

  /******************** Computes Euler Tour of the tree ************************/
  
  crt_edg_arry(array, edg_arry,vertices,TH);  

  on_one_thread {
    euler_time = get_seconds()-start_time2;
  }

  on_one_thread {
    start_time2=get_seconds();
  } 

  /************** List Ranking to store leaves in order in Array A ************/
  switch(type) {
  case FULL:
    *leaves=list_ranking(lrnoelems, 4, edg_arry, TH);
    break;
  default:
    *leaves=list_ranking(lrnoelems, 8, edg_arry, TH);
    break;
  }

  on_one_thread {
    lr_time = get_seconds()-start_time2;
  }

  on_one_thread {
    start_time2 = get_seconds();
  } 

  /********************* Allocating arrays A and B ***************************/

  A = (NODE_PTR *)SWARM_malloc(sizeof(NODE_PTR)*(*leaves),TH);
  B = (NODE_PTR *)SWARM_malloc(sizeof(NODE_PTR)*(*leaves),TH);
 

  /*****************************************************************************  
  Filling array A with pointers to the leaves in the tree. The leaves are 
  put in order (thanks to numbering the leaves with listranking).
  *****************************************************************************/  

  SWARM_pardo(i,1,vertices,1) {
    if (array[i].opp==NIL)
      A[edg_arry[(i<<1)-2].prefix] = &array[i];
  }
  SWARM_Barrier();

  on_one_thread {
    leafarray_time= get_seconds()-start_time2;
  }

  on_one_thread{
    start_time2 = get_seconds();
  }

  /* The expression evaluation function. When express_eval is called, B is empty.*/
 
  parallelresult = express_eval(array,levels, leaves, A, B, TH);

  /*Expression Evaluation with no data copy. 
    Not used because it is inefficient as compared to the 
    express_eval function which has data copy */

  /* parallelresult=express_eval_nodatacopy(array,*leaves,levels,A,TH); */
 
  on_one_thread {
    rake_time = get_seconds()-start_time2;
    /* algo_time = get_seconds()-start_time1; */
    algo_time = euler_time+lr_time+leafarray_time+rake_time;
    speedup=serial_time/algo_time;

    if (serialresult!=parallelresult) {
      fprintf(fp,"%d %d %d ERROR: Serial Result: %d  Parallel Result: %d\n",type,THREADS,levels,serialresult,parallelresult);
    }
    else {
      /*fprintf(fp,"Serial Result: %ld  Parallel Result: %ld\n",serialresult,parallelresult);*/
      fprintf(fp,"%d %d %d %11.7f  %11.7f  %11.7f  %11.7f  %11.7f  %11.7f %11.7f\n",type,THREADS,levels,serial_time,euler_time,lr_time,leafarray_time,rake_time,algo_time,speedup);
    }  
  }  

#endif

  fflush(fp);
  SWARM_free(array,TH);
  SWARM_free(A,TH);
  SWARM_free(B,TH);
  SWARM_free(edg_arry,TH); 
  SWARM_free(leaves,TH);
  fclose(fp);

  SWARM_Barrier();
  SWARM_done(TH);
}


static void init_array(NODE *array, DATA vertices, int levels, 
		int ran, int type, THREADED)
{

  DATA i,hit,parent,internal_nodes,in=1;
  long r;
 
  /* if you change this enum change it in main too */
  enum TYPE {FULL, CAT, RANDOM};

  switch(type) {
  case FULL:
    internal_nodes=(in<<(levels-1))-1;
    on_one_thread {
      array[0].a=1;
      array[0].b=0;
      array[0].c=NIL;
      array[0].opp=crt_opp(ran,TH);
      array[0].index_left=1;
      array[0].index_right=2;
      array[0].index_parent=0;
      array[0].LC=0;
      array[0].index_me=0;
    }
    SWARM_pardo (i, 1, vertices, 1) {
      internal(array, i, vertices, internal_nodes, ran, TH);
    }
    
    SWARM_pardo (i, 1, vertices-1, 2) {
      array[i].LC=1;
      array[i+1].LC=0;
    }
    break;
   
  case CAT:
    on_one_thread {
      array[0].a=1;
      array[0].b=0;
      array[0].c=NIL;
      array[0].opp=crt_opp(ran,TH);
      array[0].index_left=1;
      array[0].index_right=2;
      array[0].index_parent=0;
      array[0].LC=0;
      array[0].index_me=0;

      array[1].a=1;
      array[1].b=0;
      array[1].c=NIL;
      array[1].opp=crt_opp(ran,TH);
      array[1].index_left=3;
      array[1].index_right=4;
      array[1].index_parent=0;
      array[1].index_me=1;  
      array[1].LC=1;

      array[2].a=1;
      array[2].b=0;
      array[2].c=2%4 + 1;
      array[2].opp=NIL;
      array[2].index_left=NIL;
      array[2].index_right=NIL;
      array[2].index_parent=0;
      array[2].index_me=2;
      array[2].LC=0;

    }
    pardo(i,3,vertices-2,2) {
      array[i].a=1;
      array[i].b=0;
      array[i].c=NIL;
      array[i].opp=crt_opp(ran,TH);
      array[i].index_left=2+i;
      array[i].index_right=3+i;
      array[i].index_parent=i-2;
      array[i].index_me=i;  
      array[i].LC=1;

      array[i+1].a=1;
      array[i+1].b=0;
      array[i+1].c=(i+1)%4 + 1;
      array[i+1].opp=NIL;
      array[i+1].index_left=NIL;
      array[i+1].index_right=NIL;
      array[i+1].index_parent=i-2;
      array[i+1].index_me=i+1;
      array[i+1].LC=0;
    }
    on_one_thread {
      array[vertices-2].a=1;
      array[vertices-2].b=0;
      array[vertices-2].c=(vertices-2)%4+1;
      array[vertices-2].opp=NIL;
      array[vertices-2].index_left=NIL;
      array[vertices-2].index_right=NIL;
      array[vertices-2].index_parent=vertices-4;
      array[vertices-2].index_me=vertices-2;
      array[vertices-2].LC=1;

      array[vertices-1].a=1;
      array[vertices-1].b=0;
      array[vertices-1].c=(vertices-1)%4+1;
      array[vertices-1].opp=NIL;
      array[vertices-1].index_left=NIL;
      array[vertices-1].index_right=NIL;
      array[vertices-1].index_parent=vertices-4;
      array[vertices-1].index_me=vertices-1;
      array[vertices-1].LC=0;
    }
    break;

  case RANDOM:
    on_one_thread {
      array[0].a=1;
      array[0].b=0;
      array[0].c=NIL;
      array[0].opp=crt_opp(ran,TH);
      array[0].index_left=1;
      array[0].index_right=2;
      array[0].index_parent=0;
      array[0].LC=0;
      array[0].index_me=0;

      array[1].a=1;
      array[1].b=0;
      array[1].c=1%4+1;
      array[1].opp=NIL;
      array[1].index_left=NIL;
      array[1].index_right=NIL;
      array[1].index_parent=0;
      array[1].index_me=1;  
      array[1].LC=1;

      array[2].a=1;
      array[2].b=0;
      array[2].c=2%4 + 1;
      array[2].opp=NIL;
      array[2].index_left=NIL;
      array[2].index_right=NIL;
      array[2].index_parent=0;
      array[2].index_me=2;
      array[2].LC=0;

      for(i=3;i<vertices-1;i=i+2) {
	r = random();
	hit=((int)r)%(i-2)+2;
	parent=array[hit].index_parent;
	if (array[parent].index_left == hit)
	  array[parent].index_left=i;
	else
	  array[parent].index_right=i;

	array[i].a=1;
	array[i].b=0;
	array[i].c=NIL;
	array[i].opp=crt_opp(ran,TH);
	array[i].index_me=i;
	array[i].LC=array[hit].LC;
	array[i].index_parent=parent;

	if (array[hit].LC) {
	  array[i].index_left=hit;
	  array[i].index_right=i+1;
	}
	else {
	  array[i].index_left=i+1;
	  array[i].index_right=hit;
	}

	array[hit].index_parent=i;

	array[i+1].a=1;
	array[i+1].b=0;
	array[i+1].c=(i+1)%4 + 1;
	array[i+1].opp=NIL;
	array[i+1].index_left=NIL;
	array[i+1].index_right=NIL;
	array[i+1].index_parent=i;
	array[i+1].index_me=i+1;
	if (array[hit].LC)
	  array[i+1].LC=0;
	else
	  array[i+1].LC=1;
      }

    }
    break;

  default:
    fprintf(stdout,"Error: Not a defined structure type\n");
    fflush(stdout);
    exit(-1);
  }

  SWARM_Barrier();  
} /* end init array*/

/* this function creates both internal and leaf nodes for full binary trees */

static void internal(NODE *array, DATA positon, DATA vertices, DATA internal_nodes, 
	      int ran, THREADED)
{
  if (positon>=internal_nodes && positon<vertices) { /* Ifleaf */
    array[positon].a=1;
    array[positon].b=0;
    array[positon].c=positon%4 + 1;
    array[positon].opp=NIL;
    array[positon].index_left=NIL;
    array[positon].index_right=NIL;
    array[positon].index_parent=(positon-1)/2;
    array[positon].index_me=positon;
  }
  else { /* If internal nodes*/
    array[positon].a=1;
    array[positon].b=0;
    array[positon].c=NIL;
    array[positon].opp=crt_opp(ran,TH);
    array[positon].index_left=2*positon+1;
    array[positon].index_right=2*positon+2;
    array[positon].index_parent=(positon-1)/2;
    array[positon].index_me=positon;
  }  
}


static DATA express_eval(NODE *array,int levels, DATA *leaves, NODE_PTR *A, NODE_PTR *B, 
		  THREADED) {

  int i, flag=1;
  DATA j,final_value,newleaves=0;
  NODE *n1,*left,*right;
  
  for (i=0; i<levels - 1; i++) {
    if (flag) {
      /*don't want the first or last leaf */
      SWARM_pardo(j,1,*leaves-1, 2) { 
	if (A[j]->LC==TRUE)
	  Rake_LC(j, A, array);
      }
      SWARM_Barrier();
	  
      SWARM_pardo(j,1,*leaves-1, 2) { 
	if (A[j]->LC==FALSE)
	  Rake_RC(j, A, array);
      } 
      SWARM_Barrier();

      newleaves=((*leaves-2)>>1)+1;	  
      /* copy the non-raked leaves from A to B array. first and last leaf special*/
      SWARM_pardo(j,1,newleaves, 1) {
	B[j]=A[(j<<1)];
      }
      SWARM_Barrier();

      /* copy the first and last leaf.*/
      on_one_thread {
	B[0]=A[0];
	B[newleaves]=A[*leaves-1];
	*leaves=newleaves+1;
      }
      SWARM_Barrier();
	  
    } /* end (flag) */
      
    if (!flag) {
      /*don't want the first or last leaf */
      SWARM_pardo(j,1,*leaves-1, 2) { 
	if (B[j]->LC==TRUE)
	  Rake_LC(j, B, array);
      }	  
      SWARM_Barrier();
	
      SWARM_pardo(j,1,*leaves-1, 2) { 
	if (B[j]->LC==FALSE)
	  Rake_RC(j, B, array);
      }	  
      SWARM_Barrier();
	
      newleaves=((*leaves-2)>>1)+1;	   
      /*copy the non-raked leaves from B array to A array. first and last leaf special*/
      SWARM_pardo(j,1,newleaves, 1) {
	A[j]=B[(j<<1)];
      }
      SWARM_Barrier();

      /*copy the first and the last leaf*/
      on_one_thread {
	A[0]=B[0];
	A[newleaves]=B[*leaves-1];
	*leaves=newleaves+1;
      }
      SWARM_Barrier();
    } /*end (!flag) */

    flag=1-flag;
  }/* end for */

  n1=array;
  left=n1+n1->index_left;
  right=n1+n1->index_right;

  if (n1->opp==MULT)
    final_value= (left->a * left->c + left->b) * (right->a * right->c + right->b); 
  else   /*operator == addition*/
    final_value= (left->a * left->c + left->b) + (right->a * right->c + right->b); 
    
  return final_value;
    
} /* end express_eval*/

static DATA express_eval_nodatacopy(NODE *array,DATA leaves,int levels,NODE_PTR *A,THREADED){
  int i;
  DATA itrincr,j,final_value;
  NODE *n1;
  NODE *left,*right;

  for (i=0; i<levels - 1; i++) {
    itrincr=(1<<i);
    pardo(j,itrincr,leaves-1,(itrincr<<1)) {
      if (A[j]->LC==TRUE)
        Rake_LC(j, A, array);
    }
    SWARM_Barrier();
    pardo(j,itrincr,leaves-1,(itrincr<<1)) {
      if (A[j]->LC==FALSE)
        Rake_RC(j, A, array);
    }
    SWARM_Barrier();
  }

  n1=array;
  left=n1+n1->index_left;
  right=n1+n1->index_right;

  if (n1->opp==MULT)
    final_value= (left->a * left->c + left->b) * (right->a * right->c + right->b);
  else   /*operator == addition*/
    final_value= (left->a * left->c + left->b) + (right->a * right->c + right->b);

  return final_value;

} /* end express_eval_nodatacopy*/




static void Rake_LC(DATA i, NODE_PTR *A, NODE *array)
{
  DATA parent, parent_a, parent_b, index_me, gparent, index_sib, sib_a, sib_b;
  NODE *parentnode,*mynode,*sibnode;

  parent=A[i]->index_parent;
  parentnode=array+parent;
  parent_a=parentnode->a;
  parent_b=parentnode->b;

  index_me=A[i]->index_me;
  mynode=array+index_me;

  gparent=parentnode->index_parent;

  index_sib=parentnode->index_right;
  sibnode=array+index_sib;
  sib_a=sibnode->a;
  sib_b=sibnode->b;

  
  if (parentnode->opp == MULT) {
    sibnode->a=parent_a * (mynode->a * mynode->c + mynode->b) * sib_a;
    sibnode->b=parent_a * (mynode->a * mynode->c + mynode->b) * sib_b + parent_b; 
  }
  else {   /*operator == ADDITION */
    sibnode->a=parent_a * sib_a;
    sibnode->b=parent_a * (mynode->a * mynode->c + mynode->b + sib_b )+ parent_b;
  }

  sibnode->LC=parentnode->LC;

  if (parentnode->LC)
    array[gparent].index_left=index_sib;
  else
    array[gparent].index_right=index_sib;

  sibnode->index_parent=gparent;

} /* end Rake_LC */

static void Rake_RC(DATA i, NODE_PTR *A, NODE *array)
{
  DATA parent, parent_a, parent_b, index_me, gparent, index_sib, sib_a, sib_b;
  NODE *parentnode,*mynode,*sibnode;

  parent=A[i]->index_parent;
  parentnode=array+parent;
  parent_a=parentnode->a;
  parent_b=parentnode->b;

  index_me=A[i]->index_me;
  mynode=array+index_me;

  gparent=array[parent].index_parent;
  
  index_sib=array[parent].index_left;
  sibnode=array+index_sib;
  sib_a=sibnode->a;
  sib_b=sibnode->b;

  if (parentnode->opp == MULT) {
    sibnode->a=parent_a * (mynode->a * mynode->c + mynode->b) * sib_a;
    sibnode->b=parent_a * (mynode->a * mynode->c + mynode->b) * sib_b + parent_b; 
  }
  else { /*operator == ADDITION	*/
    sibnode->a=parent_a * sib_a;
    sibnode->b=parent_a * (mynode->a * mynode->c + mynode->b + sib_b )+ parent_b;
  }

  sibnode->LC=parentnode->LC;

  if (parentnode->LC)
    array[gparent].index_left=index_sib;
  else
    array[gparent].index_right=index_sib;

  sibnode->index_parent=gparent;
  
} /* end Rake_RC */


static void crt_edg_arry(NODE *array, list_t  *edg_arry, DATA vertices, THREADED)
{
  DATA i,parent;
  list_t *e1,*e2;
  NODE *n1;

  pardo (i,1,vertices, 1) { /*excludes root of the tree*/
    e1=edg_arry+(i<<1)-1;
    e2=e1-1;
    n1=array+i;

    parent=n1->index_parent;

    if (n1->opp==NIL) { /* if leaf */
      e2->prefix = 0;
      e1->prefix = 1;
	  
      e2->succ = (i<<1)-1;
    }
    else { /* if not leaf or root */
      e2->prefix = 0;
      e1->prefix = 0;
	  
      e2->succ = ((n1->index_left)<<1) - 2;
    }

    if (n1->LC) /* if left child*/
      e1->succ = ((array[parent].index_right)<<1) - 2;
      
    else if (array[parent].index_parent!=parent) /* Right child and parent is not root*/
      e1->succ = (parent<<1)-1;
      
    else if (array[parent].index_parent==parent) 
      e1->succ = -1;

  }
  SWARM_Barrier();
}

/* This function makes all the leaves have a value of 1 for Debugging */

static void make_c_1(NODE *array, DATA vertices,THREADED)
{
  DATA i;
  pardo(i,1,vertices,1) {
    if (array[i].index_left==NIL && array[i].index_right==NIL) {
      array[i].c=1;
    }
  }
  SWARM_Barrier();
}

/* This functions assigns the value of the operator for the internal nodes */

static int crt_opp(int ran, THREADED)
{
  /* if you change this change in main too */ 
  enum RAN {PSEUDORAN, ALLMULT, ALLADD};

  SWARM_srandomBit((MYTHREAD+17)*53, TH);
  
  switch(ran) {
  case PSEUDORAN:
    return SWARM_randomBit(TH);
  case ALLMULT:
    return 1;
  case ALLADD:
    return 0;
  default:
    printf("Error: Random type not defined\n");
    fflush(stdout);
  }
  return 0;
}

/* this function is for the serial version of the problem. */

static DATA serial_express_eval(NODE* array,DATA vertices){

  struct stack{
    DATA top;
    DATA *item;
  } s;

  DATA arrsize,i;
  DATA currentindex,*expression,index=0,serialresult;
  arrsize=vertices*sizeof(DATA);

  s.item = (DATA *) malloc(arrsize);
  expression = (DATA *) malloc(arrsize);

  s.top=-1;
  currentindex=0;

  do {
    while (currentindex!=-1) {
      s.item[++s.top]=currentindex;
      currentindex= array[currentindex].index_left;
    }
    if (s.top!=-1) {
      currentindex=s.item[s.top--];

      if (array[currentindex].opp!=NIL)
	expression[index++]=(DATA) array[currentindex].opp;
      else
	expression[index++]=array[currentindex].c;

      currentindex=array[currentindex].index_right;
    }
  } while(s.top!=-1 || currentindex !=-1);

  serialresult=expression[0];
  for (i=1;i<vertices;i=i+2) {
    if (expression[i]==0)
      serialresult=serialresult+expression[i+1];
    else
      serialresult=serialresult*expression[i+1];
  }

  return serialresult;
}


int main(int argc, char **argv) 
{
  SWARM_Init(&argc,&argv);
  SWARM_Run((void *)Rake_main);
  SWARM_Finalize();
  return 0;
}

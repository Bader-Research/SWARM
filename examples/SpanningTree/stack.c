#include "stack.h"

void push(int a, int * stack, int *top)
{
   (*top)++;
   stack[*top]=a;
}

#if 0
int pop(int * stack, int *top,int bottom)
{
   int a ;
   if(*top <=bottom || *top<0) return -1;
   else{
  	 a= stack[*top];
   	(*top)--;
   	return a;
   }

}
#endif

int pop(int * stack, int * top, int bottom)
{
   int a = stack[*top];
   (*top)--;
    return a;
}

int is_empty(int * stack, int * top, int * bottom)
{
  if(bottom){
   if((*top)<=(*bottom) || (*top)<0) return 1;
   else return 0;
  } else {
    if((*top)==-1) return 1;
    else return 0;
  }
}


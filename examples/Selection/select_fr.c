#include "smp_sel.h"


void partition_with_two_fr_i(int *A, int *B, int M, int s0, int s1,
			  int *a0, int *a1) {
  int *Aptr, *eptr;
  int *a0ptr, *a1ptr, *a2ptr;
  int *D;

  if (M>0) {
    D = (int *)malloc(M*sizeof(int));
    assert_malloc(D);

    a0ptr = B;
    a1ptr = D;
    a2ptr = B+ M-1;
  
    *a0 = 0;
    Aptr = A;
    eptr = Aptr + M;
    while (Aptr < eptr) {
      if (*Aptr<s0) {
	*a0ptr++ = *Aptr;
      }
      else {
	if (*Aptr<=s1) {
	  *a1ptr++ = *Aptr;
	}
	else {
	  *a2ptr-- = *Aptr;
	}
      }
      Aptr++;
    }

    *a0 = a0ptr-B;
    *a1 = a1ptr-D;

#if 0
    if ((B + *a0 + *a1 - 1) != a2ptr) {
      fprintf(outfile,"a0: %8d a1: %8d a2: %8d sum: %8d  M: %8d\n",
	      *a0, *a1, A+M-a2ptr, *a0+*a1+(A+M-a2ptr), M);
    }
#endif
  
 #if 0
    bcopy(D, a0ptr, *a1*sizeof(int));
#else
    bcopy(B,       A            , *a0*sizeof(int));
    bcopy(D,       A + *a0      , *a1*sizeof(int));
    bcopy(a2ptr+1, A + *a0 + *a1, (M - (*a0+*a1))*sizeof(int));
#endif
    

    free(D);
  }
  else {
    *a0 = 0;
    *a1 = 0;
  }


  return; 
}

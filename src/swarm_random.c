#include <swarm_random.h>

#ifndef HAVE_SPRNG
/*
 * rrandom.c:
 *
 * An improved random number generation package.  In addition to the standard
 * rand()/srand() like interface, this package also has a special state info
 * interface.  The iinitstate() routine is called with a seed, an array of
 * bytes, and a count of how many bytes are being passed in; this array is
 * then initialized to contain information for random number generation with
 * that much state information.  Good sizes for the amount of state
 * information are 32, 64, 128, and 256 bytes.  The state can be switched by
 * calling the ssetstate() routine with the same array as was initiallized
 * with iinitstate().  By default, the package runs with 128 bytes of state
 * information and generates far better random numbers than a linear
 * congruential generator.  If the amount of state information is less than
 * 32 bytes, a simple linear congruential R.N.G. is used.
 *
 * Internally, the state information is treated as an array of longs; the
 * zeroeth element of the array is the type of R.N.G. being used (small
 * integer); the remainder of the array is the state information for the
 * R.N.G.  Thus, 32 bytes of state information will give 7 longs worth of
 * state information, which will allow a degree seven polynomial.  (Note:
 * the zeroeth word of state information also has some other information
 * stored in it -- see ssetstate() for details).
 * 
 * The random number generation technique is a linear feedback shift register
 * approach, employing trinomials (since there are fewer terms to sum up that
 * way).  In this approach, the least significant bit of all the numbers in
 * the state table will act as a linear feedback shift register, and will
 * have period 2^deg - 1 (where deg is the degree of the polynomial being
 * used, assuming that the polynomial is irreducible and primitive).  The
 * higher order bits will have longer periods, since their values are also
 * influenced by pseudo-random carries out of the lower bits.  The total
 * period of the generator is approximately deg*(2**deg - 1); thus doubling
 * the amount of state information has a vast influence on the period of the
 * generator.  Note: the deg*(2**deg - 1) is an approximation only good for
 * large deg, when the period of the shift register is the dominant factor.
 * With deg equal to seven, the period is actually much longer than the
 * 7*(2**7 - 1) predicted by this formula.
 */

/*
 * For each of the currently supported random number generators, we have a
 * break value on the amount of state information (you need at least this
 * many bytes of state info to support this random number generator), a degree
 * for the polynomial (actually a trinomial) that the R.N.G. is based on, and
 * the separation between the two lower order coefficients of the trinomial.
 */
#define TYPE_0          0               /* linear congruential */
#define BREAK_0         8
#define DEG_0           0
#define SEP_0           0

#define TYPE_1          1               /* x**7 + x**3 + 1 */
#define BREAK_1         32
#define DEG_1           7
#define SEP_1           3

#define TYPE_2          2               /* x**15 + x + 1 */
#define BREAK_2         64
#define DEG_2           15
#define SEP_2           1

#define TYPE_3          3               /* x**31 + x**3 + 1 */
#define BREAK_3         128
#define DEG_3           31
#define SEP_3           3

#define TYPE_4          4               /* x**63 + x + 1 */
#define BREAK_4         256
#define DEG_4           63
#define SEP_4           1
#endif

void SWARM_random_init(THREADED) {

#ifdef HAVE_SPRNG
  THSPRNG = init_sprng(SPRNG_LCG,
		       MYTHREAD, THREADS,
		       make_sprng_seed(),SPRNG_DEFAULT);
  return;
#else
#if defined(SOLARIS)
  return;
#else
  long *lp;
  
  THRAND.randtbl = (long *)malloc((DEG_3 + 1)*sizeof(long));
  assert_malloc(THRAND.randtbl);

  lp = THRAND.randtbl;
  *lp++ = TYPE_3;
  *lp++ = 0x9a319039; *lp++ = 0x32d9c024; *lp++ = 0x9b663182;
  *lp++ = 0x5da1f342; *lp++ = 0xde3b81e0; *lp++ = 0xdf0a6fb5;
  *lp++ = 0xf103bc02; *lp++ = 0x48f340fb; *lp++ = 0x7449e56b;
  *lp++ = 0xbeb1dbb0; *lp++ = 0xab5c5918; *lp++ = 0x946554fd;
  *lp++ = 0x8c2e680f; *lp++ = 0xeb3d799f; *lp++ = 0xb11ee0b7;
  *lp++ = 0x2d436b86; *lp++ = 0xda672e2a; *lp++ = 0x1588ca88;
  *lp++ = 0xe369735d; *lp++ = 0x904f35f7; *lp++ = 0xd7158fd6;
  *lp++ = 0x6fa6f051; *lp++ = 0x616e6b96; *lp++ = 0xac94efdc;
  *lp++ = 0x36413f93; *lp++ = 0xc622c298; *lp++ = 0xf5a42ab8;
  *lp++ = 0x8a88d77b; *lp++ = 0xf5ad9d0e; *lp++ = 0x8999220b;
  *lp++ = 0x27fb47b9; *lp = 0;

  THRAND.fptr = &(THRAND.randtbl[SEP_3 + 1]);
  THRAND.rptr = &(THRAND.randtbl[1]);
  THRAND.state = &(THRAND.randtbl[1]);
  THRAND.rand_type = TYPE_3;
  THRAND.rand_deg = DEG_3;
  THRAND.rand_sep = SEP_3;
  THRAND.end_ptr = &(THRAND.randtbl[DEG_3 + 1]);
  return;
#endif
#endif

}

void SWARM_random_destroy(THREADED) {

#ifdef HAVE_SPRNG
  free_sprng(THSPRNG);
  return;
#else
#if defined(SOLARIS)
  return;
#else
  free(THRAND.randtbl);
  return;
#endif
#endif

}

void SWARM_srandom(unsigned int x, THREADED) {

#ifdef HAVE_SPRNG
  return;
#else
#if defined(SOLARIS)
  ti->xi[0] = MYTHREAD+(x*x*x);
  ti->xi[1] = (MYTHREAD+17)*x;
  ti->xi[2] = x + (ti->xi[0] * ti->xi[0]) + (ti->xi[1] * ti->xi[1]);
  return;
#else
  register int i;
  
  if (THRAND.rand_type == TYPE_0)
    THRAND.state[0] = x;
  else {
    THRAND.state[0] = x;
    for (i = 1; i < THRAND.rand_deg; i++)
      THRAND.state[i] = 1103515245 * THRAND.state[i - 1] + 12345;
    THRAND.fptr = &(THRAND.state[THRAND.rand_sep]);
    THRAND.rptr = &(THRAND.state[0]);
    for (i = 0; i < 10 * THRAND.rand_deg; i++)
      (void)SWARM_random(TH);
  }
  return;
#endif
#endif

}

long SWARM_random(THREADED) {

#ifdef HAVE_SPRNG
  return isprng(THSPRNG);
#else
#if defined(SOLARIS)
  return(nrand48(ti->xi));
#else
  long i;

  *(THRAND.fptr) += *(THRAND.rptr);
  i = (*(THRAND.fptr) >> 1) & 0x7fffffff;  /* chucking least random bit */
  if (++(THRAND.fptr) >= THRAND.end_ptr) {
    THRAND.fptr = THRAND.state;
    ++(THRAND.rptr);
  } else if (++(THRAND.rptr) >= THRAND.end_ptr)
    THRAND.rptr = THRAND.state;
  
  return(i);
#endif
#endif

}

/*********************************************/
/* random bits */
/*********************************************/
#define _RANDOM_BITSTRING_SIZE 31

void SWARM_srandomBit(unsigned int x, THREADED) {
  SWARM_random_init(TH);
  SWARM_srandom(x, TH);
  ti->rc = _RANDOM_BITSTRING_SIZE+1;
  return;
}

int  SWARM_randomBit(THREADED) {
  long r;
  if (ti->rc >= _RANDOM_BITSTRING_SIZE) {
    ti->rc = 0;
    ti->rbs = SWARM_random(TH);
  }
  r = (ti->rbs) & 0x1;
  ti->rbs >>= 1;
  ti->rc++;
  return ((int)r);
}


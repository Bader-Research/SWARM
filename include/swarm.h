#ifndef _SWARM_H
#define _SWARM_H

#ifdef __cplusplus
extern "C" {
#endif 

#ifdef WIN32
#include <swarm_config_win.h>
#else
#include <swarm_config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <pthread.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <math.h>

#ifdef TIME_WITH_SYS_TIME 
#include <sys/time.h>       /* struct timeval */
#endif


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#ifdef HAVE_SPRNG_SPRNG_H
#include <sprng/sprng.h>
#else
#ifdef HAVE_SPRNG_H
#include <sprng.h>
#endif 
#endif

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif


#ifndef SWARM_DLL
#ifdef WIN32
#ifdef SWARM_BUILD
/* DLL export */
#define SWARM_DLL __declspec(dllexport)
#else
/* EXE import */
#define SWARM_DLL __declspec(dllimport)
#endif
#else
#define SWARM_DLL
#endif
#endif

enum reduce_tag {MAX, MIN, SUM, PROD, LAND, BAND, LOR, BOR, LXOR, BXOR};

extern SWARM_DLL FILE* SWARM_outfile;

#define DEFAULT_PRI (PRI_OTHER_MIN+PRI_OTHER_MAX)/2
#define pNULL       (NULL)
#define SWARM_done(x) SWARM_Cleanup(x); pthread_exit(pNULL); return 0;

#define pthread_mb_np() asm("mb;")

#define MYTHREAD     (ti->mythread)
#define THREADED     uthread_info_t *ti
#define TH           ti
#define on_one_thread if (MYTHREAD == 0)
#define on_thread(k)  if (MYTHREAD == (k))
#define on_one        on_one_thread
#define BIND_TH       0
#define ERR_TH        0
#define CACHELOG      7
#define NOSHARE(x)    ((x)<<CACHELOG)

#define THARGC       (ti->argc)
#define THARGV       (ti->argv)
#define EXENAME      (ti->argv[0])

#define THRAND       (ti->rand)

/* max and min are predefined under stdlib.h in Visual Studio */
#undef min		
#undef max
#define max(a,b)   ((a) > (b) ? (a) : (b))
#define min(a,b)   ((a) < (b) ? (a) : (b))

#ifdef HAVE_LIBSPRNG
#define THSPRNG     (ti->randgen)
#else
#if defined(SOLARIS)
#else
typedef struct {
  long *randtbl;
  long *fptr;
  long *rptr;
  long *state;
  int rand_type;
  int rand_deg;
  int rand_sep;
  long *end_ptr;
} rrandom_info_t;
#endif
#endif


#ifndef MAXTHREADS
extern int          MAXTHREADS;
#endif
extern SWARM_DLL int 	THREADS;

struct thread_inf {
  int                 mythread;   /* Thread number */
  int                 argc;      
  char              **argv;
  long                  m1;        /* used in pardo */
  long                  m2;        /* used in pardo */
  long                 blk;        /* used in pardo */
#ifdef HAVE_LIBSPRNG
  int                *randgen;    /* SPRNG generator pointer */
#else
#if defined(SOLARIS)
  unsigned short      xi[3];      /* used in rrandom_th nrand48() */
#else
  rrandom_info_t      rand;       /* used in rrandom_th */
#endif
#endif
  long                rbs;        /* used in random_bit, random_bitstring */
  short               rc;         /* used in random_bit, random_count */
  int                 udata;      /* User data */
};

typedef struct thread_inf uthread_info_t;

typedef int reduce_t;

extern uthread_info_t *uthread_info;

#define SWARM_partition_loop_across_threads(i,first,last,inc) \
        SWARM_pardo((i),(first),(last),(inc))     
#define pardo(i,first,last,inc) \
        SWARM_pardo((i),(first),(last),(inc))     
#define SWARM_pardo(i,first,last,inc)                    \
        if (((first)==0)&&((last)==THREADS)) {          \
            ti->m1 = MYTHREAD;                          \
	    ti->m2 = ti->m1 + 1;                        \
	} else {                                        \
            ti->blk = ((last)-(first))/THREADS;         \
            if (ti->blk == 0) {    		        \
              ti->m1  = (first)+MYTHREAD;               \
              ti->m2  = (ti->m1) + 1;                   \
              if ((ti->m1) >= (last))                   \
                 ti->m1 = ti->m2;                       \
	    }                                           \
            else {                                      \
              ti->m1  = (ti->blk) * MYTHREAD + (first); \
	      if (MYTHREAD < THREADS-1)                 \
	          ti->m2 = (ti->m1)+(ti->blk);          \
	      else                                      \
	          ti->m2 = (last);                      \
            }                                           \
	}                                               \
        if ((inc)>1) {                            \
            while ((ti->m1-(first)) % (inc) > 0)  \
                ti->m1 += 1;                      \
        }                                         \
	for (i=ti->m1 ; i<ti->m2 ; i+=(inc))

#define task_do(x)  (MYTHREAD == ((x) % THREADS))

#if 1
#define SWARM_Barrier() SWARM_Barrier_sync(TH)
#else
#define SWARM_Barrier() SWARM_Barrier_tree(TH)
#endif
  
void   	SWARM_DLL SWARM_Barrier_tree(THREADED);
void   	SWARM_DLL SWARM_Barrier_sync(THREADED);
void  	SWARM_DLL *SWARM_malloc(int bytes, THREADED);
void   	SWARM_DLL *SWARM_malloc_l(long bytes, THREADED);
void   	SWARM_DLL SWARM_free(void *, THREADED);

typedef pthread_mutex_t SWARM_mutex_t;
typedef pthread_mutexattr_t SWARM_mutexattr_t;
int    SWARM_mutex_init(SWARM_mutex_t **, const SWARM_mutexattr_t *, THREADED);
int    SWARM_mutex_destroy(SWARM_mutex_t *, THREADED);
#define SWARM_mutex_lock(m)    pthread_mutex_lock(m)
#define SWARM_mutex_trylock(m) pthread_mutex_trylock(m)
#define SWARM_mutex_unlock(m)  pthread_mutex_unlock(m)

int    SWARM_DLL SWARM_Bcast_i(int    myval, THREADED);
long   SWARM_DLL SWARM_Bcast_l(long   myval, THREADED);
double SWARM_DLL SWARM_Bcast_d(double myval, THREADED);
char   SWARM_DLL SWARM_Bcast_c(char   myval, THREADED);
int    SWARM_DLL *SWARM_Bcast_ip(int    *myval, THREADED);
long   SWARM_DLL *SWARM_Bcast_lp(long   *myval, THREADED);
double SWARM_DLL *SWARM_Bcast_dp(double *myval, THREADED);
char   SWARM_DLL *SWARM_Bcast_cp(char   *myval, THREADED);
int    SWARM_DLL SWARM_Reduce_i(int    myval, reduce_t op, THREADED);
long   SWARM_DLL SWARM_Reduce_l(long   myval, reduce_t op, THREADED);
double SWARM_DLL SWARM_Reduce_d(double myval, reduce_t op, THREADED);
int    SWARM_DLL SWARM_Scan_i(int    myval, reduce_t op, THREADED);
long   SWARM_DLL SWARM_Scan_l(long   myval, reduce_t op, THREADED);
double SWARM_DLL SWARM_Scan_d(double myval, reduce_t op, THREADED);

void SWARM_DLL SWARM_Init(int*, char***);
void SWARM_DLL SWARM_Run(void *);
void SWARM_DLL SWARM_Finalize(void);
void SWARM_DLL SWARM_Cleanup(THREADED);

void assert_malloc(void *ptr);
double SWARM_DLL get_seconds(void);

#ifndef HAVE_GETTIMEOFDAY
	struct timezone 
	{
		int tz_minuteswest; 
		int tz_dsttime;
	};
	struct timeval
	{
		long tv_sec;
		long tv_usec;
	};
	void gettimeofday(struct timeval* p, struct timezone* tz /* IGNORED */);
#endif

#ifndef HAVE_GETOPT_H     
	extern char *optarg;
	extern int optind;
	extern int opterr;

	#ifndef NULL
		#define NULL	0
	#endif
	#define EOF     (-1)
	#define ERR(str, chr) (opterr ? fprintf(stderr, "%s%s%c\n", argv[0], str, chr) : 0)
	int getopt(int argc, char **argv, char *opts);
#endif

#define errprnt(msg)    { fprintf(stderr,"%s: %s\n",EXENAME,msg); exit(1); }

#if 1

#ifdef SOLARIS
	#include <sys/isa_defs.h>
	#ifndef volatile
		#define volatile
	#endif
	#include <sys/errno.h>
	#if defined(__sparcv9)
		#define sun_mb_mi_np() asm("membar #MemIssue ;")
	#else
		#define sun_mb_mi_np() 
	#endif
#endif
#else
	#ifdef SOLARIS
		#define sun_mb_mi_np() asm("membar #MemIssue ;")
	#endif
#endif

#ifdef __cplusplus
}
#endif 

  
#endif

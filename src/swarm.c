
#include <swarm.h>
#ifdef RRANDOM
#include <swarm_random.h>
#endif /* RRANDOM */

#include <time.h>

#include <swarm_multicore.h>

#define  DEBUG  1
#define  INFO   1

#define MAXLEN  80

#define MAXTHREADS_DEFAULT 64

#ifndef MAXTHREADS
int     MAXTHREADS = MAXTHREADS_DEFAULT;
#endif /* MAXTHREADS */
#define DEFAULT_THREADS 2
int     THREADS;
uthread_info_t *uthread_info;
static pthread_t     *spawn_thread;

static int    _swarm_init=0;

#define MAX_GATHER 2

static int    _SWARM_bcast_i;
static long   _SWARM_bcast_l;
static double _SWARM_bcast_d;
static char   _SWARM_bcast_c;
static int    *_SWARM_bcast_ip;
static long   *_SWARM_bcast_lp;
static double *_SWARM_bcast_dp;
static char   *_SWARM_bcast_cp;


#if defined(SOLARIS)&&defined(THREADMAP)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <procfs.h>

int *_solaris_step;

void _solaris_report(THREADED) {
  int s;
  char lwpsinfoFname[MAXLEN];
  int fd_lwpsinfo;
  lwpsinfo_t lwpsinfo;

  sprintf(lwpsinfoFname, "/proc/self/lwp/%d/lwpsinfo", MYTHREAD+2);

  fd_lwpsinfo = open(lwpsinfoFname, O_RDONLY);
  if (fd_lwpsinfo == -1) {
    fprintf(stderr,"ERROR: Trying to open %s\n",lwpsinfoFname);
    perror("ERROR opening proc lwpsinfo");
  }
  s = read(fd_lwpsinfo, &lwpsinfo, sizeof(lwpsinfo_t));
  if (s != sizeof(lwpsinfo_t)) 
    fprintf(stderr,"ERROR: only read %d of %d bytes of lwpsinfo\n",
	    s,sizeof(lwpsinfo_t));
  close(fd_lwpsinfo);

  fprintf(outfile,"THREADMAP %12d\t%3d\t%3d\n",
	  _solaris_step[NOSHARE(MYTHREAD)]++,
	  MYTHREAD,lwpsinfo.pr_onpro);

  return;
}
#endif /* defined(SOLARIS)&&defined(THREADMAP) */

static _SWARM_MULTICORE_barrier_t nbar;

int    SWARM_mutex_init(SWARM_mutex_t **mutex, const SWARM_mutexattr_t *attr, THREADED) {
  int r;
  r = 0;
  *mutex = (SWARM_mutex_t *)SWARM_malloc(sizeof(SWARM_mutex_t), TH);
  on_one_thread {
    r = pthread_mutex_init(*mutex, attr);
  }
  r = SWARM_Bcast_i(r, TH);
  return r;
}

int    SWARM_mutex_destroy(SWARM_mutex_t *mutex, THREADED) {
  int r;
  r = 0;
  SWARM_Barrier();
  on_one_thread {
    r = pthread_mutex_destroy(mutex);
    free (mutex);
  }
  r = SWARM_Bcast_i(r, TH);
  return r;
}

static void SWARM_Barrier_sync_init(void) {
#if defined(SOLARIS)&&defined(THREADMAP)
  int i;
  _solaris_step = (int *)malloc(NOSHARE(THREADS)*sizeof(int));
  assert_malloc(_solaris_step);
  for (i=0 ; i<THREADS ; i++) 
    _solaris_step[NOSHARE(i)] = 0;
#endif /* defined(SOLARIS)&&defined(THREADMAP) */
  nbar = _SWARM_MULTICORE_barrier_init(THREADS);
}

static void SWARM_Barrier_sync_destroy(void) {
  _SWARM_MULTICORE_barrier_destroy(nbar);
}

void SWARM_Barrier_sync(THREADED) {
#if defined(SOLARIS)&&defined(THREADMAP)
  _solaris_report(TH);
#endif /* defined(SOLARIS)&&defined(THREADMAP) */
  _SWARM_MULTICORE_barrier_wait(nbar);
}

static volatile int up_buf[NOSHARE(MAXTHREADS_DEFAULT)][2];
static volatile int down_buf[NOSHARE(MAXTHREADS_DEFAULT)];

static void
SWARM_Barrier_tree_init(void) {
  int i;

  for (i=0 ; i<THREADS ; i++) 
    up_buf[NOSHARE(i)][0] = up_buf[NOSHARE(i)][1] = down_buf[NOSHARE(i)] = 0;
  return;
}

static void
SWARM_Barrier_tree_destroy(void) { return; }

static void
SWARM_Barrier_tree_up(THREADED) {
    
  register int myidx  = MYTHREAD;
  register int parent = (MYTHREAD - 1) / 2;
  register int odd_child = 2 * MYTHREAD + 1;
  register int even_child = 2 * MYTHREAD + 2;
  register int parity = MYTHREAD & 1;

  if (MYTHREAD == 0) {
    if (THREADS > 2) {
      while (up_buf[NOSHARE(myidx)][0] == 0 ||
	     up_buf[NOSHARE(myidx)][1] == 0) ;
    }
    else if (THREADS == 2) {
	while (up_buf[NOSHARE(myidx)][1] == 0) ;
    }
  } 
  else 
    if (odd_child >= THREADS) 
      up_buf[NOSHARE(parent)][parity]++;
    else 
      if (even_child >= THREADS) {
	while (up_buf[NOSHARE(myidx)][1] == 0) ;
	up_buf[NOSHARE(parent)][parity]++;
      } 
      else {
	while (up_buf[NOSHARE(myidx)][0] == 0 ||
	       up_buf[NOSHARE(myidx)][1] == 0) ;
	up_buf[NOSHARE(parent)][parity]++;
      }

  up_buf[NOSHARE(myidx)][0] = up_buf[NOSHARE(myidx)][1] = 0;
#ifdef SOLARIS
  sun_mb_mi_np();
#endif
  return;
}

static void
SWARM_Barrier_tree_down(THREADED) {
    
  register int myidx  = MYTHREAD;
  register int left = 2 * MYTHREAD + 1;
  register int right = 2 * MYTHREAD + 2;

  if (MYTHREAD != 0) 
    while (down_buf[NOSHARE(myidx)] == 0) ;
  
  if (left < THREADS)
    down_buf[NOSHARE(left)]++;
  if (right < THREADS)
    down_buf[NOSHARE(right)]++;

  down_buf[NOSHARE(myidx)] = 0;
#ifdef SOLARIS
  sun_mb_mi_np();
#endif
  return;
}

void
SWARM_Barrier_tree(THREADED) {
  SWARM_Barrier_tree_up(TH);
  SWARM_Barrier_tree_down(TH);
  return;
}

#ifdef SUNMMAP
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
int _fdmmap = 0;
int _offmmap = 0;
int _ps = 0;
#endif

void   *SWARM_malloc(int bytes, THREADED) {
  void *ptr;
  ptr=NULL;
  on_one_thread {
#ifdef SUNMMAP
    if (!_fdmmap) {
      _fdmmap = open("/dev/zero", O_RDWR);
      if (_fdmmap < 0)
	perror("fdmmap error");
      _ps = sysconf(_SC_PAGESIZE);
    }
    if (bytes % _ps)
      bytes = ((bytes / _ps) + 1) * _ps;
    ptr = mmap((void *)0, bytes,
	       (PROT_READ | PROT_WRITE), MAP_PRIVATE, _fdmmap, _offmmap);
    _offmmap += bytes;
    if (ptr == MAP_FAILED)
      perror("mmap failed");
#else
    ptr = malloc(bytes);
    assert_malloc(ptr);
#endif
  }
  return(SWARM_Bcast_cp(ptr, TH));
}

void   *SWARM_malloc_l(long bytes, THREADED) {
  void *ptr;
  ptr=NULL;
  on_one_thread {
#ifdef SUNMMAP
    if (!_fdmmap) {
      _fdmmap = open("/dev/zero", O_RDWR);
      if (_fdmmap < 0)
	perror("fdmmap error");
      _ps = sysconf(_SC_PAGESIZE);
    }
    if (bytes % _ps)
      bytes = ((bytes / _ps) + 1) * _ps;
    ptr = mmap((void *)0, bytes,
	       (PROT_READ | PROT_WRITE), MAP_PRIVATE, _fdmmap, _offmmap);
    _offmmap += bytes;
    if (ptr == MAP_FAILED)
      perror("mmap failed");
#else

#ifdef SOLARIS
    /*  ptr = malloc64(bytes); */
    ptr = malloc(bytes);
#else
    ptr = malloc(bytes);
#endif
    
    assert_malloc(ptr);
#endif
  }
  return(SWARM_Bcast_cp(ptr, TH));
}

void   SWARM_free(void *ptr, THREADED) {
  on_one_thread {
#ifdef SUNMMAP
    ptr = (void *)NULL;
#else
    free(ptr);
#endif
  }
}

int    SWARM_Bcast_i(int    myval, THREADED) {

  SWARM_Barrier();

  on_one_thread {
    _SWARM_bcast_i = myval;
  }

  SWARM_Barrier();
  return (_SWARM_bcast_i);
}

long   SWARM_Bcast_l(long   myval, THREADED) {

  SWARM_Barrier();

  on_one_thread {
    _SWARM_bcast_l = myval;
  }

  SWARM_Barrier();
  return (_SWARM_bcast_l);
}

double SWARM_Bcast_d(double myval, THREADED) {

  SWARM_Barrier();

  on_one_thread {
    _SWARM_bcast_d = myval;
  }

  SWARM_Barrier();
  return (_SWARM_bcast_d);
}

char   SWARM_Bcast_c(char   myval, THREADED) {

  SWARM_Barrier();

  on_one_thread {
    _SWARM_bcast_c = myval;
  }

  SWARM_Barrier();
  return (_SWARM_bcast_c);
}

int    *SWARM_Bcast_ip(int    *myval, THREADED) {

  SWARM_Barrier();

  on_one_thread {
    _SWARM_bcast_ip = myval;
  }

  SWARM_Barrier();
  return (_SWARM_bcast_ip);
}

long   *SWARM_Bcast_lp(long   *myval, THREADED) {

  SWARM_Barrier();

  on_one_thread {
    _SWARM_bcast_lp = myval;
  }

  SWARM_Barrier();
  return (_SWARM_bcast_lp);
}

double *SWARM_Bcast_dp(double *myval, THREADED) {

  SWARM_Barrier();

  on_one_thread {
    _SWARM_bcast_dp = myval;
  }

  SWARM_Barrier();
  return (_SWARM_bcast_dp);
}

char   *SWARM_Bcast_cp(char   *myval, THREADED) {

  SWARM_Barrier();

  on_one_thread {
    _SWARM_bcast_cp = myval;
  }

  SWARM_Barrier();
  return (_SWARM_bcast_cp);
}

static _SWARM_MULTICORE_reduce_d_t red_d;

double SWARM_Reduce_d(double myval, reduce_t op, THREADED) {
  return (_SWARM_MULTICORE_reduce_d(red_d, myval, op));
}


static _SWARM_MULTICORE_reduce_i_t red_i;

int SWARM_Reduce_i(int myval, reduce_t op, THREADED) {
  return (_SWARM_MULTICORE_reduce_i(red_i,myval,op));
}

static _SWARM_MULTICORE_reduce_l_t red_l;

long SWARM_Reduce_l(long myval, reduce_t op, THREADED) {
  return (_SWARM_MULTICORE_reduce_l(red_l,myval,op));
}

static _SWARM_MULTICORE_scan_i_t scan_i;

int SWARM_Scan_i(int myval, reduce_t op, THREADED) {
  return(_SWARM_MULTICORE_scan_i(scan_i,myval,op,MYTHREAD));
}

static _SWARM_MULTICORE_scan_l_t scan_l;

long SWARM_Scan_l(long myval, reduce_t op, THREADED) {
  return(_SWARM_MULTICORE_scan_l(scan_l,myval,op,MYTHREAD));
}

static _SWARM_MULTICORE_scan_d_t scan_d;

double SWARM_Scan_d(double myval, reduce_t op, THREADED) {
  return(_SWARM_MULTICORE_scan_d(scan_d,myval,op,MYTHREAD));
}

static void SWARM_print_help(char **argv)
{
     printf ("SWARM usage: %s [-t #threads] [-o outfile]\n", argv[0]);
     printf ("\t-t #threads    overrides the default number of threads\n");
     printf ("\t-o outfile     redirects standard output to outfile\n");
}

#define ERRSTRINGSIZE 512

static void SWARM_error (int lineno, const char *func, const char *format, ...)
{
    char msg[ERRSTRINGSIZE];

    va_list args;
    va_start(args, format);
    vsprintf(msg, format, args);

    fprintf (stderr, "SWARM_%s (line %d): %s\n", func, lineno, msg);
    fflush (stderr);

    SWARM_Finalize();
    
    exit (-1);

}


FILE *SWARM_outfile;
static char *SWARM_outfilename;


/*
 *
 * Parses command line arguments passed to it in argc and argv, 
 * using getopt(3c). Sets the following flags when it sees the
 * corresponding command line arguments:
 *
 *   -t  (overrides the default number of threads)
 *   -o  (redirect stdout)
 *   -h  (help)
 *
 * Returns 0 if there are no more arguments, or the index of the next
 * non-option argument.
 */
static int SWARM_read_arguments (int argc, char **argv)
{
     extern char *optarg;
     extern int optind;
     char *tail;
     int c, i;

     if (argv[0] == NULL) 
	  SWARM_error (__LINE__, "SWARM_read_arguments", 
		       "invalid argument array"); 
     
     if (argc < 1) return 0;
     
     opterr = 0;
     while ((c = getopt (argc, argv, "ht:o:")) != -1)
     {
	  switch (c)
	  {
	       
	  case 't':
	       i = (int)strtol (optarg, &tail, 0);
	       if (optarg == tail)
		    SWARM_error (__LINE__, "read_arguments", 
				 "invalid argument %s to option -t", optarg); 
	       if (i <= 0)
		    SWARM_error (__LINE__, "read_arguments", 
				 "# of threads must be greater than zero");
	       else
		    THREADS = i;
	       break;

	  case 'o':
	       SWARM_outfilename = strdup(optarg);
	       if ((SWARM_outfile = fopen (SWARM_outfilename, "w")) == NULL)
		    SWARM_error (__LINE__, "read_arguments", 
				 "unable to open outfile %s", SWARM_outfilename);
	       break;
	       
	  case 'h':
	       SWARM_print_help(argv);
	       exit(0);
	       break;

	  default:
		SWARM_error (__LINE__, "read_arguments", 
			     "`%c': unrecognized argument", c);

                break;
	   }
      }
    

     if (argv[optind] != NULL) return optind;
     else return 0;
}



static void
SWARM_get_args(int *argc, char* **argv) {
  int numarg = *argc;
  int done = 0;
  char
    *s,**argvv = *argv;
  char
    *outfilename = NULL;

  SWARM_outfile = stdout;

  while ((--numarg > 0) && !done)
    if ((*++argvv)[0] == '-')
      for (s=argvv[0]+1; *s != '\0'; s++) {
	if (*s == '-')
	  done = 1;
	else {
	  switch (*s) {
	  case 'o':
	    if (numarg <= 1) 
	      perror("output filename expected after -o (e.g. -o filename)");
	    numarg--;
	    outfilename = (char *)malloc(MAXLEN*sizeof(char));
	    strcpy(outfilename, *++argvv); 
	    SWARM_outfile = fopen(outfilename,"a+");
	    break;
	  case 't':
	    if (numarg <= 1) 
	      perror("number of threads per node expected after -t");
	    numarg--;
	    THREADS = atoi(*++argvv);

	    break;
	  case 'h':
	    fprintf(SWARM_outfile,"SWARM Options:\n");
	    fprintf(SWARM_outfile," -t <number of threads per node>\n");
	    fprintf(SWARM_outfile,"\n\n");
	    fflush(SWARM_outfile);
	    break;
	    /*	default: perror("illegal option");  */
	  }
	}
      }
  if (done) {
    *argc = numarg;
    *argv = ++argvv;
  }
  else {
    *argc = 0;
    *argv = NULL;
  }

  return;
}

#ifdef WIN32

// Avoids including windows.h
// Used for getting number of cores information on Windows machine
typedef struct _SYSTEM_INFO 
{  
	union 
	{    
		unsigned long dwOemId;    
		struct 
		{      
			unsigned short wProcessorArchitecture;      
			unsigned short wReserved;    
		};  
	};  
	unsigned long dwPageSize;  
	unsigned long lpMinimumApplicationAddress;  
	unsigned long lpMaximumApplicationAddress;  
	unsigned long* dwActiveProcessorMask;  
	unsigned long dwNumberOfProcessors;  
	unsigned long dwProcessorType;  
	unsigned long dwAllocationGranularity;  
	unsigned short wProcessorLevel;  
	unsigned short wProcessorRevision;
}SYSTEM_INFO;

void __stdcall GetSystemInfo(SYSTEM_INFO*);

#endif

static int SWARM_get_num_cores(void)
{
	int num_cores = DEFAULT_THREADS;
	
	#ifdef WIN32
		SYSTEM_INFO siSysInfo;
	#endif

	#ifdef WIN32 
 		GetSystemInfo(&siSysInfo); 
 		num_cores = siSysInfo.dwNumberOfProcessors;
	#elif HAVE_SYSCONF
		#ifdef _SC_NPROCESSORS_ONLN
    		num_cores = (int)sysconf(_SC_NPROCESSORS_ONLN);
		#endif
	#endif

	return num_cores;
}
	

#ifdef HAVE_PTHREAD_SCHED_SUPPORTED
static pthread_attr_t pattr;
#endif

void SWARM_Init(int *argc, char* **argv)
{
	
	int i;
	#ifdef HAVE_PTHREAD_SCHED_SUPPORTED
    int rc;
	#endif
    uthread_info_t *ti;
    int moreargs;

	#ifdef HAVE_PTHREAD_SCHED_SUPPORTED
    struct sched_param psched;
	#endif

    THREADS = SWARM_get_num_cores();

    SWARM_outfile  = stdout;
    SWARM_outfilename = NULL;
    
    moreargs = SWARM_read_arguments (*argc, *argv);

#if INFO
    fprintf(SWARM_outfile,"THREADS: %d\n", THREADS);
    fflush(SWARM_outfile);
#endif /*INFO */

    /*********************************/
    /* ON ONE THREAD INITIALIZATION  */
    /*********************************/

    SWARM_Barrier_sync_init();
    SWARM_Barrier_tree_init();

    red_i = _SWARM_MULTICORE_reduce_init_i(THREADS);
    red_l = _SWARM_MULTICORE_reduce_init_l(THREADS);
    red_d = _SWARM_MULTICORE_reduce_init_d(THREADS);

    scan_i = _SWARM_MULTICORE_scan_init_i(THREADS);
    scan_l = _SWARM_MULTICORE_scan_init_l(THREADS);
    scan_d = _SWARM_MULTICORE_scan_init_d(THREADS);

#ifdef HAVE_PTHREAD_SCHED_SUPPORTED

    /*********************/
    /* THREAD  SCHEDULER */
    /*********************/
    rc = pthread_attr_init(&pattr);
    if (rc)
      perror("pthread_attr_init");

    rc = pthread_attr_setschedpolicy(&pattr, SCHED_FIFO);
    if (rc)
      perror("pthread_attr_setschedpolicy");
    
    psched.sched_priority = sched_get_priority_max(SCHED_FIFO);
    
    rc = pthread_attr_setschedparam(&pattr, &psched);
    if (rc)
      perror("pthread_attr_isetschedparam");
    
    rc = pthread_attr_setinheritsched(&pattr, PTHREAD_EXPLICIT_SCHED);
    if (rc)
      perror("pthread_attr_setinheritsched");
    
#endif /* HAVE_PTHREAD_SCHED_SUPPORTED */

#if (defined(SOLARIS))
    rc = pthread_setconcurrency(THREADS+1);
    if (rc)
      perror("pthread_setconcurrency");
#endif /* defined(SOLARIS) */
    
    spawn_thread = (pthread_t *)malloc(NOSHARE(THREADS)*
				       sizeof(pthread_t));
    assert_malloc(spawn_thread);
    uthread_info = (uthread_info_t *)malloc(NOSHARE(THREADS)*
					  sizeof(uthread_info_t));
    assert_malloc(uthread_info);

    ti = uthread_info;

    for (i=0 ; i<THREADS ; i++) {
      ti->mythread   = i;

      if (moreargs > 0) 
      {
	   ti->argc       = (*argc)-moreargs;
	   ti->argv       = (*argv)+moreargs;
      }
      else 
      {
	   ti->argc       = 0;
	   ti->argv       = (char **)NULL;
      }

#ifdef RRANDOM
      SWARM_rrandom_init(ti);
#endif /* RRANDOM */

#ifdef HAVE_LIBSPRNG
      THSPRNG        = (int *)NULL;
#endif
      ti += NOSHARE(1);
    }

    _swarm_init=1;
}

void SWARM_Run(void *start_routine) 
{
     int i, rc;
     int *parg;
     uthread_info_t *ti;
     void *(*f)(void *);
     
     f = (void *(*)(void *))start_routine;
     
     if (!_swarm_init) 
     {
	  fprintf(stderr,"ERROR: SWARM_Init() not called\n");
	  perror("SWARM_Run cannot call SWARM_Init()");
     }
     
     ti = uthread_info;
     
     for (i=0 ; i<THREADS ; i++) 
     {
	  
	  rc = pthread_create(spawn_thread + NOSHARE(i),
#ifdef HAVE_PTHREAD_SCHED_SUPPORTED
			      &pattr,
#else
			      NULL,
#endif
			      f,
			      ti);

	  if (rc != 0)
	  {
	       switch (rc)
	       {
	       case EAGAIN: 
		    SWARM_error (__LINE__, "Run:pthread_create", 
				 "not enough resources to create another thread");
		    break;

	       case EINVAL:
		    SWARM_error (__LINE__, "Run:pthread_create", 
				 "invalid thread attributes");
		    break;

	       case EPERM:
		    SWARM_error (__LINE__, "Run:pthread_create", 
				 "insufficient permissions for setting scheduling parameters or policy ");
		    break;

	       default:
		    SWARM_error (__LINE__, "Run:pthread_create", "error code %d", rc);
	       }
	  }
			 
	  ti += NOSHARE(1);
     }
     
     for (i=0 ; i<THREADS ; i++)
     {
	  rc = pthread_join(spawn_thread[NOSHARE(i)], (void *)&parg);
	  if (rc != 0)
	  {
	       switch (rc)
	       {
	       case EINVAL:
		    SWARM_error (__LINE__, "Run:pthread_join", "specified thread is not joinable");
		    break;

	       case ESRCH:
		    SWARM_error (__LINE__, "Run:pthread_join", "cannot find thread");
		    break;

	       default:
		    SWARM_error (__LINE__, "Run:pthread_join", "error code %d", rc);

	       }
	  }
     }
}

void SWARM_Finalize(void) 
{
     
     /*********************************/
     /* ONE ONE THREAD DESTRUCTION    */
     /*********************************/
     
     _SWARM_MULTICORE_reduce_destroy_i(red_i);
     _SWARM_MULTICORE_reduce_destroy_l(red_l);
     _SWARM_MULTICORE_reduce_destroy_d(red_d);
     _SWARM_MULTICORE_scan_destroy_i(scan_i);
     _SWARM_MULTICORE_scan_destroy_l(scan_l);
     _SWARM_MULTICORE_scan_destroy_d(scan_d);
     
     SWARM_Barrier_sync_destroy();
     SWARM_Barrier_tree_destroy();
     
     free(uthread_info);
     free(spawn_thread);
     
     if (SWARM_outfile != NULL)
	  fclose(SWARM_outfile);
     if (SWARM_outfilename != NULL)
	  free(SWARM_outfilename);
}

void SWARM_Cleanup(THREADED) 
{
#ifdef RRANDOM
     SWARM_rrandom_destroy(TH);
#endif /* RRANDOM */
     return;
}

void assert_malloc(void *ptr) 
{
	if (ptr==NULL)
    	perror("ERROR: assert_malloc");
}

double get_seconds(void)
{
	struct timeval t;
  	struct timezone z;
  	gettimeofday(&t,&z);
  	return (double)t.tv_sec+((double)t.tv_usec/(double)1000000.0);
}


#ifndef HAVE_GETTIMEOFDAY
	/*  Windows does not support the gettimeofday. This fix from Anders
 	*  Carlsson was lifted from:
 	*  http://lists.gnu.org/archive/html/bug-gnu-chess/2004-01/msg00020.html
	*/

 	/* These are winbase.h definitions, but to avoid including
    	tons of Windows related stuff, it is reprinted here */

	typedef struct _FILETIME
	{
    	unsigned long dwLowDateTime;
    	unsigned long dwHighDateTime;
	} FILETIME;

	void __stdcall GetSystemTimeAsFileTime(FILETIME*);

	void gettimeofday(struct timeval* p, struct timezone* tz /* IGNORED */)
	{
   		union 
		{
   			long long ns100; /*time since 1 Jan 1601 in 100ns units */
   			FILETIME ft;
		} _now;

		GetSystemTimeAsFileTime( &(_now.ft) );
   		p->tv_usec=(long)((_now.ns100 / 10LL) % 1000000LL );
   		p->tv_sec= (long)((_now.ns100-(116444736000000000LL))/10000000LL);
   		return;
	}
#endif

#ifndef HAVE_GETOPT_H

	int	opterr = 1;
	int optind = 1;
	int optopt = 0;
	char *optarg = 0;

	int getopt(int argc, char **argv, char *opts)
	{
		static int sp = 1;
		register int c;
		register char *cp;

		if(sp == 1) 
		{
			if(optind >= argc || (argv[optind][0] != '+' &&
             argv[optind][0] != '-') || argv[optind][1] == '\0')
               return EOF;
			else if(strcmp(argv[optind], "--") == 0) 
			{
               optind++;
               return EOF;
			}
			/* '+' for config options, '+' should not be in the opts list */
			if (argv[optind][0] == '+') 
			{
               optarg = argv[optind++] + 1;
               return '+';
			}
		}
		optopt = c = argv[optind][sp];
		if(c == ':' || (cp=strchr(opts, c)) == NULL) 
		{
			ERR(": illegal option -- ", c);
			if(argv[optind][++sp] == '\0') 
			{
               optind++;
               sp = 1;
			}
			return '\0';
		}
		if(*++cp == ':') 
		{
			if(argv[optind][sp+1] != '\0')
               optarg = &argv[optind++][sp+1];
			else if(++optind >= argc) 
			{
               ERR(": option requires an argument -- ", c);
               sp = 1;
               return '\0';
			} 
			else
               optarg = argv[optind++];
			sp = 1;
		} 
		else 
		{
			if(argv[optind][++sp] == '\0') 
			{
               sp = 1;
               optind++;
			}
			optarg = NULL;
		}
		return c;
	}
#endif

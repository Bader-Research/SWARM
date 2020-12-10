#ifndef _SWARM_MULTICORE_H
#define _SWARM_MULTICORE_H

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <swarm.h>

typedef struct _SWARM_MULTICORE_barrier {
  pthread_mutex_t lock;
  int n_clients;
  int n_waiting;
  int phase;
  pthread_cond_t wait_cv;
} *_SWARM_MULTICORE_barrier_t;

_SWARM_MULTICORE_barrier_t _SWARM_MULTICORE_barrier_init(int n_clients);
void          _SWARM_MULTICORE_barrier_destroy(_SWARM_MULTICORE_barrier_t nbarrier);
void          _SWARM_MULTICORE_barrier_wait(_SWARM_MULTICORE_barrier_t nbarrier);

typedef struct _SWARM_MULTICORE_reduce_i_s {
  pthread_mutex_t lock;
  int n_clients;
  int n_waiting;
  int phase;
  int sum;
  int result;
  pthread_cond_t wait_cv;
} *_SWARM_MULTICORE_reduce_i_t;

_SWARM_MULTICORE_reduce_i_t _SWARM_MULTICORE_reduce_init_i(int n_clients);
void           _SWARM_MULTICORE_reduce_destroy_i(_SWARM_MULTICORE_reduce_i_t nbarrier);
int            _SWARM_MULTICORE_reduce_i(_SWARM_MULTICORE_reduce_i_t nbarrier, int val, reduce_t op);

typedef struct _SWARM_MULTICORE_reduce_l_s {
  pthread_mutex_t lock;
  int n_clients;
  int n_waiting;
  int phase;
  long sum;
  long result;
  pthread_cond_t wait_cv;
} *_SWARM_MULTICORE_reduce_l_t;

_SWARM_MULTICORE_reduce_l_t _SWARM_MULTICORE_reduce_init_l(int n_clients);
void           _SWARM_MULTICORE_reduce_destroy_l(_SWARM_MULTICORE_reduce_l_t nbarrier);
long           _SWARM_MULTICORE_reduce_l(_SWARM_MULTICORE_reduce_l_t nbarrier, long val, reduce_t op);

typedef struct _SWARM_MULTICORE_reduce_d_s {
  pthread_mutex_t lock;
  int n_clients;
  int n_waiting;
  int phase;
  double sum;
  double result;
  pthread_cond_t wait_cv;
} *_SWARM_MULTICORE_reduce_d_t;

_SWARM_MULTICORE_reduce_d_t _SWARM_MULTICORE_reduce_init_d(int n_clients);
void           _SWARM_MULTICORE_reduce_destroy_d(_SWARM_MULTICORE_reduce_d_t nbarrier);
double         _SWARM_MULTICORE_reduce_d(_SWARM_MULTICORE_reduce_d_t nbarrier, double val, reduce_t op);

typedef struct _SWARM_MULTICORE_scan_i_s {
  pthread_mutex_t lock;
  int n_clients;
  int n_waiting;
  int phase;
  int *result;
  pthread_cond_t wait_cv;
} *_SWARM_MULTICORE_scan_i_t;

_SWARM_MULTICORE_scan_i_t   _SWARM_MULTICORE_scan_init_i(int n_clients);
void           _SWARM_MULTICORE_scan_destroy_i(_SWARM_MULTICORE_scan_i_t nbarrier);
int            _SWARM_MULTICORE_scan_i(_SWARM_MULTICORE_scan_i_t nbarrier, int val, reduce_t op,int th_index);


typedef struct _SWARM_MULTICORE_scan_l_s {
  pthread_mutex_t lock;
  int n_clients;
  int n_waiting;
  int phase;
  long *result;
  pthread_cond_t wait_cv;
} *_SWARM_MULTICORE_scan_l_t;

_SWARM_MULTICORE_scan_l_t   _SWARM_MULTICORE_scan_init_l(int n_clients);
void           _SWARM_MULTICORE_scan_destroy_l(_SWARM_MULTICORE_scan_l_t nbarrier);
long           _SWARM_MULTICORE_scan_l(_SWARM_MULTICORE_scan_l_t nbarrier, long val, reduce_t op,int th_index);

typedef struct _SWARM_MULTICORE_scan_d_s {
  pthread_mutex_t lock;
  int n_clients;
  int n_waiting;
  int phase;
  double *result;
  pthread_cond_t wait_cv;
} *_SWARM_MULTICORE_scan_d_t;

_SWARM_MULTICORE_scan_d_t   _SWARM_MULTICORE_scan_init_d(int n_clients);
void           _SWARM_MULTICORE_scan_destroy_d(_SWARM_MULTICORE_scan_d_t nbarrier);
double         _SWARM_MULTICORE_scan_d(_SWARM_MULTICORE_scan_d_t nbarrier, double val, reduce_t op,int th_index);

typedef struct _SWARM_MULTICORE_spin_barrier {
  int n_clients;
  pthread_mutex_t lock;
  int n_waiting;
  int phase;
} *_SWARM_MULTICORE_spin_barrier_t;

_SWARM_MULTICORE_spin_barrier_t _SWARM_MULTICORE_spin_barrier_init(int n_clients);
void           _SWARM_MULTICORE_spin_barrier_destroy(_SWARM_MULTICORE_spin_barrier_t sbarrier);
void           _SWARM_MULTICORE_spin_barrier_wait(_SWARM_MULTICORE_spin_barrier_t sbarrier);

#endif

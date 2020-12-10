#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <errno.h>
#include <swarm.h>
#include <swarm_multicore.h>

_SWARM_MULTICORE_barrier_t _SWARM_MULTICORE_barrier_init(int n_clients) {
  _SWARM_MULTICORE_barrier_t nbarrier = (_SWARM_MULTICORE_barrier_t)
    malloc(sizeof(struct _SWARM_MULTICORE_barrier));
  assert_malloc(nbarrier);
  
  if (nbarrier != NULL) {
    nbarrier->n_clients = n_clients;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 0;
    pthread_mutex_init(&nbarrier->lock, NULL);
    pthread_cond_init(&nbarrier->wait_cv, NULL);
  }
  return(nbarrier);
}

void _SWARM_MULTICORE_barrier_destroy(_SWARM_MULTICORE_barrier_t nbarrier) {
  pthread_mutex_destroy(&nbarrier->lock);
  pthread_cond_destroy(&nbarrier->wait_cv);
  free(nbarrier);
}

void _SWARM_MULTICORE_barrier_wait(_SWARM_MULTICORE_barrier_t nbarrier) {
  int my_phase;

  pthread_mutex_lock(&nbarrier->lock);
  my_phase = nbarrier->phase;
  nbarrier->n_waiting++;
  if (nbarrier->n_waiting == nbarrier->n_clients) {
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 1 - my_phase;
    pthread_cond_broadcast(&nbarrier->wait_cv);
  }
  while (nbarrier->phase == my_phase) {
    pthread_cond_wait(&nbarrier->wait_cv, &nbarrier->lock);
  }
  pthread_mutex_unlock(&nbarrier->lock);
}

_SWARM_MULTICORE_reduce_i_t _SWARM_MULTICORE_reduce_init_i(int n_clients) {
  _SWARM_MULTICORE_reduce_i_t nbarrier = (_SWARM_MULTICORE_reduce_i_t)
    malloc(sizeof(struct _SWARM_MULTICORE_reduce_i_s));
  assert_malloc(nbarrier);
  
  if (nbarrier != NULL) {
    nbarrier->n_clients = n_clients;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 0;
    nbarrier->sum       = 0;
    pthread_mutex_init(&nbarrier->lock, NULL);
    pthread_cond_init(&nbarrier->wait_cv, NULL);
  }
  return(nbarrier);
}

void _SWARM_MULTICORE_reduce_destroy_i(_SWARM_MULTICORE_reduce_i_t nbarrier) {
  pthread_mutex_destroy(&nbarrier->lock);
  pthread_cond_destroy(&nbarrier->wait_cv);
  free(nbarrier);
}

int _SWARM_MULTICORE_reduce_i(_SWARM_MULTICORE_reduce_i_t nbarrier, int val, reduce_t op) {
  int my_phase;

  pthread_mutex_lock(&nbarrier->lock);
  my_phase = nbarrier->phase;
  if (nbarrier->n_waiting==0) {
    nbarrier->sum = val;
  }
  else {
    switch (op) {
    case MIN:  nbarrier->sum  = min(nbarrier->sum,val);  break;
    case MAX : nbarrier->sum  = max(nbarrier->sum,val);  break;
    case SUM : nbarrier->sum += val;  break;
    default  : perror("ERROR: _SWARM_MULTICORE_reduce_i() Bad reduction operator");
    }
  }
  nbarrier->n_waiting++;
  if (nbarrier->n_waiting == nbarrier->n_clients) {
    nbarrier->result    = nbarrier->sum;
    nbarrier->sum       = 0;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 1 - my_phase;
    pthread_cond_broadcast(&nbarrier->wait_cv);
  }
  while (nbarrier->phase == my_phase) {
    pthread_cond_wait(&nbarrier->wait_cv, &nbarrier->lock);
  }
  pthread_mutex_unlock(&nbarrier->lock);
  return(nbarrier->result);
}

_SWARM_MULTICORE_reduce_l_t _SWARM_MULTICORE_reduce_init_l(int n_clients) {
  _SWARM_MULTICORE_reduce_l_t nbarrier = (_SWARM_MULTICORE_reduce_l_t)
    malloc(sizeof(struct _SWARM_MULTICORE_reduce_l_s));
  assert_malloc(nbarrier);

  if (nbarrier != NULL) {
    nbarrier->n_clients = n_clients;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 0;
    nbarrier->sum       = 0;
    pthread_mutex_init(&nbarrier->lock, NULL);
    pthread_cond_init(&nbarrier->wait_cv, NULL);
  }
  return(nbarrier);
}

void _SWARM_MULTICORE_reduce_destroy_l(_SWARM_MULTICORE_reduce_l_t nbarrier) {
  pthread_mutex_destroy(&nbarrier->lock);
  pthread_cond_destroy(&nbarrier->wait_cv);
  free(nbarrier);
}

long _SWARM_MULTICORE_reduce_l(_SWARM_MULTICORE_reduce_l_t nbarrier, long val, reduce_t op) {
  int my_phase;

  pthread_mutex_lock(&nbarrier->lock);
  my_phase = nbarrier->phase;
  if (nbarrier->n_waiting==0) {
    nbarrier->sum = val;
  }
  else {
    switch (op) {
    case MIN:  nbarrier->sum  = min(nbarrier->sum,val);  break;
    case MAX : nbarrier->sum  = max(nbarrier->sum,val);  break;
    case SUM : nbarrier->sum += val;  break;
    default  : perror("ERROR: _SWARM_MULTICORE_reduce_l() Bad reduction operator");
    }
  }
  nbarrier->n_waiting++;
  if (nbarrier->n_waiting == nbarrier->n_clients) {
    nbarrier->result    = nbarrier->sum;
    nbarrier->sum       = 0;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 1 - my_phase;
    pthread_cond_broadcast(&nbarrier->wait_cv);
  }
  while (nbarrier->phase == my_phase) {
    pthread_cond_wait(&nbarrier->wait_cv, &nbarrier->lock);
  }
  pthread_mutex_unlock(&nbarrier->lock);
  return(nbarrier->result);
}

_SWARM_MULTICORE_reduce_d_t _SWARM_MULTICORE_reduce_init_d(int n_clients) {
  _SWARM_MULTICORE_reduce_d_t nbarrier = (_SWARM_MULTICORE_reduce_d_t)
    malloc(sizeof(struct _SWARM_MULTICORE_reduce_d_s));
  assert_malloc(nbarrier);

  if (nbarrier != NULL) {
    nbarrier->n_clients = n_clients;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 0;
    nbarrier->sum       = 0.000001;
    pthread_mutex_init(&nbarrier->lock, NULL);
    pthread_cond_init(&nbarrier->wait_cv, NULL);
  }
  return(nbarrier);
}

void _SWARM_MULTICORE_reduce_destroy_d(_SWARM_MULTICORE_reduce_d_t nbarrier) {
  pthread_mutex_destroy(&nbarrier->lock);
  pthread_cond_destroy(&nbarrier->wait_cv);
  free(nbarrier);
}

double _SWARM_MULTICORE_reduce_d(_SWARM_MULTICORE_reduce_d_t nbarrier, double val, reduce_t op) {
  int my_phase;

  pthread_mutex_lock(&nbarrier->lock);
  my_phase = nbarrier->phase;
  if (nbarrier->n_waiting==0) {
    nbarrier->sum = val;
  }
  else {
    switch (op) {
    case MIN:  nbarrier->sum  = min(nbarrier->sum,val);  break;
    case MAX : nbarrier->sum  = max(nbarrier->sum,val);  break;
    case SUM : nbarrier->sum += val;  break;
    default  : perror("ERROR: _SWARM_MULTICORE_reduce_i() Bad reduction operator");
    }
  }
  nbarrier->n_waiting++;
  if (nbarrier->n_waiting == nbarrier->n_clients) {
    nbarrier->result    = nbarrier->sum;
    nbarrier->sum       = 0.0;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 1 - my_phase;
    pthread_cond_broadcast(&nbarrier->wait_cv);
  }
  while (nbarrier->phase == my_phase) {
    pthread_cond_wait(&nbarrier->wait_cv, &nbarrier->lock);
  }
  pthread_mutex_unlock(&nbarrier->lock);
  return(nbarrier->result);
}

_SWARM_MULTICORE_scan_i_t _SWARM_MULTICORE_scan_init_i(int n_clients) {
  _SWARM_MULTICORE_scan_i_t nbarrier = (_SWARM_MULTICORE_scan_i_t)
    malloc(sizeof(struct _SWARM_MULTICORE_scan_i_s));
  assert_malloc(nbarrier);

  if (nbarrier != NULL) {
    nbarrier->n_clients = n_clients;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 0;
    nbarrier->result    = (int *)malloc(n_clients*sizeof(int));
    assert_malloc(nbarrier->result);
    pthread_mutex_init(&nbarrier->lock, NULL);
    pthread_cond_init(&nbarrier->wait_cv, NULL);
  }
  return(nbarrier);
}

void _SWARM_MULTICORE_scan_destroy_i(_SWARM_MULTICORE_scan_i_t nbarrier) {
  pthread_mutex_destroy(&nbarrier->lock);
  pthread_cond_destroy(&nbarrier->wait_cv);
  free(nbarrier->result);
  free(nbarrier);
}

int _SWARM_MULTICORE_scan_i(_SWARM_MULTICORE_scan_i_t nbarrier, int val, reduce_t op,int th_index) {
  int my_phase,i,temp;

  pthread_mutex_lock(&nbarrier->lock);
  my_phase = nbarrier->phase;
  nbarrier->result[th_index]  = val;

  nbarrier->n_waiting++;
  if (nbarrier->n_waiting == nbarrier->n_clients) { /* get the prefix result in result array*/
    switch (op) {
    case MIN : temp = nbarrier->result[0];
      for(i = 1; i < nbarrier->n_clients;i++) {
         temp  = min(nbarrier->result[i],temp);
         nbarrier->result[i] = temp;
      }  
      break;
    case MAX : temp = nbarrier->result[0];
      for(i = 1; i < nbarrier->n_clients;i++) {
         temp  = max(nbarrier->result[i],temp);
         nbarrier->result[i] = temp;
      }  
      break;
    case SUM : 
      for(i = 1; i < nbarrier->n_clients;i++) 
         nbarrier->result[i] += nbarrier->result[i-1];
      break;
    default  : perror("ERROR: _SWARM_MULTICORE_scan_i() Bad reduction operator");
    }
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 1 - my_phase;
    pthread_cond_broadcast(&nbarrier->wait_cv);
  }
  while (nbarrier->phase == my_phase) {
    pthread_cond_wait(&nbarrier->wait_cv, &nbarrier->lock);
  }
  pthread_mutex_unlock(&nbarrier->lock);
  return(nbarrier->result[th_index]);
}

_SWARM_MULTICORE_scan_l_t _SWARM_MULTICORE_scan_init_l(int n_clients) {
  _SWARM_MULTICORE_scan_l_t nbarrier = (_SWARM_MULTICORE_scan_l_t)
    malloc(sizeof(struct _SWARM_MULTICORE_scan_l_s));
  assert_malloc(nbarrier);

  if (nbarrier != NULL) {
    nbarrier->n_clients = n_clients;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 0;
    nbarrier->result    = (long *)malloc(n_clients*sizeof(long));
    assert_malloc(nbarrier->result);
    pthread_mutex_init(&nbarrier->lock, NULL);
    pthread_cond_init(&nbarrier->wait_cv, NULL);
  }
  return(nbarrier);
}

void _SWARM_MULTICORE_scan_destroy_l(_SWARM_MULTICORE_scan_l_t nbarrier) {
  pthread_mutex_destroy(&nbarrier->lock);
  pthread_cond_destroy(&nbarrier->wait_cv);
  free(nbarrier->result);
  free(nbarrier);
}

long _SWARM_MULTICORE_scan_l(_SWARM_MULTICORE_scan_l_t nbarrier, long val, reduce_t op, int th_index) {
  int my_phase,i;
  long temp;

  pthread_mutex_lock(&nbarrier->lock);
  my_phase = nbarrier->phase;
  nbarrier->result[th_index] = val; 

  nbarrier->n_waiting++;
  if (nbarrier->n_waiting == nbarrier->n_clients) {/*get the prefix*/
    switch (op) {
    case MIN : temp = nbarrier->result[0];
      for(i = 1; i < nbarrier->n_clients;i++) {
         temp  = min(nbarrier->result[i],temp);
         nbarrier->result[i] = temp;
      }  
      break;
    case MAX : temp = nbarrier->result[0];
      for(i = 1; i < nbarrier->n_clients;i++) {
         temp  = max(nbarrier->result[i],temp);
         nbarrier->result[i] = temp;
      }  
      break;
    case SUM : 
      for(i = 1; i < nbarrier->n_clients;i++) 
         nbarrier->result[i] += nbarrier->result[i-1];
      break;
    default  : perror("ERROR: _SWARM_MULTICORE_scan_i() Bad reduction operator");
    }
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 1 - my_phase;
    pthread_cond_broadcast(&nbarrier->wait_cv);
  }
  while (nbarrier->phase == my_phase) {
    pthread_cond_wait(&nbarrier->wait_cv, &nbarrier->lock);
  }
  pthread_mutex_unlock(&nbarrier->lock);
  return(nbarrier->result[th_index]);
}

_SWARM_MULTICORE_scan_d_t _SWARM_MULTICORE_scan_init_d(int n_clients) {
  _SWARM_MULTICORE_scan_d_t nbarrier = (_SWARM_MULTICORE_scan_d_t)
    malloc(sizeof(struct _SWARM_MULTICORE_scan_d_s));
  assert_malloc(nbarrier);

  if (nbarrier != NULL) {
    nbarrier->n_clients = n_clients;
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 0;
    nbarrier->result    = (double *)malloc(n_clients*sizeof(double));
    assert_malloc(nbarrier->result);
    pthread_mutex_init(&nbarrier->lock, NULL);
    pthread_cond_init(&nbarrier->wait_cv, NULL);
  }
  return(nbarrier);
}

void _SWARM_MULTICORE_scan_destroy_d(_SWARM_MULTICORE_scan_d_t nbarrier) {

  pthread_mutex_destroy(&nbarrier->lock);
  pthread_cond_destroy(&nbarrier->wait_cv);
  free(nbarrier->result);
  free(nbarrier);
}

double _SWARM_MULTICORE_scan_d(_SWARM_MULTICORE_scan_d_t nbarrier, double val, reduce_t op,int th_index) {
  int my_phase,i;
  double temp;

  pthread_mutex_lock(&nbarrier->lock);
  my_phase = nbarrier->phase;
  nbarrier->result[th_index] = val;
  nbarrier->n_waiting++;
  if (nbarrier->n_waiting == nbarrier->n_clients) {
    switch (op) {
    case MIN : temp = nbarrier->result[0];
      for(i = 1; i < nbarrier->n_clients;i++) {
         temp  = min(nbarrier->result[i],temp);
         nbarrier->result[i] = temp;
      }  
      break;
    case MAX : temp = nbarrier->result[0];
      for(i = 1; i < nbarrier->n_clients;i++) {
         temp  = max(nbarrier->result[i],temp);
         nbarrier->result[i] = temp;
      }  
      break;
    case SUM : 
      for(i = 1; i < nbarrier->n_clients;i++) 
         nbarrier->result[i] += nbarrier->result[i-1];
      break;
    default  : perror("ERROR: _SWARM_MULTICORE_scan_i() Bad reduction operator");
    }
    nbarrier->n_waiting = 0;
    nbarrier->phase     = 1 - my_phase;
    pthread_cond_broadcast(&nbarrier->wait_cv);
  }
  while (nbarrier->phase == my_phase) {
    pthread_cond_wait(&nbarrier->wait_cv, &nbarrier->lock);
  }
  pthread_mutex_unlock(&nbarrier->lock);
  return(nbarrier->result[th_index]);
}



_SWARM_MULTICORE_spin_barrier_t _SWARM_MULTICORE_spin_barrier_init(int n_clients) {
  _SWARM_MULTICORE_spin_barrier_t sbarrier = (_SWARM_MULTICORE_spin_barrier_t)
    malloc(sizeof(struct _SWARM_MULTICORE_spin_barrier));
  assert_malloc(sbarrier);

  if (sbarrier != NULL) {
    sbarrier->n_clients = n_clients;
    sbarrier->n_waiting = 0;
    sbarrier->phase     = 0;
    pthread_mutex_init(&sbarrier->lock, NULL);
  }
  return(sbarrier);
}

void _SWARM_MULTICORE_spin_barrier_destroy(_SWARM_MULTICORE_spin_barrier_t sbarrier) {
  pthread_mutex_destroy(&sbarrier->lock);
  free(sbarrier);
}

void _SWARM_MULTICORE_spin_barrier_wait(_SWARM_MULTICORE_spin_barrier_t sbarrier) {
  int my_phase;

  while (pthread_mutex_trylock(&sbarrier->lock) == EBUSY) ;
  my_phase = sbarrier->phase;
  sbarrier->n_waiting++;
  if (sbarrier->n_waiting == sbarrier->n_clients) {
    sbarrier->n_waiting = 0;
    sbarrier->phase     = 1 - my_phase;
  }
  pthread_mutex_unlock(&sbarrier->lock);

  while (sbarrier->phase == my_phase) ;
}


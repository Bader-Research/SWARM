#ifndef _LOCK_H_
#define _LOCK_H_

#include <pthread.h>
#include "swarm.h"

typedef pthread_mutex_t* LOCK_t;

LOCK_t lock_init(THREADED);
void lock_it(LOCK_t lock);
void unlock_it(LOCK_t lock);
void lock_destroy(LOCK_t lock,THREADED);

#endif



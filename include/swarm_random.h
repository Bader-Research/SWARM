#ifndef _SWARM_RANDOM_H
#define _SWARM_RANDOM_H

#include <swarm.h>

void SWARM_DLL SWARM_random_init(THREADED);
void SWARM_DLL SWARM_random_destroy(THREADED);
void SWARM_DLL SWARM_srandom(unsigned int, THREADED);
long SWARM_DLL SWARM_random(THREADED);

void SWARM_DLL SWARM_srandomBit(unsigned int, THREADED);
int  SWARM_DLL SWARM_randomBit(THREADED);

#endif


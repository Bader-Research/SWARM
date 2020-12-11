# SWARM: Software and Algorithms for Running on Multicore

SWARM has been introduced as an open source parallel programming
framework. It is a library of primitives that fully exploit the
multicore processors. SWARM is built on POSIX threads that allow the
user to use either the already developed primitives or direct thread
primitives. SWARM has constructs for parallelization, restricting
control of threads, allocation and de-allocation of shared memory, and
communication primitives for synchronization, replication and
broadcasting. The framework has been successfully used to implement
efficient parallel versions of primitive algorithms. Viz. List
ranking, Prefix sums, Symmetry breaking, etc. SWARM is available for
Linux, FreeBSD, Solaris, AIX, and Microsoft Windows.

References:

D.A. Bader, V.N. Kanade, and K. Madduri, "SWARM: A Parallel
Programming Framework for Multi-Core Processors," First Workshop on
Multithreaded Architectures and Applications (MTAAP), Long Beach, CA,
March 30, 2007.

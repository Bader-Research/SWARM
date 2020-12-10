## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include 

TESTS = swarm_test

check_PROGRAMS = swarm_test

swarm_test_SOURCES = swarm_test.c  

swarm_test_DYNAMIC_LIBS = swarm

swarm_test_LDADD = -L$(top_srcdir)/src -l$(swarm_test_DYNAMIC_LIBS) 

swarm_test_LDFLAGS = -rpath @libdir@


CLEANFILES = *~

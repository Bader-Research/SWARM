## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include 

pkglib_PROGRAMS = SpanningTree

SpanningTree_SOURCES = elim_2d.c init.c lock.c span.c span_CREW.c span_breadth.c span_breadth_B.c span_breadth_B1.c span_breadth_B2.c span_e2d.c stack.c

noinst_HEADERS = graph.h lock.h stack.h

SpanningTree_DYNAMIC_LIBS = swarm

SpanningTree_LDADD = -lm -L$(top_srcdir)/src -l$(SpanningTree_DYNAMIC_LIBS) 

SpanningTree_LDFLAGS = -rpath @libdir@


CLEANFILES = *~

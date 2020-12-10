## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include 

#TESTS = isort64

pkglib_PROGRAMS = isort64

isort64_SOURCES = create_input.c  main.c  nas_r.c radixsort.c

noinst_HEADERS = create_input.h  nas_r.h  radixsort.h

isort64_DYNAMIC_LIBS = swarm

isort64_LDADD = -lm -L$(top_srcdir)/src -l$(isort64_DYNAMIC_LIBS) 

isort64_LDFLAGS = -rpath @libdir@


CLEANFILES = *~

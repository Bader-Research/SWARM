## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include 

pkglib_PROGRAMS = Listrank

Listrank_SOURCES = listrank.c main.c

noinst_HEADERS = listrank.h

Listrank_DYNAMIC_LIBS = swarm

Listrank_LDADD = -lm -L$(top_srcdir)/src -l$(Listrank_DYNAMIC_LIBS) 

Listrank_LDFLAGS = -rpath @libdir@


CLEANFILES = *~

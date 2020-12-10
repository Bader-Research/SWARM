## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include 

pkglib_PROGRAMS = Rake

Rake_SOURCES = rake.c listrank.c

noinst_HEADERS = listrank.h

Rake_DYNAMIC_LIBS = swarm

Rake_LDADD = -lm -L$(top_srcdir)/src -l$(Rake_DYNAMIC_LIBS) 

Rake_LDFLAGS = -rpath @libdir@


CLEANFILES = *~

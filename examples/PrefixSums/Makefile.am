## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include 

pkglib_PROGRAMS = PrefixSums

PrefixSums_SOURCES = prefixSums.c

#noinst_HEADERS =

PrefixSums_DYNAMIC_LIBS = swarm

PrefixSums_LDADD = -lm -L$(top_srcdir)/src -l$(PrefixSums_DYNAMIC_LIBS) 

PrefixSums_LDFLAGS = -rpath @libdir@


CLEANFILES = *~

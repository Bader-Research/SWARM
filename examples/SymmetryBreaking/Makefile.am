## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include 

pkglib_PROGRAMS = SymmetryBreaking

SymmetryBreaking_SOURCES = symbreak.c

noinst_HEADERS = symbreak.h

SymmetryBreaking_DYNAMIC_LIBS = swarm

SymmetryBreaking_LDADD = -lm -L$(top_srcdir)/src -l$(SymmetryBreaking_DYNAMIC_LIBS) 

SymmetryBreaking_LDFLAGS = -rpath @libdir@


CLEANFILES = *~

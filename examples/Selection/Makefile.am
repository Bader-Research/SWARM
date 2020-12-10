## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include 

pkglib_PROGRAMS = Selection

Selection_SOURCES = inputs_select.c select_fr.c select_randsamp.c select_seq.c smp_sel.c sorting.c 

noinst_HEADERS = inputs_select.h select_randsamp.h select_seq.h smp_sel.h sorting.h

Selection_DYNAMIC_LIBS = swarm

Selection_LDADD = -lm -L$(top_srcdir)/src -l$(Selection_DYNAMIC_LIBS) 

Selection_LDFLAGS = -rpath @libdir@


CLEANFILES = *~

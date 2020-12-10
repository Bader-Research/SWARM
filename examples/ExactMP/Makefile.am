## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include 

pkglib_PROGRAMS = ExactMP

ExactMP_SOURCES = addOrder.c aggGreedy.c eckDay.c frameWork.c main.c nj.c oBnB.c oDispTree.c oTaskGen.c preProc.c randOpt.c tbr.c uBnB.c uCompCost.c uDataStr.c uPLB.c uTaskGen.c

noinst_HEADERS = header.h

ExactMP_DYNAMIC_LIBS = swarm

ExactMP_LDADD = -lm -L$(top_srcdir)/src -l$(ExactMP_DYNAMIC_LIBS) 

ExactMP_LDFLAGS = -rpath @libdir@


CLEANFILES = *~

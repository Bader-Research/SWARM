## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

## library to be built using libtool. The `pkglib_' prefix means that
## the library will be installed in ${prefix}/lib/${package}. Use the
## `lib_' prefix to install libraries in ${prefix}/lib.
#pkglib_LTLIBRARIES = libswarm.la
lib_LTLIBRARIES = libswarm.la

## list of files to compile in order to produce the library
libswarm_la_SOURCES = swarm_random.c swarm_multicore.c swarm.c

## we want to link with pthread, and also with sprng, if it was
## detected by configure
libswarm_la_LDFLAGS = -release $(VERSION) @MINGW_LDFLAGS@ @SPRNG_LDFLAGS@ @SPRNG_LIBS@

## add the sprng include directory to the CFLAGS
libswarm_la_CFLAGS = -DSWARM_BUILD @SPRNG_CFLAGS@ 

## Include files that are going to be a part of the distribution.  The
## `pkginclude_' prefix means that these files will be installed in
## ${prefix}/include/${package}. Use the `include_' prefix to install
## header files in ${prefix}/include.
#pkginclude_HEADERS = swarm.h
#include_HEADERS = $(wildcard $(top_srcdir)/include/*.h)
include_HEADERS = $(top_srcdir)/include/swarm.h $(top_srcdir)/include/swarm_random.h

## Include files that are NOT going to be a part of the distribution.
noinst_HEADERS = $(top_srcdir)/include/swarm_multicore.h

## To ensure that in Windows the SWARM DLL can be called by a
## Microsoft Visual C++ program, an import library (.lib) must be
## generated, using Microsoft's LIB tool:
%.lib: %.dll
	lib /def:$<-def /machine:i386 /out:$@



CLEANFILES = *~ 

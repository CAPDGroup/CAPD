# -*- Makefile -*-

AM_CPPFLAGS = -I${top_srcdir}/src -I${top_srcdir}/src/rts \
              -I${top_srcdir}/src/fi_lib -I${top_srcdir}/src/asm

pkglibexec_PROGRAMS = example lexample io inewton linewton rungekutta trace

LDADD = ../src/libcxsc.la

example_SOURCES = example.cpp
lexample_SOURCES = lexample.cpp
io_SOURCES = io.cpp
inewton_SOURCES = inewton.cpp
linewton_SOURCES = linewton.cpp
rungekutta_SOURCES = rungekutta.cpp
trace_SOURCES = trace.cpp

# -*- Makefile -*-

AM_CPPFLAGS = -I${top_srcdir}/src -I${top_srcdir}/src/rts \
              -I${top_srcdir}/src/fi_lib -I${top_srcdir}/src/asm

noinst_PROGRAMS = benchmark timing

LDADD = ../src/libcxsc.la

benchmark_SOURCES = benchmark.cpp
timing_SOURCES = timing.cpp

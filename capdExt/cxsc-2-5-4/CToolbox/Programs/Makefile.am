# -*- Makefile -*-

AM_CPPFLAGS = -DNO_CXSC_STATUS -I${top_srcdir}/CToolbox/Modules \
              -I${top_srcdir}/src -I${top_srcdir}/src/rts \
              -I${top_srcdir}/src/fi_lib -I${top_srcdir}/src/asm

pkglibexec_PROGRAMS = \
	cpz_ex ddf_ex fastlss_ex gop_ex gop1_ex hess_ex jac_ex lop_ex lss_ex \
	nlfz_ex nlss_ex rpe_ex xev_ex1 xev_ex2 xev_ex3

LDADD = ../../src/libcxsc.la

cpz_ex_SOURCES = cpz_ex.cpp
ddf_ex_SOURCES = ddf_ex.cpp
gop_ex_SOURCES = gop_ex.cpp
gop1_ex_SOURCES = gop1_ex.cpp
hess_ex_SOURCES = hess_ex.cpp
jac_ex_SOURCES = jac_ex.cpp
lop_ex_SOURCES = lop_ex.cpp
lss_ex_SOURCES = lss_ex.cpp
nlfz_ex_SOURCES = nlfz_ex.cpp
nlss_ex_SOURCES = nlss_ex.cpp
rpe_ex_SOURCES = rpe_ex.cpp
xev_ex1_SOURCES = xev_ex1.cpp
xev_ex2_SOURCES = xev_ex2.cpp
xev_ex3_SOURCES = xev_ex3.cpp
fastlss_ex_SOURCES = fastlss_ex.cpp

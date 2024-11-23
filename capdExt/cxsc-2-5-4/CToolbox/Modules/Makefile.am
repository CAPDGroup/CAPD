# -*- Makefile -*-

AM_CPPFLAGS = -DNO_CXSC_STATUS -I${top_srcdir}/src -I${top_srcdir}/src/rts \
              -I${top_srcdir}/src/fi_lib -I${top_srcdir}/src/asm

noinst_LTLIBRARIES = libCToolbox.la

libCToolbox_la_SOURCES = \
	ci_util.cpp cipoly.cpp cpoly.cpp cpzero.cpp ddf_ari.cpp expreval.cpp \
	fastlss.cpp gop.cpp gop1.cpp grad_ari.cpp hess_ari.cpp linsys.cpp \
	lop.cpp lop_ari.cpp lst1_ari.cpp lst_ari.cpp matinv.cpp nlfzero.cpp \
	nlinsys.cpp rev_simp.cpp rpeval.cpp rpoly.cpp set_ari.cpp xi_ari.cpp

install-data-local:
	mkdir -p "${DESTDIR}${pkgincludedir}"
	${INSTALL} -pm0644 *.hpp *.inl "${DESTDIR}${pkgincludedir}/"

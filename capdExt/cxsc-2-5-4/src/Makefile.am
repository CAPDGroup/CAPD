# -*- Makefile -*-

SUBDIRS = rts fi_lib ../CToolbox/Modules .

AM_CPPFLAGS = -DNO_CXSC_STATUS -I${srcdir}/rts -I${srcdir}/fi_lib \
              -I${srcdir}/asm

lib_LTLIBRARIES = libcxsc.la
noinst_LTLIBRARIES = libcxsc_inl.la

libcxsc_inl_la_SOURCES = \
	cdot.cpp cidot.cpp cimath.cpp cinterval.cpp complex.cpp dot.cpp \
	dotio.cpp except.cpp idot.cpp imath.cpp interval.cpp ioflags.cpp \
	l_cimath.cpp l_cinterval.cpp l_cmath.cpp l_complex.cpp l_imath.cpp \
	l_interval.cpp l_real.cpp l_rmath.cpp lx_cinterval.cpp \
	lx_civector.cpp lx_complex.cpp lx_interval.cpp lx_ivector.cpp \
	lx_real.cpp real.cpp realio.cpp rmath.cpp

libcxsc_la_SOURCES = \
	cimatrix.cpp civector.cpp cmatrix.cpp cvector.cpp imatrix.cpp \
	intmatrix.cpp intvector.cpp ivector.cpp l_imatrix.cpp l_ivector.cpp \
	l_rmatrix.cpp l_rvector.cpp rmatrix.cpp rvector.cpp scimatrix.cpp \
	scivector.cpp scmatrix.cpp scvector.cpp simatrix.cpp sivector.cpp \
	srmatrix.cpp srvector.cpp

if HOST_DARWIN
if USE_ASM_X86_64
libcxsc_la_SOURCES += asm/r_ari_x86_64_mac.s
endif
if USE_ASM_X86
libcxsc_la_SOURCES += asm/r_ari_x86_mac.s
endif
else !HOST_DARWIN
if USE_ASM_X86_64
libcxsc_la_SOURCES += asm/r_ari_x86_64.s
endif
if USE_ASM_X86
libcxsc_la_SOURCES += asm/r_ari_x86.s
endif
endif !HOST_DARWIN

libcxsc_la_LIBADD  = rts/librts.la fi_lib/libfi_lib.la \
                     ../CToolbox/Modules/libCToolbox.la libcxsc_inl.la
libcxsc_la_LDFLAGS = -version-number 2:5:3

install-data-local:
	mkdir -p "${DESTDIR}${pkgincludedir}"
	install -pm0644 *.hpp *.inl *.h rts/*.h fi_lib/*.hpp asm/r_ari.h \
		../cxscconf.h "${DESTDIR}${pkgincludedir}/"

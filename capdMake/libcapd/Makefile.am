include ${capdMake}/make/common_makefile.mkf
ACLOCAL_AMFLAGS= -I ../ac_aux/m4

COMMON_LDFLAGS=-version-info ${CAPD_VERSION_INFO}
COMMON_LIBS=

if CAPD_NO_UNDEFINED
   COMMON_LDFLAGS+= -no-undefined
endif


pkgconfigdir = @PKG_CONFIG_DIR@

bin_SCRIPTS=
pkgconfig_DATA=
lib_LTLIBRARIES=


if HAVE_CAPD

noinst_LTLIBRARIES=libcapdweaksymbols.la
libcapdweaksymbols_la_SOURCES=main_entry.cpp
libcapdweaksymbols_la_CXXFLAGS= ${AM_CXXFLAGS} ${capd_CFLAGS} -I${top_srcdir}/../../capdAlg/include -I${top_srcdir}/../../capdAux/include
libcapdweaksymbols_la_LDFLAGS=

bin_SCRIPTS+= capd-config
pkgconfig_DATA+= capd.pc
lib_LTLIBRARIES+=libcapd.la

libcapd_la_SOURCES=
libcapd_la_LIBADD=${capd_LIBS} ${COMMON_LIBS}
libcapd_la_LDFLAGS=${COMMON_LDFLAGS} -weak_library libcapdweaksymbols.la
libcapd_la_DEPENDENCIES=libcapdweaksymbols.la ${filter %.a %.la,${libcapd_la_LIBADD}}

endif

if HAVE_CAPD_GUI

bin_SCRIPTS+= capd-gui-config
pkgconfig_DATA+= capd-gui.pc
lib_LTLIBRARIES+=libcapd-gui.la

libcapd_gui_la_SOURCES=
libcapd_gui_la_LIBADD=${capd_gui_LIBS} ${COMMON_LIBS}
libcapd_gui_la_LDFLAGS= ${COMMON_LDFLAGS} -weak_library libcapdweaksymbols.la
libcapd_gui_la_DEPENDENCIES=libcapdweaksymbols.la ${filter %.a %.la,${libcapd_gui_la_LIBADD}}

endif

if HAVE_MPCAPD

bin_SCRIPTS+= mpcapd-config
pkgconfig_DATA+= mpcapd.pc
lib_LTLIBRARIES+=libmpcapd.la

libmpcapd_la_SOURCES=
libmpcapd_la_LIBADD=${mpcapd_LIBS} ${COMMON_LIBS}
libmpcapd_la_LDFLAGS=${COMMON_LDFLAGS} -weak_library libcapdweaksymbols.la
libmpcapd_la_DEPENDENCIES=libcapdweaksymbols.la ${filter %.a %.la,${libmpcapd_la_LIBADD}}

endif

if HAVE_MPCAPD_GUI

bin_SCRIPTS+= mpcapd-gui-config
pkgconfig_DATA+= mpcapd-gui.pc
lib_LTLIBRARIES+=libmpcapd-gui.la

libmpcapd_gui_la_SOURCES=
libmpcapd_gui_la_LIBADD=${mpcapd_gui_LIBS} ${COMMON_LIBS}
libmpcapd_gui_la_LDFLAGS=${COMMON_LDFLAGS}  -weak_library libcapdweaksymbols.la
libmpcapd_gui_la_DEPENDENCIES=libcapdweaksymbols.la ${filter %.a %.la,${libmpcapd_gui_la_LIBADD}}

endif

all: chmodScripts

chmodScripts:
	chmod +x ${top_capd_build}/bin/*

.PHONY: chmodScripts

.PHONY: lib

lib: chmodScripts

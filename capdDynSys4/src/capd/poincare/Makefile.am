include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapdpoincare.la
libcapdpoincare_la_SOURCES= AbstractSection.cpp TimeMap.cpp PoincareMap.cpp
libcapdpoincare_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapdpoincare_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capdpoincare_CFLAGS} ${capdAlg_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

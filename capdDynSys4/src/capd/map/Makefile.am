include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapdmap.la
libcapdmap_la_SOURCES= Function.cpp  Map.cpp Parser.cpp
libcapdmap_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapdmap_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capdmap_CFLAGS} ${capdAlg_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

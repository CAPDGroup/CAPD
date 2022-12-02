include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapdautodiff.la
libcapdautodiff_la_SOURCES= DagIndexer.cpp NodeType.cpp
libcapdautodiff_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapdautodiff_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capdautodiff_CFLAGS} ${capdAlg_CFLAGS} ${capdAux_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

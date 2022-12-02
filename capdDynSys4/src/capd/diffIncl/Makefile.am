include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapddiffIncl.la
libcapddiffIncl_la_SOURCES=DiffInclusion.cpp  DiffInclusionCW.cpp  DiffInclusionLN.cpp  InclRect2Set.cpp   MultiMap.cpp
libcapddiffIncl_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapddiffIncl_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capddiffIncl_CFLAGS} ${capdAlg_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

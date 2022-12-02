include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapdcovrel.la
libcapdcovrel_la_SOURCES=GraphicsHSet3D.cpp  GraphicsTripleSet.cpp  HSet2D.cpp  HSet3D.cpp  HSetMD.cpp  HSetND.cpp  TripleSet.cpp
libcapdcovrel_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapdcovrel_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capdcovrel_CFLAGS} ${capdAlg_CFLAGS} ${capdAux_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

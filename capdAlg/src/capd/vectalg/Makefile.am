include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapdvectalg.la
libcapdvectalg_la_SOURCES=Container.cpp  DVector.cpp  IVector.cpp LDVector.cpp   Makefile.am Norm.cpp DMatrix.cpp    IMatrix.cpp  LDMatrix.cpp Multiindex.cpp  ZVector.cpp
libcapdvectalg_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapdvectalg_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capdvectalg_CFLAGS} ${capdAux_CFLAGS}  ${capdExt_CFLAGS}




pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

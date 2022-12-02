include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapdgeomset.la
libcapdgeomset_la_SOURCES=AffineSet.cpp  DoubletonSet.cpp  MatrixAffineSet.cpp  MatrixDoubletonSet.cpp TripletonSet.cpp
libcapdgeomset_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapdgeomset_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capdgeomset_CFLAGS} ${capdAlg_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

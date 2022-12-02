include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapddynset.la
libcapddynset_la_SOURCES = C0BallSet.cpp C0AffineSet.cpp C0DoubletonSet.cpp C0HOSet.cpp C0TripletonSet.cpp CoordReorganize.cpp C0FlowballSet.cpp FactorReorganization.cpp C1AffineSet.cpp C1DoubletonSet.cpp C1HOSet.cpp C11Rect2Set.cpp C2DoubletonSet.cpp CnRect2Set.cpp CnDoubletonSet.cpp
libcapddynset_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapddynset_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capddynset_CFLAGS} ${capdAlg_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

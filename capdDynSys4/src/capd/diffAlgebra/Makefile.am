include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapddiffAlgebra.la
libcapddiffAlgebra_la_SOURCES= C1TimeJet.cpp C2TimeJet.cpp CnTimeJet.cpp  Curve.cpp Hessian.cpp SolutionCurve.cpp
libcapddiffAlgebra_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapddiffAlgebra_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capddiffAlgebra_CFLAGS} ${capdAlg_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

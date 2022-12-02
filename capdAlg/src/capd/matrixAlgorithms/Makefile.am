include ${capdMake}/make/common_makefile.mkf


noinst_LTLIBRARIES=libcapdmatrixAlgorithms.la
libcapdmatrixAlgorithms_la_SOURCES= floatMatrixAlgorithms.cpp  intMatrixAlgorithms.cpp PARIConvert.h PARIInterface.cpp PARIInterface.hpp SmithFormFactory.cpp SmithFormFactory.hpp PARISmithForm.cpp PARISmithForm.hpp
libcapdmatrixAlgorithms_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapdmatrixAlgorithms_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capdmatrixAlgorithms_CFLAGS} ${capdAux_CFLAGS} ${capdExt_CFLAGS} ${mpcapdAlg_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =

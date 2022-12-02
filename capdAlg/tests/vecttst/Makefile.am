include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS=matrixMDTest vectMDTest vecttst
TESTS=matrixMDTest vectMDTest

if HAVE_BOOST
   check_PROGRAMS+=vectalgTest
   TESTS+=vectalgTest
   check_PROGRAMS+=C11vectalgTest
   TESTS+=C11vectalgTest
endif

matrixMDTest_SOURCES=matrixMDTest.cpp
matrixMDTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdAlg_CFLAGS}
matrixMDTest_libs=${capdAlg_LIBS}
matrixMDTest_LDADD=${matrixMDTest_libs}
matrixMDTest_DEPENDENCIES=${filter %.a %.la,${matrixMDTest_libs}}

vectMDTest_SOURCES=vectMDTest.cpp
vectMDTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdAlg_CFLAGS}
vectMDTest_libs=${capdAlg_LIBS}
vectMDTest_LDADD=${vectMDTest_libs}
vectMDTest_DEPENDENCIES=${filter %.a %.la,${vectMDTest_libs}}

vecttst_SOURCES=vecttst.cpp
vecttst_CXXFLAGS= ${AM_CXXFLAGS} ${capdAlg_CFLAGS}
vecttst_libs=${capdAlg_LIBS}
vecttst_LDADD=${vecttst_libs}
vecttst_DEPENDENCIES=${filter %.a %.la,${vecttst_libs}}

vectalgTest_SOURCES=vectalgTestMain.cpp StaticArrayInitTest.cpp
vectalgTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdAlg_CFLAGS} ${BOOST_CPPFLAGS}
vectalgTest_libs=${capdAlg_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
vectalgTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
vectalgTest_LDADD=${vectalgTest_libs}
vectalgTest_DEPENDENCIES=${filter %.a %.la,${vectalgTest_libs}}

C11vectalgTest_SOURCES=C11vectalgTestMain.cpp C11vectalgInitTest.cpp 
C11vectalgTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdAlg_CFLAGS} ${BOOST_CPPFLAGS} ${CAPD_CXXFLAGS_CXX11} 
C11vectalgTest_libs=${capdAlg_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
C11vectalgTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
C11vectalgTest_LDADD=${C11vectalgTest_libs}
C11vectalgTest_DEPENDENCIES=${filter %.a %.la,${C11vectalgTest_libs}}

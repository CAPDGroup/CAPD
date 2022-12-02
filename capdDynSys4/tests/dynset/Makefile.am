include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS=AffineSet C0Pped2SetTest Rect2SetTest \
               CoordReorganization    FactorReorganization\
            CenteredAffineSetTest DoubletonSetTest C0RectSetTest C0PpedSetTest

TESTS=${check_PROGRAMS}

CAPD_AM_CXXFLAGS+=${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
LDADD=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
LDFLAGS+=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}

#-lboost_unit_test_framework

AffineSet_SOURCES= AffineSetTest.cpp AffineSetCommonTest.hpp AffineSetSpecific.hpp
AffineSet_DEPENDENCIES=${filter %.a %.la,${LDADD}}

C0Pped2SetTest_SOURCES= C0Pped2SetTest.cpp
C0Pped2SetTest_DEPENDENCIES=${filter %.a %.la,${LDADD}}

Rect2SetTest_SOURCES= Rect2SetTest.cpp
Rect2SetTest_DEPENDENCIES=${filter %.a %.la,${LDADD}}

CoordReorganization_SOURCES= CoordReorganization.cpp
CoordReorganization_DEPENDENCIES=${filter %.a %.la,${LDADD}}

FactorReorganization_SOURCES= FactorReorganization.cpp FactorReorganizationSpecific.hpp
FactorReorganization_DEPENDENCIES=${filter %.a %.la,${LDADD}}

CenteredAffineSetTest_SOURCES= CenteredAffineSetTest.cpp CenteredAffineSetSpecific.hpp
CenteredAffineSetTest_DEPENDENCIES=${filter %.a %.la,${LDADD}}

DoubletonSetTest_SOURCES= DoubletonSetTest.cpp CenteredDoubletonSpecific.hpp DoubletonSetSpecific.hpp
DoubletonSetTest_DEPENDENCIES=${filter %.a %.la,${LDADD}}

C0RectSetTest_SOURCES=C0RectSetTest.cpp
C0RectSetTest_DEPENDENCIES=${filter %.a %.la,${LDADD}}

C0PpedSetTest_SOURCES= C0PpedSetTest.cpp
C0PpedSetTest_DEPENDENCIES=${filter %.a %.la,${LDADD}}

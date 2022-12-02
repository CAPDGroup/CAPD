include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS=
TESTS=
if HAVE_BOOST
   check_PROGRAMS+=TexWriterTest
   TESTS+=TexWriterTest
endif

TexWriterTest_SOURCES=TexWriterTest.cpp 
TexWriterTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdAlg_CFLAGS} ${BOOST_CPPFLAGS}
TexWriterTest_libs=${capdAlg_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
TexWriterTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
TexWriterTest_LDADD=${TexWriterTest_libs}
TexWriterTest_DEPENDENCIES=${filter %.a %.la,${TexWriterTest_libs}}

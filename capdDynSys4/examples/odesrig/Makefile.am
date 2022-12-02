include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=ITimeMapExample ITimeMapEncloseTrajectoryExample ITimeMapSolutionCurveExample ITimeMapVariationalEquationsExample ITimeMapMonodromyMatrixCurveExample

ITimeMapExample_SOURCES=ITimeMapExample.cpp
ITimeMapExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
ITimeMapExample_libs=${capdDynSys_LIBS}
ITimeMapExample_LDADD=${ITimeMapExample_libs}
ITimeMapExample_DEPENDENCIES=${filter %.a %.la,${ITimeMapExample_libs}}

ITimeMapEncloseTrajectoryExample_SOURCES=ITimeMapEncloseTrajectoryExample.cpp
ITimeMapEncloseTrajectoryExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
ITimeMapEncloseTrajectoryExample_libs=${capdDynSys_LIBS}
ITimeMapEncloseTrajectoryExample_LDADD=${ITimeMapEncloseTrajectoryExample_libs}
ITimeMapEncloseTrajectoryExample_DEPENDENCIES=${filter %.a %.la,${ITimeMapEncloseTrajectoryExample_libs}}

ITimeMapSolutionCurveExample_SOURCES=ITimeMapSolutionCurveExample.cpp
ITimeMapSolutionCurveExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
ITimeMapSolutionCurveExample_libs=${capdDynSys_LIBS}
ITimeMapSolutionCurveExample_LDADD=${ITimeMapSolutionCurveExample_libs}
ITimeMapSolutionCurveExample_DEPENDENCIES=${filter %.a %.la,${ITimeMapSolutionCurveExample_libs}}

ITimeMapVariationalEquationsExample_SOURCES=ITimeMapVariationalEquationsExample.cpp
ITimeMapVariationalEquationsExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
ITimeMapVariationalEquationsExample_libs=${capdDynSys_LIBS}
ITimeMapVariationalEquationsExample_LDADD=${ITimeMapVariationalEquationsExample_libs}
ITimeMapVariationalEquationsExample_DEPENDENCIES=${filter %.a %.la,${ITimeMapVariationalEquationsExample_libs}}

ITimeMapMonodromyMatrixCurveExample_SOURCES=ITimeMapMonodromyMatrixCurveExample.cpp
ITimeMapMonodromyMatrixCurveExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
ITimeMapMonodromyMatrixCurveExample_libs=${capdDynSys_LIBS}
ITimeMapMonodromyMatrixCurveExample_LDADD=${ITimeMapMonodromyMatrixCurveExample_libs}
ITimeMapMonodromyMatrixCurveExample_DEPENDENCIES=${filter %.a %.la,${ITimeMapMonodromyMatrixCurveExample_libs}}

EXTRA_DIST=odes_rigorous.dox odes_rigorous_timemap.dox odesvar_rigorous.dox 

#CMakeLists.txt Configuration File, 2012-12-11
cmake_minimum_required(VERSION 2.4)

project(cxsc C CXX)

#set(CMAKE_BUILD_TYPE release)

# the current C-XSC version
set(VER_MAJOR "2")
set(VER_MINOR "5")
set(VER_PATCH "4")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "C-XSC and the C++ Toolbox for Verified Computation")
set(CPACK_PACKAGE_VENDOR "Universitaet Wuppertal, Germany")
set(CPACK_PACKAGE_VERSION_MAJOR "${VER_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${VER_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${VER_PATCH}")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/docu/COPYING")
include(CPack)


# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


include_directories(
  # where cxscconf.h is
  ${CMAKE_CURRENT_BINARY_DIR}

  # the C-XSC sources
  ${CMAKE_SOURCE_DIR}/CToolbox/Modules
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/fi_lib
  ${CMAKE_SOURCE_DIR}/src/rts
  ${CMAKE_SOURCE_DIR}/src/asm
)

# Check for intel compiler
if (${CMAKE_C_COMPILER} MATCHES "icc") 
  set(USING_INTEL_COMPILER TRUE)
endif ()

# Rounding mode options
if(NOT WIN32)
  if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64" OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86" OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64")
    option(ROUNDING_ASM "Use assembler routines for rounding" ON)
    if(NOT ROUNDING_ASM)
      option(ROUNDING_HARD "Use hardware support for IEEE 759 arithmetic" ON)
      if(ROUNDING_HARD)
        option(ROUNDING_QUICK "Use the fast but critical hardware IEEE arithmetic" OFF)
      endif(ROUNDING_HARD)
    endif(NOT ROUNDING_ASM)
  endif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64" OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86" OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64")
else(NOT WIN32)
  option(ROUNDING_HARD "Use hardware support for IEEE 759 arithmetic" ON)
  if(ROUNDING_HARD)
    option(ROUNDING_QUICK "Use the fast but critical hardware IEEE arithmetic" OFF)
  endif(ROUNDING_HARD)
endif(NOT WIN32)

# Option for shared or static library (either/or)
if(NOT WIN32)
  if(APPLE)
    option(BUILD_SHARED "Build the shared version of the library" OFF)
  else(APPLE)
    option(BUILD_SHARED "Build the shared version of the library" ON)
  endif(APPLE)
endif(NOT WIN32)

# option(BUILD_ALL_FILES "Build all .c and .cpp files instead of only those specified by the CMakeLists.txt in the subdirectories" OFF)


# set the rounding options in cxscconf.h. The C96 options are only used
# in the combination "SunOS" with a GNU compiler (as it is done in the
# install_cxsc file). Since the #if preprocessor directive treats undefined
# macros as zero, we can omit setting the others to "0".
if(NOT WIN32)
  if(ROUNDING_ASM)
    set(ROUND_ASM 1)
  else(ROUNDING_ASM)
    if(ROUNDING_HARD)
      if(ROUNDING_QUICK)
        set(ROUND_C99_QUICK 1)
      else(ROUNDING_QUICK)
        set(ROUND_C99_SAVE 1)
      endif(ROUNDING_QUICK)
    endif(ROUNDING_HARD)
  endif(ROUNDING_ASM)
endif(NOT WIN32)

# System settings for the RTS (big/little endian, 64 Bit)
if(APPLE)

  # powerpc
  if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "powerpc")
    set(SUN4_GNU_C 1)
  # Default (e.g. i386)
  else(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "powerpc")
    set(IBM_LINUX_C 1)
    set(GNU_X86_64 1)   
    set(CXSC_x86_64 1)
  endif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "powerpc")

elseif(UNIX)

  # Linux systems
  if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    # alpha
    if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "alpha")
      set(DEC_ALPHA_C 1)
    # x86_64
    elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
      set(IBM_LINUX_C 1)
      set(GNU_X86_64 1)
      set(CXSC_x86_64 1)
    elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64")
      set(CXSC_PPC64 1)
      set(GNU_X86_64 1)
    # Default (e.g. i386, ia64)
    else(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "alpha")
      set(IBM_LINUX_C 1)
    endif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "alpha")

  # AIX
  elseif(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    set(IBM_RS_6000_C 1)
    set(GNU_X86_64 1)

  # SunOS
  elseif(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    set(SUN4_GNU_C 1)

  # Default if nothing else matches. Used for example in cygwin
  else(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(IBM_LINUX_C 1)
  endif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

elseif(MSVC)

  set(IBM_LINUX_C 1)
  set(WINDOWS_X86_32 1)

elseif(MINGW)

  set(IBM_LINUX_C 1)
  set(WINDOWS_MINGW 1)

else(APPLE)

  set(IBM_LINUX_C 1)

endif(APPLE)

# Create cxscconf.h file (will be created in the build directory)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cxscconf.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/cxscconf.h
)

# add the list of files to compile (excluding example and test programs)
include(${CMAKE_SOURCE_DIR}/src/CMakeLists.txt)
include(${CMAKE_SOURCE_DIR}/src/fi_lib/CMakeLists.txt)
include(${CMAKE_SOURCE_DIR}/src/rts/CMakeLists.txt)
include(${CMAKE_SOURCE_DIR}/CToolbox/Modules/CMakeLists.txt)


# set the variables for all header files
file(GLOB cxsc_headers
  "${CMAKE_SOURCE_DIR}/src/*.h"
  "${CMAKE_SOURCE_DIR}/src/*.hpp"
  "${CMAKE_SOURCE_DIR}/src/*.inl"
)
file(GLOB rts_headers
  "${CMAKE_SOURCE_DIR}/src/rts/*.h"
)
file(GLOB fi_lib_headers
  "${CMAKE_SOURCE_DIR}/src/fi_lib/*.hpp"
)
file(GLOB asm_headers
   "${CMAKE_SOURCE_DIR}/src/asm/*.h"
)
file(GLOB ctoolbox_headers
  "${CMAKE_SOURCE_DIR}/CToolbox/Modules/*.hpp"
  "${CMAKE_SOURCE_DIR}/CToolbox/Modules/*.inl"
)

#Set special options depending on compiler (optimization, error free transformations)
if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86" AND CMAKE_COMPILER_IS_GNUCXX)
  #Compiler is GCC, System is 32 Bit x86
  if(ROUNDING_ASM)
    #Full optimizations, FP computations using SSE
    set(CMAKE_CXX_FLAGS "-O3 -mfpmath=sse -msse2")
    set(CMAKE_CC_FLAGS "-O3 -mfpmath=sse -msse2")
  else(ROUNDING_ASM)
    #Without assembler: Reduced optimizations, FP computations using SSE
    set(CMAKE_CXX_FLAGS "-O1 -finline-functions -mfpmath=sse -msse2")
    set(CMAKE_CC_FLAGS "-O1 -finline-functions -mfpmath=sse -msse2")
  endif(ROUNDING_ASM)
elseif(CMAKE_COMPILER_IS_GNUCXX)
  #Compiler is GCC, 64 Bit x86 System
  if(ROUNDING_ASM)
    #Full optimizations, SSE registers used by default
    set(CMAKE_CXX_FLAGS "-O3")
    set(CMAKE_CC_FLAGS "-O3")
  else(ROUNDING_ASM)
    #Reduced optimizations, SSE registers used by default
    set(CMAKE_CXX_FLAGS "-O1 -finline-functions")
    set(CMAKE_CC_FLAGS "-O1 -finline-functions")
  endif(ROUNDING_ASM)
elseif(USING_INTEL_COMPILER)
  #Intel Compiler: Full optimizations, reduced FP optimzations for correct rounding and error free transformations
  set(CMAKE_CXX_FLAGS "-O3 -fp-model source")
  set(CMAKE_CC_FLAGS "-O3 -fp-model source")
else(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86" AND CMAKE_COMPILER_IS_GNUCXX)
  #Default compiler options
  if(ROUNDING_ASM)
    #With assembler rounding: Full optimizations
    set(CMAKE_CXX_FLAGS "-O1")
    set(CMAKE_CC_FLAGS "-O1")
  elseif(ROUNDING_ASM)
    #Without assembler rounding: Reduced optimizations
    set(CMAKE_CXX_FLAGS "-O1")
    set(CMAKE_CC_FLAGS "-O1")
  endif(ROUNDING_ASM)
endif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86" AND CMAKE_COMPILER_IS_GNUCXX)


#Add static or shared library to build
if(NOT WIN32)
  if(BUILD_SHARED)
    add_library(cxsc SHARED
      ${cxsc_src}
      ${rts_src}
      ${fi_lib_src}
      ${ctoolbox_src}
    )
  else(BUILD_SHARED)
    add_library(cxsc STATIC
      ${cxsc_src}
      ${rts_src}
      ${fi_lib_src}
      ${ctoolbox_src}
    )
  endif(BUILD_SHARED)
else(NOT WIN32)
  # on Windows only a static version of the library can be built.
  add_library(cxsc STATIC
    ${cxsc_src}
    ${rts_src}
    ${fi_lib_src}
    ${ctoolbox_src}
  )
endif(NOT WIN32)

# set the SO version of the library file
if(BUILD_SHARED)
  set_target_properties(cxsc PROPERTIES
    VERSION ${VER_MAJOR}.${VER_MINOR}.${VER_PATCH}
    SOVERSION ${VER_MAJOR}
  )
endif(BUILD_SHARED)

# build the examples
add_subdirectory(examples)
add_subdirectory(CToolbox/Programs)
add_subdirectory(CToolbox/test)
add_subdirectory(CToolbox/test/test_programms)

# installation of the header files
install(FILES
  ${cxsc_headers}
  ${rts_headers}
  ${asm_headers}
  ${fi_lib_headers}
  ${ctoolbox_headers}
  ${CMAKE_CURRENT_BINARY_DIR}/cxscconf.h
DESTINATION include)

# installation of the documentation
install(FILES changelog README DESTINATION doc)

# install the library
install(TARGETS cxsc LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

# the source_group command groups the files in IDEs, e.g. in
# Microsoft Visual C++.
source_group("Source Files" FILES ${cxsc_src})
source_group("Source Files\\rts" FILES ${rts_src})
source_group("Source Files\\fi_lib" FILES ${fi_lib_src})
source_group("Source Files\\CToolbox" FILES ${ctoolbox_src})
source_group("Header Files" FILES ${cxsc_headers})
source_group("Header Files\\rts" FILES ${rts_headers})
source_group("Header Files\\asm" FILES ${asm_headers})
source_group("Header Files\\fi_lib" FILES ${fi_lib_headers})
source_group("Header Files\\CToolbox" FILES ${ctoolbox_headers})

#CTests for Testing CXSC Library of correctness
SET (CTOOLBOX_PROGRAMS_DIR ${CMAKE_CURRENT_BINARY_DIR}/CToolbox/test/test_programms) 

SET (CTOOLBOX_INPUT_DIR ${CMAKE_SOURCE_DIR}/CToolbox/test/test_input) 
SET (CTOOLBOX_VLG_EXP_DIR ${CMAKE_SOURCE_DIR}/CToolbox/test/test_vgl_exp) 

#make test, runs some test programms
ENABLE_TESTING()

#cpz_test_1
ADD_TEST(cpz_test_1 ${CTOOLBOX_PROGRAMS_DIR}/cpz_test 
  ${CTOOLBOX_INPUT_DIR}/cpz_ex.in1
  ${CTOOLBOX_VLG_EXP_DIR}/cpz_ex.vgl1
  ${CTOOLBOX_VLG_EXP_DIR}/cpz_ex.exp1)

#cpz_test_2
ADD_TEST(cpz_test_2 ${CTOOLBOX_PROGRAMS_DIR}/cpz_test 
  ${CTOOLBOX_INPUT_DIR}/cpz_ex.in2
  ${CTOOLBOX_VLG_EXP_DIR}/cpz_ex.vgl2 
  ${CTOOLBOX_VLG_EXP_DIR}/cpz_ex.exp2)

#gop1_test
ADD_TEST(gop1_test ${CTOOLBOX_PROGRAMS_DIR}/gop1_test 
  ${CTOOLBOX_INPUT_DIR}/gop1_ex.in
  ${CTOOLBOX_VLG_EXP_DIR}/gop1_ex.vgl
  ${CTOOLBOX_VLG_EXP_DIR}/gop1_ex.exp)

#gop_test
ADD_TEST(gop_test ${CTOOLBOX_PROGRAMS_DIR}/gop_test 
  ${CTOOLBOX_INPUT_DIR}/gop_ex.in
  ${CTOOLBOX_VLG_EXP_DIR}/gop_ex.vgl
  ${CTOOLBOX_VLG_EXP_DIR}/gop_ex.exp)
  
#hess_test
ADD_TEST(hess_test ${CTOOLBOX_PROGRAMS_DIR}/hess_test 
  ${CTOOLBOX_INPUT_DIR}/hess_ex.in
  ${CTOOLBOX_VLG_EXP_DIR}/hess_ex.vgl
  ${CTOOLBOX_VLG_EXP_DIR}/hess_ex.exp)
  
#jac_test
ADD_TEST(jac_test ${CTOOLBOX_PROGRAMS_DIR}/jac_test 
  ${CTOOLBOX_INPUT_DIR}/jac_ex.in
  ${CTOOLBOX_VLG_EXP_DIR}/jac_ex.vgl
  ${CTOOLBOX_VLG_EXP_DIR}/jac_ex.exp)
  
#lop_test_1
ADD_TEST(lop_test_1 ${CTOOLBOX_PROGRAMS_DIR}/lop_test
  ${CTOOLBOX_INPUT_DIR}/lop_ex.in1
  ${CTOOLBOX_VLG_EXP_DIR}/lop_ex.vgl1
  ${CTOOLBOX_VLG_EXP_DIR}/lop_ex.exp1)

#lop_test_2
ADD_TEST(lop_test_2 ${CTOOLBOX_PROGRAMS_DIR}/lop_test
  ${CTOOLBOX_INPUT_DIR}/lop_ex.in2
  ${CTOOLBOX_VLG_EXP_DIR}/lop_ex.vgl2
  ${CTOOLBOX_VLG_EXP_DIR}/lop_ex.exp2)
  
#lss_test_1
ADD_TEST(lss_test_1 ${CTOOLBOX_PROGRAMS_DIR}/lss_test
  ${CTOOLBOX_INPUT_DIR}/lss_ex.in1
  ${CTOOLBOX_VLG_EXP_DIR}/lss_ex.vgl1
  ${CTOOLBOX_VLG_EXP_DIR}/lss_ex.exp1)
  
#lss_test_2
ADD_TEST(lss_test_2 ${CTOOLBOX_PROGRAMS_DIR}/lss_test
  ${CTOOLBOX_INPUT_DIR}/lss_ex.in2
  ${CTOOLBOX_VLG_EXP_DIR}/lss_ex.vgl2
  ${CTOOLBOX_VLG_EXP_DIR}/lss_ex.exp2)
  
#nlfz_test
ADD_TEST(nlfz_test ${CTOOLBOX_PROGRAMS_DIR}/nlfz_test
  ${CTOOLBOX_INPUT_DIR}/nlfz_ex.in
  ${CTOOLBOX_VLG_EXP_DIR}/nlfz_ex.vgl
  ${CTOOLBOX_VLG_EXP_DIR}/nlfz_ex.exp)
  
#/nlss_test_1
ADD_TEST(nlss_test_1 ${CTOOLBOX_PROGRAMS_DIR}/nlss_test
  ${CTOOLBOX_INPUT_DIR}/nlss_ex.in1
  ${CTOOLBOX_VLG_EXP_DIR}/nlss_ex.vgl1
  ${CTOOLBOX_VLG_EXP_DIR}/nlss_ex.exp1)

#/nlss_test_2
ADD_TEST(nlss_test_2 ${CTOOLBOX_PROGRAMS_DIR}/nlss_test
  ${CTOOLBOX_INPUT_DIR}/nlss_ex.in2
  ${CTOOLBOX_VLG_EXP_DIR}/nlss_ex.vgl2
  ${CTOOLBOX_VLG_EXP_DIR}/nlss_ex.exp2)
  
#/rpe_test_1
ADD_TEST(rpe_test_1 ${CTOOLBOX_PROGRAMS_DIR}/rpe_test
  ${CTOOLBOX_INPUT_DIR}/rpe_ex.in1
  ${CTOOLBOX_VLG_EXP_DIR}/rpe_ex.vgl1
  ${CTOOLBOX_VLG_EXP_DIR}/rpe_ex.exp1)

#/rpe_test_2
ADD_TEST(rpe_test_2 ${CTOOLBOX_PROGRAMS_DIR}/rpe_test
  ${CTOOLBOX_INPUT_DIR}/rpe_ex.in2
  ${CTOOLBOX_VLG_EXP_DIR}/rpe_ex.vgl2
  ${CTOOLBOX_VLG_EXP_DIR}/rpe_ex.exp2)
  
#/xev_test_1
ADD_TEST(xev_test_1 ${CTOOLBOX_PROGRAMS_DIR}/xev_test1
  ${CTOOLBOX_INPUT_DIR}/xev_ex1.in
  ${CTOOLBOX_VLG_EXP_DIR}/xev_ex1.vgl
  ${CTOOLBOX_VLG_EXP_DIR}/xev_ex1.exp)

#/xev_test_2
ADD_TEST(xev_test_2 ${CTOOLBOX_PROGRAMS_DIR}/xev_test2
  ${CTOOLBOX_INPUT_DIR}/xev_ex2.in
  ${CTOOLBOX_VLG_EXP_DIR}/xev_ex2.vgl
  ${CTOOLBOX_VLG_EXP_DIR}/xev_ex2.exp)

#/xev_test_3
ADD_TEST(xev_test_3 ${CTOOLBOX_PROGRAMS_DIR}/xev_test3
  ${CTOOLBOX_INPUT_DIR}/xev_ex3.in
  ${CTOOLBOX_VLG_EXP_DIR}/xev_ex3.vgl
  ${CTOOLBOX_VLG_EXP_DIR}/xev_ex3.exp)

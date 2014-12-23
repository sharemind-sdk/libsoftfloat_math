#
# This file is a part of the Sharemind framework.
# Copyright (C) Cybernetica AS
#
# All rights are reserved. Reproduction in whole or part is prohibited
# without the written consent of the copyright owner. The usage of this
# code is subject to the appropriate license agreement.
#

CMAKE_MINIMUM_REQUIRED(VERSION "2.8.12")
PROJECT("SHAREMIND_LIBSOFTFLOAT_MATH" "C")

INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/config.local" OPTIONAL)
INCLUDE("${CMAKE_CURRENT_BINARY_DIR}/config.local" OPTIONAL)

ADD_COMPILE_OPTIONS(
  "-std=c99" "-O2" "-fstrict-aliasing" "-fwrapv" "-fPIC" "-Wextra" "-Wall"
  "-Wstrict-aliasing=1" "-Wabi" "-Wbad-function-cast" "-Wc++-compat" "-Wformat"
  "-Wswitch-default" "-Wunused" "-Wfloat-equal" "-Wshadow" "-Wpointer-arith"
  "-Wcast-qual" "-Wstrict-prototypes" "-Wconversion" "-Wsign-conversion"
  "-Wunused-macros"
  "$<$<NOT:$<STREQUAL:$<CONFIGURATION>,Release>>:-ggdb>"
  "$<$<NOT:$<STREQUAL:$<CONFIGURATION>,Release>>:-fno-omit-frame-pointer>"
  "$<$<STREQUAL:$<CONFIGURATION>,Release>:-DNDEBUG>"
  "$<$<STREQUAL:$<CONFIGURATION>,Release>:-fomit-frame-pointer>"
)

INCLUDE(CheckCCompilerFlag)
MACRO(FLAGS_MAGIC flags)
    CHECK_C_COMPILER_FLAG("-${flags}" SHAREMIND_HAS_CFLAGS__${flags})
    IF(SHAREMIND_HAS_CFLAGS__${flags})
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -${flags}")
    ENDIF(SHAREMIND_HAS_CFLAGS__${flags})
ENDMACRO(FLAGS_MAGIC)
FLAGS_MAGIC("Weverything")
FLAGS_MAGIC("Wlogical-op")

FIND_PACKAGE(SharemindLibsoftfloat REQUIRED)
FIND_PACKAGE(SharemindVmGlobals REQUIRED)

INCLUDE_DIRECTORIES(
    ${SharemindLibsoftfloat_INCLUDE_DIRS}
    ${SharemindVmGlobals_INCLUDE_DIRS}
)

ADD_DEFINITIONS(
    ${SharemindLibsoftfloat_DEFINITIONS}
    ${SharemindVmGlobals_DEFINITIONS}
    "-D__STDC_LIMIT_MACROS"
    "-D__STDC_CONSTANT_MACROS"
)

################################################################################
# libsoftfloat_math
################################################################################

SET(SHAREMIND_LIBSOFTFLOAT_MATH_HEADERS
    "3rdparty/libsoftfloat_math/sf_abs.h"
    "3rdparty/libsoftfloat_math/sf_ceil.h"
    "3rdparty/libsoftfloat_math/sf_erf.h"
    "3rdparty/libsoftfloat_math/sf_exp.h"
    "3rdparty/libsoftfloat_math/sf_floor.h"
    "3rdparty/libsoftfloat_math/sf_log.h"
    "3rdparty/libsoftfloat_math/sf_scalbn.h"
    "3rdparty/libsoftfloat_math/sf_sin.h"
)

SET(SHAREMIND_LIBSOFTFLOAT_MATH_SOURCES
    "3rdparty/libsoftfloat_math/internal/__sf_cos.c"
    "3rdparty/libsoftfloat_math/internal/__sf_rem_pio2.c"
    "3rdparty/libsoftfloat_math/internal/__sf_rem_pio2_large.c"
    "3rdparty/libsoftfloat_math/internal/__sf_sin.c"
    "3rdparty/libsoftfloat_math/sf_abs.c"
    "3rdparty/libsoftfloat_math/sf_ceil.c"
    "3rdparty/libsoftfloat_math/sf_erf.c"
    "3rdparty/libsoftfloat_math/sf_exp.c"
    "3rdparty/libsoftfloat_math/sf_floor.c"
    "3rdparty/libsoftfloat_math/sf_log.c"
    "3rdparty/libsoftfloat_math/sf_scalbn.c"
    "3rdparty/libsoftfloat_math/sf_sin.c"
)

ADD_LIBRARY("libsoftfloat_math" STATIC ${SHAREMIND_LIBSOFTFLOAT_MATH_SOURCES}
                                       ${SHAREMIND_LIBSOFTFLOAT_MATH_HEADERS})

TARGET_LINK_LIBRARIES("libsoftfloat_math"
    ${SharemindLibsoftfloat_LIBRARIES}
    ${SharemindVmGlobals_LIBRARIES}
)

SET(LIBSOFTFLOAT_MATH_VERSION "0.1.0.0")
SET_TARGET_PROPERTIES("libsoftfloat_math" PROPERTIES
                      OUTPUT_NAME "sharemind_softfloat_math"
                      VERSION "${LIBSOFTFLOAT_MATH_VERSION}"
                      SOVERSION "0.1")

INSTALL(TARGETS "libsoftfloat_math"
        ARCHIVE DESTINATION "lib")

INSTALL(FILES ${SHAREMIND_LIBSOFTFLOAT_MATH_HEADERS}
        DESTINATION "include/sharemind/3rdparty/libsoftfloat_math/")

# Install cmake files:
FILE(GLOB_RECURSE "sharemind_libsoftfloat_math_in_files"
     "${CMAKE_CURRENT_SOURCE_DIR}/*.in")
ADD_CUSTOM_TARGET("include_sharemind_libsoftfloat_math_in_files_in_IDE" SOURCES
                  ${sharemind_libsoftfloat_math_in_files})
SET(LIBSOFTFLOAT_MATH_INSTALL_LINK_LIBRARIES
    "-L${CMAKE_INSTALL_PREFIX}/lib"
    "-lsharemind_softfloat_math"
)

SET(LIBSOFTFLOAT_MATH_INSTALL_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")

CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/SharemindLibsoftfloatMathConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/SharemindLibsoftfloatMathConfig.cmake" @ONLY)
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/SharemindLibsoftfloatMathConfigVersion.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/SharemindLibsoftfloatMathConfigVersion.cmake"
    @ONLY)
INSTALL(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SharemindLibsoftfloatMathConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SharemindLibsoftfloatMathConfigVersion.cmake"
    DESTINATION "lib/SharemindLibsoftfloatMath")

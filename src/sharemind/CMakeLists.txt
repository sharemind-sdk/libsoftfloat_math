#
# Copyright (c) 2015, Cybernetica AS
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   * Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#
#   * Neither the name of Cybernetica AS nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL CYBERNETICA AS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

CMAKE_MINIMUM_REQUIRED(VERSION "2.8.12")
PROJECT("SHAREMIND_LIBSOFTFLOAT_MATH" "C")

SET(LIBSOFTFLOAT_MATH_SOVERSION "0.1")
SET(LIBSOFTFLOAT_MATH_VERSION "${LIBSOFTFLOAT_MATH_SOVERSION}.0")

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
FLAGS_MAGIC("Wno-padded")

FIND_PACKAGE(SharemindCHeaders REQUIRED)
FIND_PACKAGE(SharemindLibsoftfloat REQUIRED)

INCLUDE_DIRECTORIES(
    ${SharemindCHeaders_INCLUDE_DIRS}
    ${SharemindLibsoftfloat_INCLUDE_DIRS}
)

ADD_DEFINITIONS(
    ${SharemindCHeaders_DEFINITIONS}
    ${SharemindLibsoftfloat_DEFINITIONS}
    "-D__STDC_LIMIT_MACROS"
    "-D__STDC_CONSTANT_MACROS"
)

################################################################################
# libsoftfloat_math
################################################################################

FILE(GLOB SHAREMIND_LIBSOFTFLOAT_MATH_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libsoftfloat_math/*.h")
FILE(GLOB_RECURSE SHAREMIND_LIBSOFTFLOAT_MATH_HEADERS_P
     "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libsoftfloat_math/*.h")
LIST(REMOVE_ITEM SHAREMIND_LIBSOFTFLOAT_MATH_HEADERS_P
     ${SHAREMIND_LIBSOFTFLOAT_MATH_HEADERS})
FILE(GLOB_RECURSE SHAREMIND_LIBSOFTFLOAT_MATH_SOURCES
     "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libsoftfloat_math/*.c")

ADD_LIBRARY("libsoftfloat_math" STATIC ${SHAREMIND_LIBSOFTFLOAT_MATH_SOURCES}
                                       ${SHAREMIND_LIBSOFTFLOAT_MATH_HEADERS}
                                       ${SHAREMIND_LIBSOFTFLOAT_MATH_HEADERS_P})

TARGET_LINK_LIBRARIES("libsoftfloat_math"
    ${SharemindCHeaders_LIBRARIES}
    ${SharemindLibsoftfloat_LIBRARIES}
)

SET_TARGET_PROPERTIES("libsoftfloat_math" PROPERTIES
                      OUTPUT_NAME "sharemind_softfloat_math"
                      VERSION "${LIBSOFTFLOAT_MATH_VERSION}"
                      SOVERSION "${LIBSOFTFLOAT_MATH_SOVERSION}")

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

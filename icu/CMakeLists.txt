cmake_minimum_required(VERSION 3.21)
project(icu-superbuild LANGUAGES C CXX)

include(ExternalProject)

set(ICU_SRC_DIR ${CMAKE_SOURCE_DIR}/icu/src/icu4c)
set(ICU_HOST_BUILD ${CMAKE_BINARY_DIR}/icu-host)
set(ICU_TARGET_BUILD ${CMAKE_BINARY_DIR}/icu-target)
set(ICU_INSTALL_DIR ${CMAKE_BINARY_DIR}/icu-install)

set(CMAKE_C_FLAGS "-Dchar16_t=uint16_t")

if(MSVC)
  set(CMAKE_CXX_FLAGS "/utf-8")
endif()

# -------------------------------------------------------
# Stage 1: Host ICU (build tools only)
# -------------------------------------------------------

ExternalProject_Add(icu-host
    SOURCE_DIR ${ICU_SRC_DIR}
    CONFIGURE_COMMAND
        <SOURCE_DIR>/source/configure
        --disable-shared
        --enable-static
        --disable-tests
    BUILD_COMMAND make -j
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
)

# -------------------------------------------------------
# Stage 2: Target ICU (static libraries)
# -------------------------------------------------------

# Detect cross-compiling
set(ICU_HOST_ARG "")
set(ICU_ENV "")

if(CMAKE_CROSSCOMPILING)
    set(ICU_HOST_ARG --with-cross-build=${ICU_HOST_BUILD})

    set(ICU_ENV
        CC=${CMAKE_C_COMPILER}
        CXX=${CMAKE_CXX_COMPILER}
    )
endif()

ExternalProject_Add(icu-target
    DEPENDS icu-host
    SOURCE_DIR ${ICU_SRC_DIR}

    CONFIGURE_COMMAND
        ${ICU_ENV}
        <SOURCE_DIR>/source/configure
        ${ICU_HOST_ARG}
        --enable-static
        --disable-shared
        --disable-tools
        --disable-tests
        --with-data-packaging=static
        --prefix=${ICU_INSTALL_DIR}

    BUILD_COMMAND make -j
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE 1
)

# -------------------------------------------------------
# Imported targets for consumers
# -------------------------------------------------------

add_library(icuuc STATIC IMPORTED GLOBAL)
add_library(icui18n STATIC IMPORTED GLOBAL)
add_library(icudata STATIC IMPORTED GLOBAL)

add_dependencies(icuuc icu-target)
add_dependencies(icui18n icu-target)
add_dependencies(icudata icu-target)

set_target_properties(icuuc PROPERTIES
    IMPORTED_LOCATION ${ICU_INSTALL_DIR}/lib/libicuuc.a
)

set_target_properties(icui18n PROPERTIES
    IMPORTED_LOCATION ${ICU_INSTALL_DIR}/lib/libicui18n.a
)

set_target_properties(icudata PROPERTIES
    IMPORTED_LOCATION ${ICU_INSTALL_DIR}/lib/libicudata.a
)

set(ICU_INCLUDE_DIR ${ICU_INSTALL_DIR}/include CACHE PATH "ICU include directory")

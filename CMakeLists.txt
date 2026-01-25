cmake_minimum_required(VERSION 3.21)

project(gpkg VERSION 0.9.18 LANGUAGES C CXX)

# Force build in proper directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "In-source builds are not supported.\n"
    "Remove CMakeCache.txt and CMakeFiles/, then run:\n"
    "  cmake -S . -B build"
  )
endif()


# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE true)
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_CURRENT_SOURCE_DIR}/c_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/cxx_overrides.cmake)

set(GPKG_BUILD_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
set(GPKG_INTERFACE_VERSION ${PROJECT_VERSION_MAJOR})


# Required CMake modules
include(CMakeDependentOption)
include(CheckIncludeFile)
include(CheckFunctionExists)

option(GPKG_TEST "Enable Testing?" OFF)
option(GPKG_COVERAGE "Enable Code coverage?" OFF)
option(GPKG_ICU "Enable International Components for Unicode (ICU)?" OFF)



add_subdirectory(gpkg)
add_subdirectory(shell)
add_subdirectory(sqlite)


find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(
    doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif()


# International Components for Unicode (ICU)
if(GPKG_ICU)
  add_compile_definitions(SQLITE_ENABLE_ICU=1)
  add_compile_definitions(SQLITE_ENABLE_ICU_COLLATIONS=1)
  add_subdirectory(icu)
endif()

if(GPKG_TEST)
  include(CTest)
  add_subdirectory(test)
endif()
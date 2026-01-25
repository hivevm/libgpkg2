# Required CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


#### SQLite
add_definitions(-D_HAVE_SQLITE_CONFIG_H=1)

check_include_file(stdint.h HAVE_STDINT_H)
check_include_file(inttype.h HAVE_INTTYPES_H)
check_include_file(utime.h HAVE_UTIME_H)

check_function_exists(gmtime_r HAVE_GMTIME_R)
check_function_exists(localtime_r HAVE_LOCALTIME_R)
check_function_exists(posix_fallocate HAVE_POSIX_FALLOCATE)
check_function_exists(usleep HAVE_USLEEP)
check_function_exists(fdatasync HAVE_FDATASYNC)
check_function_exists(fullfsync HAVE_FULLFSYNC)
check_function_exists(strerror_r HAVE_STRERROR_R)


set(SQLITE_ENABLE_RTREE 1)

set(SQLITE_ENABLE_FTS3 1)
set(SQLITE_ENABLE_FTS3_PARENTHESIS 1)
set(SQLITE_ENABLE_FTS5 1)

set(SQLITE_DISABLE_LFS 1)
set(SQLITE_ENABLE_COLUMN_METADATA 1)

# set(SQLITE_OMIT_LOAD_EXTENSION 1)
#### SQLite


# Common objects for static and shared library
set(GPKG_SOURCE_FILES
  binstream.c
  blobio.c
  error.c
  fp.c
  geomio.c
  gpkg.c
  gpkg_db.c
  gpkg_geom.c
  i18n.c
  sql.c
  spatialdb.c
  strbuf.c
  wkb.c
  wkt.c
)

#
# Check to determine how to do thread local storage
#
include(UseTLS)
check_tls()

include(UseLocale)
check_locale()

include_directories("${PROJECT_SOURCE_DIR}/sqlite")
list(APPEND GPKG_SOURCE_FILES 
  ${BOOST_GEOMETRY_SOURCE}
  ${PROJECT_SOURCE_DIR}/sqlite/sqlite3.c
)

# Geometry functions library from Boost.Geometry
file(GLOB_RECURSE BOOST_GEOMETRY_SOURCE
  "${PROJECT_SOURCE_DIR}/boost/boost/*"
)
include_directories("${PROJECT_SOURCE_DIR}/boost")


list(APPEND GPKG_SOURCE_FILES
  boost_geom_func.cpp
  boost_geom_io.cpp
  boost_geom_io.hpp
  boost_geom_proj.cpp
  boost_geom_proj.hpp
  boost_geom_wrapper.cpp
  boost_geom_wrapper.hpp
  boost_geometries.hpp
)

set(GPKG_GEOM_FUNC "Boost.Geometry")
set(GPKG_GEOM_FUNC ${GPKG_GEOM_FUNC} PARENT_SCOPE)

if (GPKG_GEOM_FUNC)
  list(APPEND GPKG_SOURCE_FILES "geom_func.h")
endif()


# Static and shared library version of geopackage.
add_library(gpkg_static STATIC
  ${GPKG_SOURCE_FILES}
)
set_target_properties(gpkg_static PROPERTIES
  OUTPUT_NAME "gpkg"
  LIBRARY_OUTPUT_DIRECTORY static
  COMPILE_DEFINITIONS "SQLITE_CORE=1"
)

if (NOT IOS)
  add_library(gpkg_shared SHARED
  ${GPKG_SOURCE_FILES}
)
  set_target_properties(gpkg_shared PROPERTIES
    OUTPUT_NAME "gpkg"
    LIBRARY_OUTPUT_DIRECTORY shared
  )
endif()


find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT)
    set(SQLITE_THREADSAFE 2)
    target_link_libraries(gpkg_static "${CMAKE_THREAD_LIBS_INIT}")
  if(NOT IOS)
    target_link_libraries(gpkg_shared "${CMAKE_THREAD_LIBS_INIT}")
  endif()
elseif(CMAKE_USE_WIN32_THREADS_INIT)
    set(SQLITE_THREADSAFE 2)
else()
    set(SQLITE_THREADSAFE 1)
endif()


if(UNIX)
  target_link_libraries(gpkg_static "${CMAKE_DL_LIBS}")
  if(NOT IOS)
    target_link_libraries(gpkg_shared "${CMAKE_DL_LIBS}")
  endif()
endif()


add_definitions(-DGPKG_HAVE_CONFIG_H)
configure_file("${PROJECT_SOURCE_DIR}/sqlite/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/sqlite_cfg.h")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Strip 'lib' prefix when compiling for windows
if (WIN32)
  set_target_properties(gpkg_static PROPERTIES PREFIX "")
  set_target_properties(gpkg_shared PROPERTIES PREFIX "")
endif()

if (MSVC)
  add_definitions(-DGPKG_CALL=__cdecl)
  set_target_properties(gpkg_shared PROPERTIES
    COMPILE_DEFINITIONS "SQLITE_API=__declspec\(dllexport\);GPKG_EXPORT=__declspec\(dllexport\)")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /TP")
endif()

# Compiler specific settings
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU" OR ${CMAKE_C_COMPILER_ID} MATCHES "Clang")
  add_definitions(-DGPKG_EXPORT=__attribute__\(\(visibility\("default"\)\)\))
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fvisibility=hidden -std=c99 -pedantic -Wall -Wno-unused-variable")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -std=c++11 -pedantic -Wall -Wno-unused-variable")
endif()


# Enable code coverage
if (${CMAKE_C_COMPILER_ID} MATCHES "GNU" AND GPKG_COVERAGE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
  set_target_properties(gpkg_shared PROPERTIES LINK_FLAGS "--coverage")
endif()


install(DIRECTORY ${PROJECT_SOURCE_DIR}/sqlite
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
)
install(DIRECTORY .
    DESTINATION include
    FILES_MATCHING
        PATTERN "spatialdb.h"
)


set_target_properties(gpkg_static PROPERTIES
    OUTPUT_NAME gpkg
)

install(TARGETS gpkg_static
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)


if(NOT IOS)
  set_target_properties(gpkg_shared PROPERTIES
      OUTPUT_NAME gpkg
      VERSION ${GPKG_BUILD_VERSION}
      SOVERSION ${GPKG_INTERFACE_VERSION}
  )

  install(TARGETS gpkg_shared
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib
  )
endif()

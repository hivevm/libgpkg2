# Required CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


# Geopackage shell command line tool
add_executable(GPKG_SHELL shell.c)

include_directories(
  "${PROJECT_SOURCE_DIR}/sqlite"
  "${PROJECT_SOURCE_DIR}/gpkg"
  "${PROJECT_BINARY_DIR}"
)

target_link_libraries(GPKG_SHELL
  #sqlite_static
  gpkg_static
)

find_package(Readline)
if (READLINE_FOUND)
  set (HAVE_READLINE=1)
  include_directories("${READLINE_INCLUDE_DIR}")
  target_link_libraries(GPKG_SHELL "${READLINE_LIBRARY}")
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h")
add_definitions(-DHAVE_CONFIG_H=1)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set_target_properties(
  GPKG_SHELL PROPERTIES
  OUTPUT_NAME "geopackage"
)

if(${CMAKE_C_COMPILER_ID} MATCHES "GNU" AND GPKG_COVERAGE)
  set_target_properties(GPKG_SHELL PROPERTIES LINK_FLAGS "--coverage")
endif()


# Installation rules
install(TARGETS GPKG_SHELL
    RUNTIME DESTINATION bin
)
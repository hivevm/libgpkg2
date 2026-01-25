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


add_subdirectory(sqlite)
add_subdirectory(gpkg)
if(NOT ANDROID AND NOT IOS) 
    add_subdirectory(shell)
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

find_package(Doxygen)
    if(DOXYGEN_FOUND)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
        add_custom_target(
            doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif()

# ####################################################

# # Library source files
# set(LIBRARY_SOURCES
#     sqlite/sqlite3.c
# )

# # Library header files
# set(LIBRARY_HEADERS
#     sqlite/sqlite3.h
#     sqlite/sqlite3ext.h
# )

# # Create static library
# add_library(${PROJECT_NAME} STATIC
#     ${LIBRARY_SOURCES}
#     ${LIBRARY_HEADERS}
# )

# # install(
# #     DIRECTORY include/
# #     DESTINATION include
# #     FILES_MATCHING
# #         REGEX "^boost/.*\\.(h|hpp)$"
# # )

# # Set the public include directories for build and install
# target_include_directories(${PROJECT_NAME}
#     PUBLIC
#         $<BUILD_INTERFACE:
#             ${CMAKE_CURRENT_SOURCE_DIR}/boost
#             ${CMAKE_CURRENT_SOURCE_DIR}/include
#         >
#         $<INSTALL_INTERFACE:include>
# )

# # C/C++ Version definitions
# target_compile_definitions(gpkg
#   PRIVATE
#     GPKG_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
#     GPKG_VERSION_MINOR=${PROJECT_VERSION_MINOR}
#     GPKG_VERSION_PATCH=${PROJECT_VERSION_PATCH}
# )

# # Platform-specific settings
# if(ANDROID)
#     target_compile_definitions(${PROJECT_NAME} PRIVATE PLATFORM_ANDROID)
# elseif(IOS)
#     target_compile_definitions(${PROJECT_NAME} PRIVATE PLATFORM_IOS)
#     # Remove FRAMEWORK settings - we want a static library, not a framework
#     # set_target_properties(${PROJECT_NAME} PROPERTIES
#     #     FRAMEWORK TRUE
#     #     MACOSX_FRAMEWORK_IDENTIFIER org.hivevm.mylibrary
#     # )
# elseif(WIN32)
#     target_compile_definitions(${PROJECT_NAME} PRIVATE PLATFORM_WINDOWS)
# elseif(APPLE)
#     target_compile_definitions(${PROJECT_NAME} PRIVATE PLATFORM_MACOS)
# elseif(UNIX)
#     target_compile_definitions(${PROJECT_NAME} PRIVATE PLATFORM_LINUX)
# endif()


# # Installation rules
# install(TARGETS ${PROJECT_NAME}
#     ARCHIVE DESTINATION lib
#     LIBRARY DESTINATION lib
#     FRAMEWORK DESTINATION lib
# )
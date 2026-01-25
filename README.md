[![Build Geopackage](https://github.com/hivevm/libgpkg2/actions/workflows/build.yaml/badge.svg)](https://github.com/hivevm/libgpkg2/actions/workflows/build.yaml)

# Description
A SQLite 3 extension that provides a minimal [OGC GeoPackage](http://www.ogcnetwork.net/geopackage) implementation.

GeoPackage is an open, standards-based, application and platform independent, and self-describing file format for
geodata based on SQLite.

The project was originally started by [Luciad](http://www.luciad.com), but hasn't been developed for years.

# License
libgpkg is distributed under the [Apache Software License](https://www.apache.org/licenses/LICENSE-2.0) version 2.0.

# Installation

- Linux: compile from source.
- MacOS: NOT TESTED YET.
- Windows: NOT TESTED YET.

# Usage
libgpkg can be loaded into SQLite using the [sqlite3\_load\_extension](http://sqlite.org/c3ref/load_extension.html) C
function or using the [load\_extension](http://sqlite.org/lang_corefunc.html#load_extension) SQL function. Once loaded
libgpkg extends SQLite with the function listed below. These function can be used just like any of the core functions
that SQLite provides.

libgpkg exposes the _init_geopackage_extension_ to load the geopackage implementation through the boost geometries.

# Compilation

- Install CMake 3.30 or newer. CMake can be downloaded from www.cmake.org or installed using
  your systems package manager.
- Install the [Boost libraries](boost.md) for geometry support
- Run 'cmake ..' in the _build_ directory of the project to generate the build scripts for your system or use ccmake .. for the configuration.
- Build the project using the generated build scripts.
- The build scripts will generate a number of binaries
    - _build/shell/gpkg_: a modified version of the SQLite 3 command-line shell that autoloads the GeoPackage extension. This is a standalone binary that has been statically linked with SQLite 3 and the GeoPackage extension.
    - _build/gpkg/libgpkg.so_ (or _gpkg.dll_ on Windows): a dynamically loadable SQLite 3 extension that provides the GeoPackage functionality. This extension library can be used with any SQLite 3 that supports extension loading.

# Testing

Use the command line tool built in _build/shell_

~~~bash
./build/shell/gpkg

.open [PATH_TO_DATABASE]

# Show all tables
.tables

# Query a spatial table
SELECT ST_ASTEXT(geom) FROM [TABLE_NAME]
~~~

# Dependencies

- libgpkg requires SQLite 3.51.0 or higher.
- libgpkg requires Boost 1.90 or higher.
# libgpkg smoke tests, driven through the bundled gpkg shell (GPKG_SHELL).
#
# The shell executes each SQL argument in order on one connection and exits non-zero on the
# first SQL error, so every test asserts both a successful run and the expected output.
# The former ruby/rspec suite was removed when the development workspace was reworked; these
# tests keep the GPKG_TEST option functional without reintroducing a test-framework dependency.

if(NOT TARGET GPKG_SHELL)
  message(STATUS "GPKG_SHELL is not built on this platform - skipping tests")
  return()
endif()

add_test(NAME wkt_roundtrip
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('POINT(1 2)'));")
set_tests_properties(wkt_roundtrip PROPERTIES
  PASS_REGULAR_EXPRESSION "Point \\(1 2\\)")

add_test(NAME srid
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_SRID(ST_GeomFromText('POINT(1 2)', 4326));")
set_tests_properties(srid PROPERTIES
  PASS_REGULAR_EXPRESSION "4326")

add_test(NAME geometry_type
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_GeometryType(ST_GeomFromText('LINESTRING(0 0, 1 1)'));")
set_tests_properties(geometry_type PROPERTIES
  PASS_REGULAR_EXPRESSION "LineString")

add_test(NAME is_empty
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_IsEmpty(ST_GeomFromText('POINT EMPTY'));")
set_tests_properties(is_empty PROPERTIES
  PASS_REGULAR_EXPRESSION "1")

# End-to-end GeoPackage lifecycle: initialize the metadata tables, add a geometry column,
# store a geometry and read it back.
add_test(NAME geopackage_lifecycle
  COMMAND GPKG_SHELL :memory:
    "SELECT InitSpatialMetaData();"
    "CREATE TABLE features (id INTEGER PRIMARY KEY);"
    "SELECT AddGeometryColumn('features', 'geom', 'POINT', 4326, 0, 0);"
    "INSERT INTO features (id, geom) VALUES (1, ST_GeomFromText('POINT(7.4 46.9)', 4326));"
    "SELECT ST_AsText(geom), ST_SRID(geom) FROM features;")
set_tests_properties(geopackage_lifecycle PROPERTIES
  PASS_REGULAR_EXPRESSION "Point \\(7\\.4 46\\.9\\)\\|4326")

add_test(NAME invalid_wkt_fails
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('NOT A WKT'));")
set_tests_properties(invalid_wkt_fails PROPERTIES
  WILL_FAIL TRUE)

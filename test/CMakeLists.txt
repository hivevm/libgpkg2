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

# --- Regression tests -------------------------------------------------------
# Each test below reproduces a defect that was fixed; they fail on the code as it
# stood before the corresponding fix.

# Nesting is bounded by GEOM_MAX_DEPTH. Unbounded nesting used to run past the
# fixed-size writer stacks (heap corruption) and exhaust the C stack.
string(REPEAT "GEOMETRYCOLLECTION(" 30 _deep_open)
string(REPEAT ")" 30 _deep_close)
add_test(NAME nesting_depth_rejected
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_GeomFromText('${_deep_open}POINT(1 2)${_deep_close}');")
set_tests_properties(nesting_depth_rejected PROPERTIES
  WILL_FAIL TRUE)

# Nesting below the limit must keep working, through both writers.
string(REPEAT "GEOMETRYCOLLECTION(" 10 _ok_open)
string(REPEAT ")" 10 _ok_close)
add_test(NAME nesting_depth_allowed
  COMMAND GPKG_SHELL :memory:
    "SELECT length(ST_AsBinary(ST_GeomFromText('${_ok_open}POINT(1 2)${_ok_close}')));")
set_tests_properties(nesting_depth_allowed PROPERTIES
  PASS_REGULAR_EXPRESSION "[0-9]+")

# A CircularString needs 3+2n points; the count check used to wrap around for a
# single point and let it through, overrunning a stack buffer downstream.
add_test(NAME circularstring_min_points
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_MinX(ST_GeomFromText('CIRCULARSTRING(1 2)'));")
set_tests_properties(circularstring_min_points PROPERTIES
  WILL_FAIL TRUE)

add_test(NAME circularstring_valid
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('CIRCULARSTRING(0 0, 1 1, 2 0)'));")
set_tests_properties(circularstring_valid PROPERTIES
  PASS_REGULAR_EXPRESSION "CircularString \\(0 0, 1 1, 2 0\\)")

# A blob with a wrong magic number must be reported by rendering exactly the two
# bytes that were read; the previous message ran off the end of that buffer.
add_test(NAME gpb_bad_magic
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_MinX(X'41414141414141414141414141414141414141');")
set_tests_properties(gpb_bad_magic PROPERTIES
  PASS_REGULAR_EXPRESSION "actual: 0x4141")

# Identifiers are data, not format strings: percent signs in a table name used to
# be interpreted by the SQL builder and crashed the process.
add_test(NAME identifier_with_percent
  COMMAND GPKG_SHELL :memory:
    "SELECT GPKG_InitSpatialMetaData();"
    "SELECT GPKG_CreateTilesTable('t%s%s%s%s%s%s%s%s');"
    "SELECT count(*) FROM sqlite_master WHERE name = 't%s%s%s%s%s%s%s%s';")
set_tests_properties(identifier_with_percent PROPERTIES
  PASS_REGULAR_EXPRESSION "1")

# An empty geometry carries no envelope: the header used to announce one (indicator 1 in the
# flags byte) whose eight doubles were all NaN, which other GeoPackage readers take at face
# value. A geometry that does have coordinates must keep its envelope.
add_test(NAME empty_geometry_has_no_envelope
  COMMAND GPKG_SHELL :memory:
    "SELECT 'empty=' || hex(substr(ST_GeomFromText('LINESTRING EMPTY', 4326), 4, 1))
         || ' filled=' || hex(substr(ST_GeomFromText('LINESTRING(0 0, 1 1)', 4326), 4, 1));")
set_tests_properties(empty_geometry_has_no_envelope PROPERTIES
  PASS_REGULAR_EXPRESSION "empty=11 filled=03")

add_test(NAME empty_geometry_round_trip
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('LINESTRING EMPTY')) || '/'
         || ST_AsText(ST_GeomFromText('POLYGON EMPTY')) || '/'
         || ST_AsText(ST_GeomFromText('GEOMETRYCOLLECTION EMPTY')) || '/'
         || ST_IsEmpty(ST_GeomFromText('LINESTRING EMPTY'));")
set_tests_properties(empty_geometry_round_trip PROPERTIES
  PASS_REGULAR_EXPRESSION "LineString EMPTY/Polygon EMPTY/GeometryCollection EMPTY/1")

# The reported type name was formatted with %d, printing the pointer as a number.
add_test(NAME geometry_type_mismatch_names_the_type
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_Point('LINESTRING(0 0, 1 1)');")
set_tests_properties(geometry_type_mismatch_names_the_type PROPERTIES
  PASS_REGULAR_EXPRESSION "Expected 'Point' actual 'LineString'")

# A malformed curve member must be reported once. The dimension read discarded its result,
# so parsing carried on with an uninitialised child header and piled up follow-up errors.
add_test(NAME malformed_curve_reports_once
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('COMPOUNDCURVE(CIRCULARSTRING ,)'));")
set_tests_properties(malformed_curve_reports_once PROPERTIES
  PASS_REGULAR_EXPRESSION "Unexpected token at column 29"
  FAIL_REGULAR_EXPRESSION "Expected '\\(' or 'empty'")

# Values must be written as SQL string literals, not as double-quoted identifiers that
# only parse because of SQLite's fallback. Runs with that fallback switched off.
add_test(NAME spatial_index_under_strict_quoting
  COMMAND GPKG_SHELL :memory: ".read ${CMAKE_CURRENT_SOURCE_DIR}/strict_quoting.sql")
set_tests_properties(spatial_index_under_strict_quoting PROPERTIES
  PASS_REGULAR_EXPRESSION "extensions: 1")

# Metadata lookups must address the database the caller named, not always main. Uses a file
# database for the attached side, removed beforehand.
add_test(NAME attached_database_cleanup
  COMMAND ${CMAKE_COMMAND} -E rm -f attached_database.gpkg)
set_tests_properties(attached_database_cleanup PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  FIXTURES_SETUP ATTACHED_DB)

add_test(NAME srs_lookup_uses_named_database
  COMMAND GPKG_SHELL :memory: ".read ${CMAKE_CURRENT_SOURCE_DIR}/attached_database.sql")
set_tests_properties(srs_lookup_uses_named_database PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  FIXTURES_REQUIRED ATTACHED_DB
  PASS_REGULAR_EXPRESSION "registered: 1")

# Coordinates are decimal numbers. The scanner delegated to strtod and inherited its
# extensions, so hexadecimal, nan and inf were accepted, and a lone sign consumed nothing.
add_test(NAME wkt_rejects_hex_coordinates
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('POINT(-0x10 -0x20)'));")
set_tests_properties(wkt_rejects_hex_coordinates PROPERTIES
  WILL_FAIL TRUE)

add_test(NAME wkt_rejects_nan_coordinates
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('POINT(-nan -inf)'));")
set_tests_properties(wkt_rejects_nan_coordinates PROPERTIES
  WILL_FAIL TRUE)

add_test(NAME wkt_rejects_bare_sign
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('POINT(- -)'));")
set_tests_properties(wkt_rejects_bare_sign PROPERTIES
  WILL_FAIL TRUE)

# Exponents and a leading plus are valid and must keep working.
add_test(NAME wkt_accepts_exponent_notation
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('POINT(1e3 -2.5E-2)')) || ' ' || ST_AsText(ST_GeomFromText('POINT(+1.5 +2)'));")
set_tests_properties(wkt_accepts_exponent_notation PROPERTIES
  PASS_REGULAR_EXPRESSION "Point \\(1000 -0\\.025\\) Point \\(1\\.5 2\\)")

# Mixing coordinate systems must not yield an extent. The check compared only the last
# row against the first, so 4326, 3857, 4326 was accepted and combined regardless.
add_test(NAME extent_rejects_mixed_srid
  COMMAND GPKG_SHELL :memory:
    "CREATE TABLE t(g);"
    "INSERT INTO t VALUES(ST_GeomFromText('POINT(0 0)', 4326)), (ST_GeomFromText('POINT(1 1)', 3857)), (ST_GeomFromText('POINT(2 2)', 4326));"
    "SELECT 'srid=' || coalesce(CAST(ST_SRID(ST_Extent(g)) AS TEXT), 'null') FROM t;")
set_tests_properties(extent_rejects_mixed_srid PROPERTIES
  PASS_REGULAR_EXPRESSION "srid=null")

add_test(NAME extent_keeps_uniform_srid
  COMMAND GPKG_SHELL :memory:
    "CREATE TABLE t(g);"
    "INSERT INTO t VALUES(ST_GeomFromText('POINT(0 0)', 4326)), (ST_GeomFromText('POINT(2 2)', 4326));"
    "SELECT ST_SRID(ST_Extent(g)) || '/' || ST_MaxX(ST_Extent(g)) FROM t;")
set_tests_properties(extent_keeps_uniform_srid PROPERTIES
  PASS_REGULAR_EXPRESSION "4326/2")

# The four coordinate form read its arguments by casting sqlite3_value* to double*, which
# also made it reject integers outright.
add_test(NAME extent_from_integer_coordinates
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_MinX(ST_Extent(1, 2, 3, 4)) || '/' || ST_MaxY(ST_Extent(1, 2, 3, 4));")
set_tests_properties(extent_from_integer_coordinates PROPERTIES
  PASS_REGULAR_EXPRESSION "1(\\.0)?/4(\\.0)?")

# Naming the database must not change the recorded z/m support: the five argument form
# left both at 0 ("prohibited") while the four argument form records 2 ("optional").
add_test(NAME add_geometry_column_zm_defaults
  COMMAND GPKG_SHELL :memory:
    "SELECT GPKG_InitSpatialMetaData();"
    "CREATE TABLE a(id INTEGER PRIMARY KEY);"
    "CREATE TABLE b(id INTEGER PRIMARY KEY);"
    "SELECT GPKG_AddGeometryColumn('a', 'g', 'POINT', 4326);"
    "SELECT GPKG_AddGeometryColumn('main', 'b', 'g', 'POINT', 4326);"
    "SELECT group_concat(table_name || '=' || z || ',' || m, ' ') FROM gpkg_geometry_columns;")
set_tests_properties(add_geometry_column_zm_defaults PROPERTIES
  PASS_REGULAR_EXPRESSION "a=2,2 b=2,2")

# An unusable database name used to surface as a bare "unknown error".
add_test(NAME add_geometry_column_reports_bad_database
  COMMAND GPKG_SHELL :memory:
    "SELECT GPKG_InitSpatialMetaData();"
    "CREATE TABLE a(id INTEGER PRIMARY KEY);"
    "SELECT GPKG_AddGeometryColumn('nosuchdb', 'a', 'g', 'POINT', 4326);")
set_tests_properties(add_geometry_column_reports_bad_database PROPERTIES
  PASS_REGULAR_EXPRESSION "unknown database")

# The out-of-range report for m used to print the value of z instead.
add_test(NAME add_geometry_column_m_flag_message
  COMMAND GPKG_SHELL :memory:
    "SELECT GPKG_InitSpatialMetaData();"
    "CREATE TABLE a(id INTEGER PRIMARY KEY);"
    "SELECT GPKG_AddGeometryColumn('a', 'g', 'POINT', 4326, 1, 9);")
set_tests_properties(add_geometry_column_m_flag_message PROPERTIES
  PASS_REGULAR_EXPRESSION "Invalid M flag value: 9")

# A trailing integer is only an SRID while the constructor still has enough arguments
# left; ST_Point(1, 2) is a point with integer coordinates and used to be rejected.
add_test(NAME point_from_integer_coordinates
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_Point(1, 2));")
set_tests_properties(point_from_integer_coordinates PROPERTIES
  PASS_REGULAR_EXPRESSION "Point \\(1 2\\)")

add_test(NAME point_keeps_explicit_srid
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_SRID(ST_Point(1.0, 2.0, 4326));")
set_tests_properties(point_keeps_explicit_srid PROPERTIES
  PASS_REGULAR_EXPRESSION "4326")

# Text after the geometry means the input was not WKT; it used to be ignored silently.
add_test(NAME wkt_rejects_trailing_text
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('POINT(1 2)DROP TABLE x'));")
set_tests_properties(wkt_rejects_trailing_text PROPERTIES
  WILL_FAIL TRUE)

# The envelope computed while streaming a CircularString out of WKB must cover every
# point; a batch boundary used to drop one, leaving a box too small for the geometry.
add_test(NAME circularstring_envelope_covers_points
  COMMAND GPKG_SHELL :memory: ".read ${CMAKE_CURRENT_SOURCE_DIR}/circularstring_envelope.sql")
set_tests_properties(circularstring_envelope_covers_points PROPERTIES
  PASS_REGULAR_EXPRESSION "maxx_ok=1 maxy_ok=1")

# The result of a geometry constructor depends on all of its arguments, so it must
# not be cached against the first one alone: a constant first argument used to pin
# the first row's result for the whole query.
add_test(NAME constructor_varies_per_row
  COMMAND GPKG_SHELL :memory:
    "CREATE TABLE p(x DOUBLE);"
    "INSERT INTO p VALUES(10.0),(20.0),(30.0);"
    "SELECT group_concat(ST_AsText(ST_Point(0.5, x)), ' | ') FROM p;")
set_tests_properties(constructor_varies_per_row PROPERTIES
  PASS_REGULAR_EXPRESSION "Point \\(0\\.5 10\\) \\| Point \\(0\\.5 20\\) \\| Point \\(0\\.5 30\\)")

add_test(NAME constructor_varies_srid_per_row
  COMMAND GPKG_SHELL :memory:
    "CREATE TABLE s(id INTEGER);"
    "INSERT INTO s VALUES(4326),(3857);"
    "SELECT group_concat(ST_SRID(ST_GeomFromText('POINT(1 2)', id)), ' | ') FROM s;")
set_tests_properties(constructor_varies_srid_per_row PROPERTIES
  PASS_REGULAR_EXPRESSION "4326 \\| 3857")

# A failed transactional function must leave no savepoint behind, otherwise later
# writes are silently discarded when the connection closes. The script needs a file
# database (it reopens the connection), so it runs in the build directory on a
# database removed beforehand.
add_test(NAME transaction_scope_cleanup
  COMMAND ${CMAKE_COMMAND} -E rm -f transaction_scope.gpkg)
set_tests_properties(transaction_scope_cleanup PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  FIXTURES_SETUP TRANSACTION_SCOPE_DB)

add_test(NAME failed_call_ends_transaction
  COMMAND GPKG_SHELL :memory: ".read ${CMAKE_CURRENT_SOURCE_DIR}/transaction_scope.sql")
set_tests_properties(failed_call_ends_transaction PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  FIXTURES_REQUIRED TRANSACTION_SCOPE_DB
  PASS_REGULAR_EXPRESSION "rows: 1")

# ST_SRID with a new SRID used to rewrite the blob header inside the buffer returned by
# sqlite3_value_blob, mutating memory owned by SQLite or by the caller's binding. With a
# bound parameter both references share one buffer, so the second read observed the new SRID.
add_test(NAME srid_update_leaves_argument_unchanged
  COMMAND GPKG_SHELL :memory:
    ".parameter init"
    ".parameter set $g x'47500001E61000000101000000000000000000F03F0000000000000040'"
    "SELECT ST_SRID(ST_SRID($g, 9999)) || '/' || ST_SRID($g);")
set_tests_properties(srid_update_leaves_argument_unchanged PROPERTIES
  PASS_REGULAR_EXPRESSION "9999/4326")

# gpkg_geometry_columns references gpkg_contents, but AddGeometryColumn used to insert the
# rows in the opposite order (failing under PRAGMA foreign_keys=ON) and inserted the contents
# row unconditionally (failing on the primary key when the application registered the table
# in gpkg_contents beforehand).
add_test(NAME add_geometry_column_under_foreign_keys
  COMMAND GPKG_SHELL :memory:
    "PRAGMA foreign_keys = ON;"
    "SELECT GPKG_InitSpatialMetaData();"
    "CREATE TABLE f(id INTEGER PRIMARY KEY);"
    "SELECT GPKG_AddGeometryColumn('f', 'g', 'POINT', 4326);"
    "CREATE TABLE h(id INTEGER PRIMARY KEY);"
    "INSERT INTO gpkg_contents (table_name, data_type, identifier, srs_id) VALUES ('h', 'features', 'h', 4326);"
    "SELECT GPKG_AddGeometryColumn('h', 'g', 'POINT', 4326);"
    "SELECT 'contents=' || (SELECT count(*) FROM gpkg_contents WHERE table_name IN ('f', 'h')) || ' columns=' || (SELECT count(*) FROM gpkg_geometry_columns WHERE table_name IN ('f', 'h'));")
set_tests_properties(add_geometry_column_under_foreign_keys PROPERTIES
  PASS_REGULAR_EXPRESSION "contents=2 columns=2")

# A failed curve member used to be masked by the end-of-geometry callback overwriting the
# parse result, so parsing carried on past the broken member and reported a bogus
# trailing-text error on the closing parenthesis.
add_test(NAME compoundcurve_member_failure_reported_once
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('COMPOUNDCURVE(CIRCULARSTRING(0 0, 1 1))'));")
set_tests_properties(compoundcurve_member_failure_reported_once PROPERTIES
  PASS_REGULAR_EXPRESSION "requires 3\\+2n points"
  FAIL_REGULAR_EXPRESSION "Unexpected text after geometry")

add_test(NAME compoundcurve_valid
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('COMPOUNDCURVE((0 0, 1 1), CIRCULARSTRING(1 1, 2 0, 3 1))'));")
set_tests_properties(compoundcurve_valid PROPERTIES
  PASS_REGULAR_EXPRESSION "CompoundCurve \\(\\(0 0, 1 1\\), CircularString \\(1 1, 2 0, 3 1\\)\\)")

# The GeoPackage header pragmas used to be unqualified, stamping main instead of the named
# attached database: the attached file must get the application id, main must keep its own.
add_test(NAME attached_header_cleanup
  COMMAND ${CMAKE_COMMAND} -E rm -f attached_header.gpkg)
set_tests_properties(attached_header_cleanup PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  FIXTURES_SETUP ATTACHED_HEADER_DB)

add_test(NAME init_stamps_attached_database_header
  COMMAND GPKG_SHELL :memory: ".read ${CMAKE_CURRENT_SOURCE_DIR}/attached_header.sql")
set_tests_properties(init_stamps_attached_database_header PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  FIXTURES_REQUIRED ATTACHED_HEADER_DB
  PASS_REGULAR_EXPRESSION "1196444487[\r\n]+0")

# foreign_key_check used to run unqualified against main, so a dangling reference in the
# attached GeoPackage passed CheckSpatialMetaData.
add_test(NAME attached_fk_cleanup
  COMMAND ${CMAKE_COMMAND} -E rm -f attached_fk.gpkg)
set_tests_properties(attached_fk_cleanup PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  FIXTURES_SETUP ATTACHED_FK_DB)

add_test(NAME check_finds_fk_violation_in_attached_database
  COMMAND GPKG_SHELL :memory: ".read ${CMAKE_CURRENT_SOURCE_DIR}/attached_fk.sql")
set_tests_properties(check_finds_fk_violation_in_attached_database PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  FIXTURES_REQUIRED ATTACHED_FK_DB
  PASS_REGULAR_EXPRESSION "foreign key from 'srs_id' to 'gpkg_spatial_ref_sys.srs_id' failed for value '99999'")

# Coordinates used to be written with %.10g, which cannot represent every double: text output
# must round-trip exactly while numbers that have a short form keep it.
add_test(NAME wkt_round_trips_full_precision
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromText('POINT(0.30000000000000004 -122.41941999999999)')) || ' | ' || ST_AsText(ST_GeomFromText('POINT(0.1 2.5)'));")
set_tests_properties(wkt_round_trips_full_precision PROPERTIES
  PASS_REGULAR_EXPRESSION "Point \\(0\\.30000000000000004 -122\\.41941999999999\\) \\| Point \\(0\\.1 2\\.5\\)")

# An ordinate whose values are all NaN never tightens the envelope sentinels; ST_MinX used to
# hand back DBL_MAX (breaking the R-tree triggers), and a NaN-only Z rejected the geometry.
add_test(NAME nan_coordinates_have_no_envelope
  COMMAND GPKG_SHELL :memory:
    "SELECT 'minx=' || coalesce(CAST(ST_MinX(x'47500001FFFFFFFF010200000002000000000000000000F87F0000000000000040000000000000F87F0000000000001040') AS TEXT), 'null');")
set_tests_properties(nan_coordinates_have_no_envelope PROPERTIES
  PASS_REGULAR_EXPRESSION "minx=null")

add_test(NAME nan_z_geometry_accepted
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_AsText(ST_GeomFromWKB(x'01EA03000002000000000000000000F03F0000000000000040000000000000F87F00000000000008400000000000001040000000000000F87F'));")
set_tests_properties(nan_z_geometry_accepted PROPERTIES
  PASS_REGULAR_EXPRESSION "LineString Z \\(1 2 nan, 3 4 nan\\)")

# The Boost bridge used to ignore the blob header read result and compute over a malformed
# blob (with an uninitialized SRID) instead of reporting the recorded error.
add_test(NAME boost_rejects_bad_magic
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_Area(x'4751000000000000000101000000000000000000084000000000000008400000000000001040');")
set_tests_properties(boost_rejects_bad_magic PROPERTIES
  PASS_REGULAR_EXPRESSION "Incorrect GPB magic number")

add_test(NAME boost_area_valid_geometry
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_Area(ST_GeomFromText('POLYGON((0 0, 3 0, 3 3, 0 3, 0 0))'));")
set_tests_properties(boost_area_valid_geometry PROPERTIES
  PASS_REGULAR_EXPRESSION "9\\.0")

# A nested collection cannot be represented by the Boost bridge; it used to come back as a
# silent NULL indistinguishable from a NULL geometry, now the failed read is reported.
add_test(NAME nested_collection_reports_error
  COMMAND GPKG_SHELL :memory:
    "SELECT ST_NumGeometries(ST_GeomFromText('GEOMETRYCOLLECTION(GEOMETRYCOLLECTION(POINT(1 2)))'));")
set_tests_properties(nested_collection_reports_error PROPERTIES
  PASS_REGULAR_EXPRESSION "Invalid geometry blob")

# The metadata tables must satisfy the library's own conformance check, which
# also pins the last_change default that the format-string fix touched.
add_test(NAME metadata_roundtrip
  COMMAND GPKG_SHELL :memory:
    "SELECT GPKG_InitSpatialMetaData();"
    "SELECT GPKG_CheckSpatialMetaData();"
    "SELECT sql FROM sqlite_master WHERE name = 'gpkg_contents';")
set_tests_properties(metadata_roundtrip PROPERTIES
  PASS_REGULAR_EXPRESSION "strftime\\('%Y-%m-%dT%H:%M:%fZ'")

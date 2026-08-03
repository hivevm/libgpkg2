-- Regression: registering a spatial index must not rely on double-quoted string literals.
--
-- The gpkg_extensions insert wrote its values double-quoted, which SQLite only accepts
-- through its double-quoted-string fallback. On a connection that turns that fallback off
-- the statement failed with "no such column", the savepoint rolled the whole index back,
-- and GPKG_CreateSpatialIndex reported that it could not register the rtree usage.
.bail off
.dbconfig dqs_dml off
SELECT GPKG_InitSpatialMetaData();
CREATE TABLE f(id INTEGER PRIMARY KEY);
SELECT GPKG_AddGeometryColumn('f', 'geom', 'POINT', 4326);
SELECT GPKG_CreateSpatialIndex('f', 'geom', 'id');
SELECT 'extensions: ' || count(*) FROM gpkg_extensions;

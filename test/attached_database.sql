-- Regression: metadata lookups must address the database the caller named.
--
-- The check for the referenced SRS queried gpkg_spatial_ref_sys unqualified while every
-- other statement was qualified, so it always consulted main. An SRS defined only in an
-- attached GeoPackage was reported as missing, and one defined only in main would have been
-- accepted and written into the attached database as a dangling srs_id.
.bail off
ATTACH 'attached_database.gpkg' AS aux;
SELECT GPKG_InitSpatialMetaData();
SELECT GPKG_InitSpatialMetaData('aux');
INSERT INTO aux.gpkg_spatial_ref_sys (srs_name, srs_id, organization, organization_coordsys_id, definition)
  VALUES ('ETRS89 / UTM zone 33N', 32633, 'EPSG', 32633, 'PROJCS[]');
CREATE TABLE aux.f(id INTEGER PRIMARY KEY);
SELECT GPKG_AddGeometryColumn('aux', 'f', 'geom', 'POINT', 32633);
SELECT 'registered: ' || count(*) FROM aux.gpkg_geometry_columns WHERE srs_id = 32633;

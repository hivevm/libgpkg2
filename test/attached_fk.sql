-- foreign_key_check must run against the database the caller named. Unqualified it checked
-- main only, so a dangling reference in the attached GeoPackage passed validation.
ATTACH 'attached_fk.gpkg' AS aux1;
SELECT GPKG_InitSpatialMetaData('aux1');
INSERT INTO aux1.gpkg_contents (table_name, data_type, identifier, srs_id) VALUES ('ghost', 'features', 'ghost', 99999);
SELECT GPKG_CheckSpatialMetaData('aux1', 1);

-- The GeoPackage header pragmas must address the database the caller named. Unqualified they
-- always applied to main, so the attached file never received the application id and main's
-- header was overwritten.
ATTACH 'attached_header.gpkg' AS aux1;
SELECT GPKG_InitSpatialMetaData('aux1');
PRAGMA aux1.application_id;
PRAGMA main.application_id;

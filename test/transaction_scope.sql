-- Regression: a failed transactional function must not leave a savepoint behind.
--
-- The rollback path used to omit the matching RELEASE, so the connection stayed
-- inside an open transaction. Writes issued afterwards were reported as
-- successful but discarded when the connection closed.
--
-- Reopening the same file drops the connection: with the defect the INSERT below
-- is rolled back and the count is 0, otherwise it is 1.
.bail off
.open transaction_scope.gpkg
DROP TABLE IF EXISTS keep;
CREATE TABLE keep(x);
SELECT GPKG_AddGeometryColumn('nosuchtable', 'geom', 'POINT', 4326);
INSERT INTO keep VALUES(1);
.open transaction_scope.gpkg
SELECT 'rows: ' || count(*) FROM keep;

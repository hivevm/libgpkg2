# Changelog

All notable changes to this project are documented in this file. The format follows
[Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/), versioning follows
[SemVer 2.0.0](https://semver.org/spec/v2.0.0.html) — see [ADR-0005](docs/adr/0005-versioning-and-releases.md).
Update the **Unreleased** section in the same change as any user-visible modification.

Changes made before this changelog was introduced are recorded in the historical
[`ChangeLog`](ChangeLog) file inherited from the original libgpkg project.

## [Unreleased]

### Fixed

- `ST_SRID(geometry, srid)` used to rewrite the blob header in place inside the buffer returned
  by `sqlite3_value_blob`, corrupting memory owned by SQLite or by the caller (visible as other
  references to the same bound value observing the new SRID). The result is now assembled in a
  fresh buffer and the argument is left untouched.
- `AddGeometryColumn` inserted into `gpkg_geometry_columns` before the `gpkg_contents` row it
  references existed, so it always failed under `PRAGMA foreign_keys = ON`; it also inserted the
  contents row unconditionally, so it failed on the primary key for a table the application had
  already registered in `gpkg_contents`. The contents row is now created first and only when
  missing.
- The WKT parser silently discarded failures inside CompoundCurve members, accepting invalid
  curve members whenever the underlying error did not reach the error stream and reporting a
  bogus trailing-text error otherwise.
- `InitSpatialMetaData` stamped `application_id`/`user_version` onto `main` instead of the named
  attached database, and `CheckSpatialMetaData` ran `foreign_key_check`/`integrity_check`
  against the wrong database for the same reason; all three pragmas are now schema-qualified.
- `ST_AsText` formatted coordinates with ten significant digits, silently altering any
  coordinate that needs more; coordinates are now written with the fewest digits that still
  parse back to the same value, so the text round-trips exactly.
- A geometry ordinate consisting entirely of NaN values handed the envelope initialization
  sentinels (±DBL_MAX) to callers: `ST_MinX` returned 1.8e308 and broke the R-tree triggers,
  and a NaN-only Z rejected the whole geometry with a message-less error. Such ordinates are
  now reported as absent.
- The Boost bridge ignored blob-header read failures and computed results over malformed blobs
  with an uninitialized SRID; malformed blobs are now reported as errors, as is a nested
  geometry collection the bridge cannot represent (previously a silent NULL).
- The Boost geometry functions leaked their intermediate and result allocations whenever a
  Boost.Geometry algorithm threw (including the `ST_Union` aggregate's whole accumulation
  chain), and the writer reset performed an undefined-behaviour `memset` over `boost::variant`
  objects that covered only half the stack.

### Added

- Agentic-coding governance adopted from the NUC template: agent rules (`AGENTS.md`),
  specification and ADRs (`docs/`), conformance tracking against the OGC GeoPackage standard
  (`docs/CONFORMANCE.md`), repository consistency checks (`scripts/`, enforced by the Checks
  workflow), git conventions, supply-chain pinning with Dependabot, and this changelog.

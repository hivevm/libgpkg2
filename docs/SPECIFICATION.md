# Specification — libgpkg

> The **specification** of this project: the problem it solves, where it is going, and the
> vocabulary everyone (humans and agents) must use. This document is the **constitution** —
> every Architecture Decision Record in [`adr/`](adr/) derives from it and must not contradict it.

## Problem

Geodata is often exchanged in proprietary or fragmented formats that require heavyweight GIS
stacks to read and write. The [OGC GeoPackage standard](https://www.ogc.org/standards/geopackage)
solves this with an open, self-describing, SQLite-based container format — but plain SQLite knows
nothing about GeoPackage tables, geometry encoding, or spatial functions. Applications that embed
SQLite need a small, portable extension that makes a GeoPackage readable and writable without
pulling in a full GIS suite. The original Luciad `libgpkg` filled this gap but has been
unmaintained for years and predates current SQLite, Boost, and GeoPackage revisions.

## Mission

Provide a minimal, portable SQLite 3 extension that lets any SQLite-embedding application read,
write, and query OGC GeoPackage files.

## Vision

libgpkg is the drop-in GeoPackage capability for SQLite on every relevant platform — Linux,
macOS, Windows, Android, and iOS. It stays small enough to embed in mobile applications while
conforming to the GeoPackage Encoding Standard tracked in [`CONFORMANCE.md`](CONFORMANCE.md),
with geometry operations delegated to Boost.Geometry instead of hand-rolled algorithms.

## Strategy

- Revive the proven Luciad `libgpkg` code base and modernize it incrementally (current SQLite,
  CMake build, C++17/Boost.Geometry for geometry operations, optional ICU for collations).
- Keep the extension minimal: implement what the GeoPackage standard requires and what embedding
  applications concretely need — nothing speculative.
- Track conformance to the pinned GeoPackage revision explicitly in
  [`CONFORMANCE.md`](CONFORMANCE.md); deviations are recorded decisions, not accidents.
- Build for all target platforms from one CMake tree, verified by CI.

## Core Concepts & Vocabulary

- **GeoPackage** — an SQLite database file laid out according to the OGC GeoPackage Encoding
  Standard (`gpkg_contents`, `gpkg_spatial_ref_sys`, and related tables).
- **Extension (SQLite extension)** — the loadable library `libgpkg` that registers GeoPackage
  functions and behavior in a SQLite connection; not to be confused with a *GeoPackage extension*
  (an entry in `gpkg_extensions`).
- **GeoPackage Binary (GPB)** — the standard's BLOB encoding of a geometry: GP header (SRID,
  envelope, flags) followed by the well-known binary of the geometry.
- **WKB / WKT** — well-known binary / well-known text encodings of geometries, used by the
  `ST_GeomFromWKB` / `ST_GeomFromText` / `ST_AsBinary` / `ST_AsText` functions.
- **Spatial function** — an SQL function prefixed `ST_` registered by the extension and operating
  on geometry BLOBs.
- **SRS / SRID** — a spatial reference system and its integer identifier, kept in
  `gpkg_spatial_ref_sys`.
- **Spatial index** — the RTree-based index on a geometry column, per the standard's
  `gpkg_rtree_index` extension.
- **Shell** — the bundled `gpkg` command-line program: a SQLite shell statically linked with the
  extension.

## Goals / Success Criteria

1. A SQLite 3 application can load the extension (`load_extension` / `sqlite3_load_extension`)
   and create, read, and update GeoPackage files whose required tables and geometry BLOBs conform
   to the pinned revision in [`CONFORMANCE.md`](CONFORMANCE.md).
2. The `ST_*` function set documented in [`README.md`](../README.md) works on GeoPackage
   geometries, backed by Boost.Geometry.
3. The library builds from one CMake tree for Linux, macOS, Windows, Android, and iOS
   (see the CI workflows in `.github/workflows/`).
4. The bundled `gpkg` shell opens a GeoPackage and answers spatial queries out of the box.
5. Tests pass and the conformance status is kept current in [`CONFORMANCE.md`](CONFORMANCE.md).

## Non-Goals

- Not a full GIS: no rendering, no reprojection pipelines, no format conversion beyond
  GeoPackage/WKB/WKT.
- No support for non-SQLite containers or other geodata formats (Shapefile, GeoJSON, …).
- No network services — this is an embeddable library, not a server.
- No compatibility promises for behavior of the original libgpkg beyond what
  [`CONFORMANCE.md`](CONFORMANCE.md) lists as in scope.

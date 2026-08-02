# Conformance — libgpkg against the OGC GeoPackage Encoding Standard

> This document anchors the implementation to the **external source** it derives from: another
> project, a published specification, or a reference implementation. It answers three questions —
> *which* revision of that source this project is measured against, *how far* it has come, and
> *where it deliberately differs*. [`SPECIFICATION.md`](SPECIFICATION.md) remains the constitution:
> it says what this project is *for*; this file says what it *owes an outside authority*. Where the
> two collide, the specification wins and the deviation is recorded below.

## Source

- **Name:** OGC GeoPackage Encoding Standard (OGC 12-128)
- **Origin:** <https://www.ogc.org/standards/geopackage> (HTML rendering at
  <https://www.geopackage.org>)
- **Pinned revision:** GeoPackage Encoding Standard version 1.4.0
- **Last reconciled:** 2026-08-02 by Claude (initial pin at template adoption — the inventory
  below is an estimate from the code base, not yet verified unit by unit; see Open questions)
- **In scope:** the Core requirement classes (container, `gpkg_spatial_ref_sys`,
  `gpkg_contents`), the Features option (vector features with GeoPackage Binary geometry
  encoding), and the registered extensions the code base carries tables for (RTree spatial
  index, Schema, Metadata).
- **Out of scope:** Tiled Gridded Coverage Data and WKT for Coordinate Reference Systems beyond
  what `gpkg_spatial_ref_sys` minimally requires; rendering or reprojection of tile content.

**Reconciling** means: fetch the current upstream revision, diff it against the pinned one, add or
update the rows below, then move the pin. Moving the pin is therefore a reviewable change, and its
diff is exactly the answer to "what did the outside world do while we were not looking?". Reconcile
on a rhythm the source's pace warrants, and always before planning a milestone.

## Conformance inventory

One row per unit of the external source — a specification section, a feature, a test-suite entry, a
public API surface — at whatever granularity the source itself is organized in. Use the source's own
identifiers verbatim, so a diff against a newer upstream revision maps onto this table mechanically
instead of by interpretation.

**Status legend:** 🟢 conformant · 🟡 partial · 🔵 planned · ⚪ out of scope · 🔴 deviates

| Unit | Source ref | Status | Notes |
|------|------------|--------|-------|
| Core: SQLite container | Clause 1.1.1 | 🟡 partial | Container, application id, and user version handling exist; conformance against 1.4.0 not yet verified |
| Core: `gpkg_spatial_ref_sys` | Clause 1.1.2 | 🟡 partial | Table created and used by the extension; default SRS entries to be verified |
| Core: `gpkg_contents` | Clause 1.1.3 | 🟡 partial | Table created and maintained; column-level conformance to be verified |
| Features: `gpkg_geometry_columns` | Clause 2.1.5 | 🟡 partial | Table supported; geometry type constraints to be verified |
| Features: GeoPackage Binary geometry encoding | Clause 2.1.3 | 🟡 partial | GPB read/write with WKB/WKT conversion implemented (`gpkg_geom.c`, `wkb.c`, `wkt.c`) |
| Tiles: `gpkg_tile_matrix_set`, `gpkg_tile_matrix` | Clause 2.2 | 🟡 partial | Tables are created; tile pyramid behavior beyond table creation is untracked |
| Extension mechanism: `gpkg_extensions` | Clause 2.3 | 🟡 partial | Table supported |
| Extension: RTree spatial indexes | Annex F.3 (`gpkg_rtree_index`) | 🟡 partial | Spatial index creation implemented on SQLite RTree |
| Extension: Schema | Annex F.9 (`gpkg_data_columns`, `gpkg_data_column_constraints`) | 🟡 partial | Tables supported |
| Extension: Metadata | Annex F.8 (`gpkg_metadata`, `gpkg_metadata_reference`) | 🟡 partial | Tables supported |
| Attributes (non-spatial tables) | Clause 2.4 | 🔵 planned | Not yet inventoried |
| Tiled Gridded Coverage Data | OGC 17-066 | ⚪ out of scope | Not a goal of this extension |

Completeness of this table is itself a claim: a unit that exists upstream and is missing here is an
unreconciled gap, not an implicit "out of scope".

## Deviations

Deliberate departures from the source. Each entry states what the source requires, what this project
does instead, and why. A deviation that is architecture-relevant additionally needs an ADR (see
[`adr/`](adr/)) and cites it here — this list records *that* we deviate, the ADR records *why it was
allowed*. **A difference not listed here is a defect, not a decision.**

- None recorded yet — differences found during the first full reconciliation must be entered here
  or fixed.

## Open questions

Where the source is ambiguous, contradictory, or silent, and what this project assumed in the
meantime. Resolve upstream where possible — an answer from the source removes an assumption from
this list.

- The inventory above was estimated from the tables and functions present in the code
  (`gpkg/spatialdb.c` and related files) when the conformance tracking was adopted; a first full
  reconciliation against the pinned 1.4.0 text — verifying each row and its exact clause
  numbers — is still outstanding.

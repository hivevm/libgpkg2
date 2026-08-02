# ADR-0006: Conformance to an external source is tracked in `docs/CONFORMANCE.md`

- **Status:** 🟢 accepted
- **Date:** 2026-08-01
- **Deciders:** Maintainer

## Context

Some projects are not free-standing: they implement a published specification, port another
project, or reimplement a reference implementation. For those, a second authority exists next to
[`docs/SPECIFICATION.md`](../SPECIFICATION.md) — an external source that keeps moving on its own
schedule and against which "are we done?" is a measurable question rather than a matter of taste.

Without a written anchor, three failures follow predictably. The upstream revision the
implementation was derived from is remembered by nobody, so nothing can be diffed when the source
advances. Progress is claimed per feeling instead of per unit, so gaps surface late. And deliberate
departures become indistinguishable from bugs: an agent "fixing" an intentional deviation is a
recurring, expensive failure mode. Agents are especially exposed here, since they cannot recover
any of this from the code — a missing feature and a rejected feature look identical in a source
tree.

The anchor must be a plain tracked document (the template adds no toolchain), must stay honest
about scope, and must not exist at all in projects that derive from nothing.

## Decision

We will record conformance in **`docs/CONFORMANCE.md`**, a sibling of the specification, shipped as
the optional policy module `conformance` and adopted at the template's first-run module selection —
projects that derive from nothing never see the file. It carries three parts:

- **Source** — the external project or specification, its origin, and a **pinned revision**
  (version, tag, commit SHA, or document revision) plus the date of the last reconciliation, and an
  explicit in-scope/out-of-scope split. Reconciling means diffing the current upstream revision
  against the pin, updating the rows, and only then moving the pin — so the pin's diff is precisely
  the record of external progress since the last look.
- **Conformance inventory** — one row per unit of the source (specification section, feature,
  test-suite entry, API surface), at the granularity the source itself uses and under the source's
  own identifiers, with a status per row. Using upstream identifiers verbatim is what makes an
  upstream diff map onto the table mechanically. A unit that exists upstream and is absent here
  counts as an unreconciled gap, never as an implicit "out of scope".
- **Deviations** — deliberate departures, each stating what the source requires, what this project
  does instead, and why; architecture-relevant ones cite the ADR that allowed them. **A difference
  not listed is a defect, not a decision** — this sentence is the point of the whole section.

The authority chain is unchanged: the specification stays the constitution, accepted ADRs derive
from it, and `docs/CONFORMANCE.md` is *descriptive* — it states what the implementation owes an
outside authority and where it knowingly does not pay. A conflict between the source and the
specification is raised, not resolved silently; if the deviation stands, it is decided in an ADR
and listed here.

Keeping it current is part of the Definition of Done ([`AGENTS.md`](../../AGENTS.md) §5): a change
that implements or alters a derived unit updates its row in the same change.

## Alternatives considered

- **Track conformance in issues or a project board** — invisible to a coding agent reading the
  repository, and lost when the tracker changes; the anchor must live in the tree.
- **Fold it into `docs/SPECIFICATION.md`** — mixes the constitution ("what we want") with a moving
  external report ("what they have"), and forces free-standing projects to carry the concept.
- **Derive conformance from a test suite alone** — a green suite proves what is covered, not what
  the source contains; upstream units with no test map to nothing, and deviations carry no rationale.
- **One ADR per derived unit** — accurate but unusably granular; ADRs record decisions, not the
  moving state of an external work in progress.
- **No anchor, rely on commit history** — the failure this decision removes.

## Sources / Prior art

- **Protocol Implementation Conformance Statement (PICS)** — the ISO/IEC 9646 conformance-testing
  tradition of a structured statement asserting which requirements an implementation meets;
  the inventory table is a lightweight PICS proforma.
  <https://en.wikipedia.org/wiki/Protocol_implementation_conformance_statement>
- **Requirements Traceability Matrix**, as recommended by ISO/IEC/IEEE 29148 (requirements
  engineering) — one row per requirement, traced to the artifacts that satisfy it; the standard's
  own conformance clause distinguishes full from tailored conformance, which is the in-scope /
  out-of-scope split adopted here. <https://www.reqview.com/blog/requirements-traceability-matrix/>
- **Web Platform Tests / wpt.fyi implementation reports** — per-unit conformance of several
  implementations against one spec, tracked continuously rather than at release time.
  <https://testdev.tools/wpt-fyi/>
- **Technology Compatibility Kit (TCK)** — conformance asserted against an upstream-owned suite;
  the model this decision deliberately does *not* require, since most external sources ship none.
  <https://en.wikipedia.org/wiki/Technology_Compatibility_Kit>
- **Matrix `complement`** — a compliance suite maintained alongside a moving specification, an
  example of upstream identifiers as the unit of tracking.
  <https://github.com/matrix-org/complement>

## Consequences

- Positive: the derived-from revision is written down, so external progress can be diffed instead of
  guessed; progress is countable per unit; intentional deviations are protected from well-meaning
  "fixes" by humans and agents alike; a new contributor or agent learns the project's relationship
  to its source from one file.
- Negative / trade-offs: the inventory is manual and goes stale between reconciliations — it is a
  *reviewed claim*, not a measurement, and a stale pin overstates conformance; granularity has to be
  chosen per source and is hard to change later; for large sources the table grows long.
- Follow-ups: projects whose source ships a machine-checkable suite should consider deriving the
  status column from a CI run instead of maintaining it by hand — a decision for those projects, in
  their own ADR, not for this template.

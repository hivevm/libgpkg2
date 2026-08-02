# Architecture Decision Records

This directory contains all Architecture Decision Records (ADRs) for this project.
Accepted ADRs are **binding** for humans and coding agents alike (see [`AGENTS.md`](../../AGENTS.md)
in the repository root). ADRs derive from the specification in [`docs/SPECIFICATION.md`](../SPECIFICATION.md).

## Process

1. Copy [`template.md`](template.md) to `NNNN-short-title.md` (next free number).
2. Fill in context, decision, alternatives, and consequences. Set status `proposed`.
3. A human reviewer accepts or rejects the ADR. **Only humans change the status.**
4. Add the ADR to the index below, with its status shown via the colored bullet from the legend.
5. A decision is changed by a *new* ADR that supersedes the old one — never by editing an
   accepted ADR.
6. **Once this template is in use, ADRs are immutable and their numbers are permanent.** Never
   renumber, delete, or merge ADRs — other ADRs, commits (`Implements ADR-NNNN`), and code may
   reference a number. Superseded ADRs stay as historical record (status `superseded by ADR-NNNN`);
   filter active ones via the Status column. To curb sprawl, supersede — do not consolidate.
   The numbers therefore run `0001..N` without gaps, which
   [`scripts/check-docs.sh`](../../scripts/check-docs.sh) verifies.
7. **Never reference an ADR number that does not exist yet.** Every `ADR-NNNN` reference must point
   to a file that is already present in this directory. Anticipated follow-up decisions are
   described by topic (e.g., "a follow-up ADR on session storage") in the Consequences section —
   the concrete number is cited only once that ADR file exists.

## Index

**Status legend:** 🟢 accepted · 🟡 proposed · 🔴 rejected · ⚪ superseded

| ADR | Title | Status |
|-----|-------|--------|
| [0001](0001-agent-governance-model.md) | Specification + ADRs governed through a single `AGENTS.md` | 🟢 accepted |
| [0002](0002-dev-container-runtime.md) | Debian Dev Container without host Docker access | 🟢 accepted |
| [0003](0003-git-conventions.md) | Git conventions: branches, Conventional Commits, squash merge | 🟢 accepted |
| [0004](0004-secrets-and-supply-chain.md) | Secrets handling and supply-chain pinning | 🟢 accepted |
| [0005](0005-versioning-and-releases.md) | Versioning and release process | 🟢 accepted |
| [0006](0006-external-conformance-tracking.md) | Conformance to an external source is tracked in `docs/CONFORMANCE.md` | 🟢 accepted |
| [0007](0007-cpp-dev-container-and-toolchain.md) | C++ Dev Container image and CMake toolchain | 🟡 proposed |


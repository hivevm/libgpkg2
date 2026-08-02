# ADR-0005: Versioning and release process

- **Status:** 🟢 accepted
- **Date:** 2026-07-21
- **Deciders:** Maintainer

## Context

The template has no versioning or release policy: [`SECURITY.md`](../../SECURITY.md) deferred its
"Supported versions" section to the first release, there is no changelog, and no tag convention.
Agents need a mechanical rule for what to update when they change something user-visible, and the
security policy needs a support statement it can point to. The template's commit-convention seed
decision (the `git-conventions` module, where adopted) fixes Conventional Commit types, which map
naturally onto semantic versions; the pinning and secrets rules of the supply-chain module are
unaffected. There are no releases yet, so the process must stay minimal (YAGNI).

## Decision

We will version releases with [SemVer 2.0.0](https://semver.org/spec/v2.0.0.html) (pre-1.0: minor
versions may break, per SemVer item 4), maintain a
[Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/)-formatted
[`CHANGELOG.md`](../../CHANGELOG.md) whose `[Unreleased]` section is updated in the same change as
any user-visible modification, tag releases with annotated `vX.Y.Z` tags, and let a human perform
releases manually (promote `[Unreleased]`, tag, publish — process in
[`CONTRIBUTING.md`](../../CONTRIBUTING.md)). Only the latest release receives security fixes,
shipped as a new release.

## Alternatives considered

- **CalVer** — carries no compatibility signal; Conventional Commit types (where the
  git-conventions module is adopted: feat → minor, fix → patch, `!` → major) map onto SemVer,
  not onto dates.
- **Automated release tooling (semantic-release, release-please)** — a new dependency and CI
  surface for a project with zero releases. The Conventional Commit history keeps this door open;
  adopting such a tool later is its own ADR.
- **A dedicated `docs/RELEASING.md`** — a whole file for a ~10-line manual process is premature;
  a section in `CONTRIBUTING.md` suffices until the process grows.
- **Supporting multiple release lines** — maintenance branches and backports for a template
  nobody has released yet; latest-release-only is honest and cheap.

## Sources / Prior art

- [Semantic Versioning 2.0.0](https://semver.org/spec/v2.0.0.html) — including item 4 on 0.x
  development versions.
- [Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/) — changelog format and the
  `[Unreleased]` convention.
- [Conventional Commits v1.0.0](https://www.conventionalcommits.org/en/v1.0.0/#how-does-this-relate-to-semver)
  — the commit-type → SemVer mapping.

## Consequences

- Positive: `SECURITY.md` gets a concrete support statement; user-visible changes leave a
  changelog trail as they land (no reconstruction at release time); version numbers carry
  compatibility meaning; releases stay a deliberate human act (consistent with
  [`AGENTS.md`](../../AGENTS.md) §6).
- Negative / trade-offs: every user-visible PR carries a small changelog chore; manual releases
  do not scale to high release frequency (revisit with tooling via a new ADR if that happens).
- Follow-ups: release automation (changelog generation from Conventional Commits) as its own ADR
  if release frequency grows.

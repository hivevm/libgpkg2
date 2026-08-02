# Contributing

Thanks for your interest in contributing. This project is **specification- and ADR-driven**, and
the working rules in [`AGENTS.md`](AGENTS.md) are the single source of truth — they govern coding
agents and human contributors alike. This file adds only what is contributor-specific and points
to the rest instead of restating it.

## Before you start

Read, in this order:

1. [`docs/SPECIFICATION.md`](docs/SPECIFICATION.md) — the constitution: problem, goals, and vocabulary.
2. [`docs/adr/`](docs/adr/) — the Architecture Decision Records. **Accepted ADRs are binding.**
3. [`AGENTS.md`](AGENTS.md) — the working rules: principles, the ADR process, the quality bar
   (Definition of Done), and all conventions, including language and vocabulary.
4. [`docs/CONFORMANCE.md`](docs/CONFORMANCE.md) — the external specification or project this
   implementation derives from, how far it conforms, and where it deliberately differs.

Authority runs **specification → accepted ADRs → individual change**.

## Workflow

1. **Open an issue first** for anything non-trivial, so scope and intent can be agreed before work
   starts.
2. **Follow [`AGENTS.md`](AGENTS.md)**: when a decision needs an ADR and how it is reviewed (§3),
   when a change counts as done (§5), and the git and project conventions (§6). The consistency
   checks in [`scripts/`](scripts/) run in CI on every pull request — run them locally before
   pushing; they need only bash and coreutils, already in the Dev Container (only the shell lint
   additionally needs [ShellCheck](https://www.shellcheck.net) and skips itself where that is
   missing — CI runs it always).
3. **Open a pull request** against `main`, fill in the PR template, and link the issue/ADR.
   PRs are **squash-merged**, so the PR title must itself be a valid Conventional Commit subject —
   it becomes the permanent commit on `main`.

## Releases

Releases are performed manually by a maintainer (never by an agent — see
[`AGENTS.md`](AGENTS.md) §6), versioned with [SemVer 2.0.0](https://semver.org/spec/v2.0.0.html)
per [ADR-0005](docs/adr/0005-versioning-and-releases.md). Before 1.0.0, minor versions may contain
breaking changes.

1. Move the content of `## [Unreleased]` in [`CHANGELOG.md`](CHANGELOG.md) into a new
   `## [X.Y.Z] - YYYY-MM-DD` section (leave `[Unreleased]` in place, empty).
2. Commit (`chore: release vX.Y.Z`) and merge via the normal PR workflow.
3. Tag the release commit: `git tag -a vX.Y.Z -m "vX.Y.Z"` and push the tag.
4. Create the GitHub release from the changelog section: `gh release create vX.Y.Z`.

## Reporting security issues

Do **not** open a public issue for vulnerabilities. Follow [`SECURITY.md`](SECURITY.md).

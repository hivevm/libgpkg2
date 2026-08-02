# ADR-0003: Git conventions: branches, Conventional Commits, squash merge

- **Status:** 🟢 accepted
- **Date:** 2026-07-21
- **Deciders:** Maintainer

## Context

Coding agents produce many small commits and pull requests. Without fixed, machine-checkable git
conventions the history degrades quickly: inconsistent subjects, unreadable branch names, and merge
topologies that make review and later changelog work harder. The repository's philosophy is to
enforce rules as CI rather than convention (see [`AGENTS.md`](../../AGENTS.md) and the existing
documentation checks), and any check must stay pure bash + coreutils — adding a linting toolchain
would itself be an architecture-relevant dependency.

## Decision

We will name work branches `type/short-topic` (kebab-case), write commit subjects as
[Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) (types
`feat fix docs refactor perf test build ci chore revert`, optional scope, optional `!`, at most
72 characters, no trailing period), merge pull requests into `main` exclusively by squash merge
with a PR title that is itself a valid Conventional Commit subject, and prohibit force pushes and
history rewrites on `main` and any other shared branch. The conventions are enforced on pull
requests by [`scripts/check-git-conventions.sh`](../../scripts/check-git-conventions.sh); the
force-push prohibition is enforced by a GitHub ruleset on `main`, configured in the repository
settings (repository files cannot enforce it — see the setup checklist in
[`README.md`](../../README.md)).

## Alternatives considered

- **Free-form commit subjects with prose rules ("imperative mood")** — not machine-checkable with
  a bash regex; a rule that CI cannot enforce decays under agent-generated volume.
- **50-character subject limit** — the classic git recommendation, but too much friction with a
  mandatory `type(scope):` prefix; 72 is the widely accepted hard wrapping limit.
- **Merge commits or rebase merge instead of squash** — merge commits make `main` non-linear and
  hide the one-change-per-PR unit; rebase merge replays intermediate WIP commits onto `main`.
  Squash yields one Conventional-Commit line per PR, matching the "small, reviewable changes" rule.
- **Commit-lint tooling (commitlint, husky, etc.)** — a new toolchain dependency for something a
  single bash regex covers; contradicts simplicity-first.

## Sources / Prior art

- [Conventional Commits v1.0.0](https://www.conventionalcommits.org/en/v1.0.0/) — subject format
  and type vocabulary.
- [GitHub: about squash merge](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/configuring-pull-request-merges/about-merge-methods-on-github)
  — squash uses the PR title as the resulting commit subject.
- [GitHub: rulesets](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-rulesets/about-rulesets)
  — blocking force pushes and requiring status checks is a settings-level control.

## Consequences

- Positive: history on `main` is linear, one Conventional-Commit line per PR; conventions are
  CI-enforced instead of reviewed by hand; commit types later feed versioning and changelog work
  without new tooling.
- Negative / trade-offs: a strict regex occasionally rejects harmless subjects; contributors and
  agents must learn the type vocabulary; branch protection must be configured manually in the
  GitHub settings.
- Follow-ups: a supply-chain ADR (pinning and secrets rules for CI) and a release-policy ADR
  (versioning and changelog, building on the commit types chosen here).

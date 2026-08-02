# ADR-0004: Secrets handling and supply-chain pinning

- **Status:** 🟢 accepted
- **Date:** 2026-07-21
- **Deciders:** Maintainer

## Context

Coding agents read, write, and execute far more than a careful human would, so a secret that lands
in a tracked file, a commit message, or CI output spreads fast — and mutable dependency references
(a moved tag, an unpinned action) are a supply-chain attack vector that CI executes with repository
permissions. The template so far only gitignores `.env*` and forbids storing `gh` tokens
([`AGENTS.md`](../../AGENTS.md) §6); there is no lockfile rule, no pinning rule, and no scanning.
Checks must remain pure bash + coreutils — a scanning toolchain would itself be a new dependency.

## Decision

We will forbid writing secrets into tracked files, commit messages, ADRs, logs, or CI output
(secrets live in environment variables, gitignored `.env*` files, or GitHub Actions secrets, and a
leaked secret is rotated immediately — deleting it from the tip is not remediation); rely on
GitHub's built-in secret scanning with push protection, enabled in the repository settings; require
that when the chosen toolchain has a lockfile it is committed and CI installs from it; pin every
GitHub Action to a full commit SHA with a trailing version comment, enforced by
[`scripts/check-workflow-pins.sh`](../../scripts/check-workflow-pins.sh); and let Dependabot
([`.github/dependabot.yml`](../../.github/dependabot.yml)) keep those pins current.

## Alternatives considered

- **gitleaks (or similar) as a CI scanning step** — real scanning coverage, but a new pinned
  dependency with rule maintenance and false-positive tuning, for a template that ships no code
  yet. A project can add it later via its own ADR; GitHub's built-in scanning covers the baseline.
- **A bash-regex secret scan in `scripts/`** — dependency-free but weak: a handful of token
  patterns produce false confidence while missing most real leaks. Worse than relying on GitHub's
  maintained pattern set.
- **Pinning actions to major tags (`@v7`)** — the previous state. Tags are mutable: whoever
  controls (or compromises) the action repository can repoint them, as in the `tj-actions/changed-files`
  compromise (March 2025). A commit SHA is immutable; the version comment keeps it readable.
- **No update automation for pins** — SHA pins rot silently and manual updates stop happening.
  Dependabot is GitHub-native (no new toolchain) and updates SHA + comment together.

## Sources / Prior art

- [GitHub: secret scanning and push protection](https://docs.github.com/en/code-security/secret-scanning/introduction/about-secret-scanning)
- [GitHub: security hardening for GitHub Actions](https://docs.github.com/en/actions/security-for-github-actions/security-guides/security-hardening-for-github-actions)
  — recommends pinning third-party actions to a full commit SHA.
- CISA/GitHub advisories on the `tj-actions/changed-files` compromise (CVE-2025-30066) — a
  retagged action exfiltrating CI secrets, the concrete case for SHA pinning.
- [GitHub: Dependabot version updates for GitHub Actions](https://docs.github.com/en/code-security/dependabot/working-with-dependabot/keeping-your-actions-up-to-date-with-dependabot)

## Consequences

- Positive: immutable action references, CI-enforced; a baseline secrets policy that agents can
  follow mechanically; lockfile rule is language-agnostic and ready for whatever toolchain a
  project adds; pin updates arrive as reviewable Dependabot PRs.
- Negative / trade-offs: SHA-pinned `uses:` lines are less readable (mitigated by the mandatory
  version comment); secret scanning and push protection are repository settings that must be
  enabled manually (README setup checklist); Dependabot PRs add review traffic.
- Follow-ups: a release-policy ADR (versioning, changelog, and support statement); per-project
  scanning tooling (e.g. gitleaks) can be proposed in its own ADR when a project's risk warrants it.

# Security Policy

## Reporting a vulnerability

Please report security vulnerabilities **privately** — do not open a public issue or pull request.

- Use GitHub's [private vulnerability reporting](https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing-information-about-vulnerabilities/privately-reporting-a-security-vulnerability)
  (**Security → Report a vulnerability**) if enabled, or
- email the maintainer: <!-- TODO: add a security contact address -->.

Please include enough detail to reproduce the issue (affected version/commit, steps, and impact).
We aim to acknowledge reports within a reasonable time frame and will coordinate a fix and
disclosure with you.

## Dev Container & agent execution

This template runs coding agents inside the Dev Container defined in
[`.devcontainer/devcontainer.json`](.devcontainer/devcontainer.json). Two properties shape its
security posture:

- **The Dev Container has no access to the host container engine.** The host Docker/Podman socket is
  **not** mounted into the container (see [ADR-0002](docs/adr/0002-dev-container-runtime.md)),
  so code or agents running inside cannot control the host engine. Host containers are managed from a
  host-side VS Code extension (see the README), keeping that capability outside the container's reach.
  The container is still not a strong security boundary, so run only agents and code you trust in it.
- **Git and GitHub writes require explicit human approval.** The agent must get a go-ahead before
  each commit, push, or `gh` action ([`AGENTS.md`](AGENTS.md) §6). For Claude Code this is enforced,
  not just documented: `.claude/settings.json` prompts on `git add`,
  `git commit`, `git push`, and `gh`. Authentication uses `gh`'s web flow with no stored tokens.

## Secrets & supply chain

The binding rules — no secrets in tracked files, commits, or CI output; rotate anything that
leaks; lockfiles committed; every GitHub Action pinned to a full commit SHA (CI-enforced by
[`scripts/check-workflow-pins.sh`](scripts/check-workflow-pins.sh)) — live in
[`AGENTS.md`](AGENTS.md) §7 (decided in [ADR-0004](docs/adr/0004-secrets-and-supply-chain.md)).
GitHub's **secret scanning with push protection** complements them; it is a repository setting,
listed in the setup checklist in [`README.md`](README.md).

## Supported versions

Security fixes are delivered on the current development line only — older states are not patched.
Only the **latest release** (the highest `vX.Y.Z` tag — see [`CHANGELOG.md`](CHANGELOG.md) and
[ADR-0005](docs/adr/0005-versioning-and-releases.md)) receives security fixes, delivered as a new
release; before 1.0.0 there are no support guarantees beyond that.

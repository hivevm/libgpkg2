# Agent Guide

> **Single source of truth for all coding agents.** Most agents read this file natively; the
> concrete agent set and how each is wired to read these rules is recorded in
> [ADR-0001](docs/adr/0001-agent-governance-model.md) — not here. Where an agent needs a pointer
> (e.g. Claude Code's `.claude/CLAUDE.md`), never duplicate rules into it — add them here.

Human-facing setup — prerequisites, the Dev Container, and the project description —
lives in [`README.md`](README.md). Do not duplicate that information here.

## 1. Principles

These three principles govern every rule below. When they tension with speed or cleverness, they win.

- **Simplicity first.** Build the simplest thing that satisfies the specification. Prefer fewer moving
  parts; add a dependency, an abstraction, or a layer of indirection only for a concrete, present need,
  never a speculative future one (YAGNI). Removing code is progress. When complexity is genuinely
  warranted it is architecture-relevant — justify it in an ADR (§3).
- **Reflection.** Think before and after acting. Before: is this the simplest path that actually serves
  the specification's goals? After: did it, and at what cost? Weigh alternatives and consequences
  instead of committing to the first solution, and make your reasoning explicit rather than silent.
- **Critical stance.** Take nothing at face value — not the human's framing, not your own prior output,
  not the existing code. Verify claims against the specification, the ADRs, the code, and authoritative
  sources; surface conflicts, risks, and uncertainty instead of smoothing them over. Disagree when the
  evidence warrants and say why, and flag what you could not confirm.

## 2. Start here

Before any non-trivial work, read:

1. [`docs/SPECIFICATION.md`](docs/SPECIFICATION.md) — the **specification** (problem, goals, core
   concepts, vocabulary, success criteria).
2. [`docs/adr/`](docs/adr/) — the Architecture Decision Records. **Accepted ADRs are binding**.
3. [`docs/CONFORMANCE.md`](docs/CONFORMANCE.md) — what this project owes the **external source** it
   derives from: the pinned upstream revision, the per-unit conformance inventory, and the
   deliberate deviations. A difference from the source that is listed there is a decision, not a
   defect — never "fix" one without an ADR that reverses it.

## 3. ADR rules

Authority runs **specification → accepted ADRs → task**: the specification in
[`docs/SPECIFICATION.md`](docs/SPECIFICATION.md) is the constitution, accepted ADRs derive from it,
and every task respects both. The full ADR mechanics — numbering, template, lifecycle, and the
index — live in [`docs/adr/README.md`](docs/adr/README.md).

1. **Create an ADR before any architecture-relevant decision** — adding a dependency or framework,
   designing or changing a public interface, choosing a persistence/synchronization strategy or a
   protocol/data format, or anything that constrains future technology choices (non-exhaustive). Copy
   [`docs/adr/template.md`](docs/adr/template.md), set status `proposed`, then **stop and ask for
   human review** before implementing.
2. **Never violate an accepted ADR.** If a task conflicts with one, do not silently work around it —
   propose a new ADR (status `proposed`) that supersedes it.
3. **You may set status `proposed` only; only a human reviewer changes it.** Never edit an accepted
   ADR — change a decision with a *new* ADR that supersedes it. **ADR numbers are permanent — never
   renumber, delete, or merge ADRs;** superseded ones remain as historical record.
4. **The specification wins.** If the specification and an ADR conflict, raise the conflict — do not
   choose silently.
5. **Use the project vocabulary** from the specification consistently in code, comments, and documentation.
6. **Calibrate — not everything needs an ADR.** ADRs are for decisions that are *costly to reverse* or
   *constrain future choices*. Routine, local, easily reversible work does **not** need one: implementing
   within an accepted ADR, bug fixes, refactorings that preserve public interfaces, tests, docs,
   formatting, or a dev-only tool that no shipped code depends on. Rule of thumb: if a change locks in
   nothing and could be undone in a single follow-up commit, skip the ADR. When genuinely unsure, prefer
   a short `proposed` ADR over a silent decision.
7. **Keep the ADR index current.** When you add, supersede, or change the status of an ADR, update the
   index in [`docs/adr/README.md`](docs/adr/README.md) **in the same change** — never as a separate
   afterthought.

## 4. Working style

- **Reply to the human in their own language.** The human may write in *any* language; always answer
  in that same language. Everything else stays English: write **all artifacts** (code, comments,
  commits, docs, ADRs, PRs) in **English**, without exception.
- Prefer **small, reviewable changes**; for larger ones, plan and explore the codebase before
  editing. Reference the relevant ADR(s) in commit messages and pull request descriptions
  (e.g., `Implements ADR-NNNN`).
- **Work interactively and iteratively**: proceed in small steps, surface your reasoning, and seek
  feedback early rather than delivering large changes at once.
- **Research before any conceptual design.** For every conceptual or design decision, first research
  the state of the art and established solutions — including external/web sources — instead of relying
  on assumptions or memory alone. Capture the relevant findings and cite them (in the ADR when one
  applies).
- When in doubt about scope or intent, ask the human before implementing.

## 5. Quality bar & Definition of Done

The exact build, test, and lint commands live in the **Build, Test & Run** section of
[`README.md`](README.md); this section defines *when* a change is done, not *how* to run the tools.

- **A change is done only when it builds, its tests pass, and linters/formatters are clean** —
  locally and in CI. Never hand off or propose merging red.
- **New behaviour ships with tests; bug fixes ship with a regression test** that fails before the fix
  and passes after. If a change is genuinely untestable, say so and explain why.
- **Never weaken the suite to make it pass.** Do not delete, skip, or loosen assertions to go green;
  fix the code or, if a test is genuinely wrong, correct it and explain the reasoning.
- **Treat a flaky test as a defect, not noise.** Do not paper over it with retries or by re-running
  until green — surface it and fix the root cause.
- **Keep the diff releasable.** No commented-out dead code, stray debug output, or `TODO` left as a
  substitute for a decision; unfinished work is tracked as an issue or a `proposed` ADR, not hidden in
  the tree.
- **A user-visible change updates the `[Unreleased]` section of [`CHANGELOG.md`](CHANGELOG.md)**
  in the same change ([ADR-0005](docs/adr/0005-versioning-and-releases.md)); internal-only changes
  (refactorings, tests, CI) do not.
- **A change to a derived unit updates its row in [`docs/CONFORMANCE.md`](docs/CONFORMANCE.md)**
  in the same change ([ADR-0006](docs/adr/0006-external-conformance-tracking.md)); a departure from
  the external source is either entered under *Deviations* with its reason or it is a bug. Never
  move the pinned upstream revision without reconciling the inventory against it.

## 6. Project rules

- **Git writes need explicit human approval, every time.** The agent may run git write operations
  (`add`, `commit`, `push`) and use the GitHub CLI (`gh`) to perform them, but must obtain the human's
  explicit go-ahead immediately before each commit or push. Never commit or push autonomously, and treat
  each approval as single-use — it does not carry over to the next commit or push.
- **GitHub interaction happens on explicit instruction only.** The Dev Container provides the GitHub
  CLI (`gh`); the agent may use it to interact with GitHub (e.g. pull requests, issues, releases) only
  when a human explicitly asks, and — as with git writes — each instruction is single-use and never
  implies the next.
- **Authenticate `gh` through its web flow.** Run `gh auth login` and choose *Login with a web
  browser*; a human then enters the displayed one-time code at <https://github.com/login/device> to
  authorize. Never request, store, or hard-code personal access tokens.
- **Branches and commits follow fixed conventions**
  ([ADR-0003](docs/adr/0003-git-conventions.md)). Work branches are named `type/short-topic`
  (kebab-case); commit subjects are [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/)
  — `type(scope)?: description` with types `feat fix docs refactor perf test build ci chore revert`,
  at most 72 characters, no trailing period. Enforced on pull requests by
  [`scripts/check-git-conventions.sh`](scripts/check-git-conventions.sh); run it locally before pushing.
- **Pull requests are squash-merged into `main`**, and the PR title must itself be a valid
  Conventional Commit subject — squash merge makes it the permanent history line on `main`.
- **Never force-push or rewrite history on `main` or any shared branch.** Force-pushing your own,
  not-yet-shared PR branch is fine. On `main` this is enforced by a GitHub ruleset configured in the
  repository settings (see the setup checklist in [`README.md`](README.md)) — repository files
  cannot enforce it.
- **Releases are human-only.** Versions follow SemVer with annotated `vX.Y.Z` tags
  ([ADR-0005](docs/adr/0005-versioning-and-releases.md)); the release steps live in
  [`CONTRIBUTING.md`](CONTRIBUTING.md). Agents never tag or publish a release.
- Build, test, and run commands live in the **Build, Test & Run** section of
  [`README.md`](README.md) — the single source for both humans and agents.

## 7. Secrets & supply chain

Decided in [ADR-0004](docs/adr/0004-secrets-and-supply-chain.md); see also [`SECURITY.md`](SECURITY.md).

- **Never write secrets into tracked files, commit messages, ADRs, logs, or CI output** — no
  tokens, API keys, passwords, or `.env` contents. Secrets live in environment variables,
  gitignored `.env*` files, or GitHub Actions secrets.
- **A leaked secret is compromised the moment it lands in git.** Rotate it first; deleting it from
  the tip of the branch is not remediation.
- **When the toolchain has a lockfile, it is committed** and CI installs from it (frozen install);
  manifest and lockfile change together in the same commit.
- **Pin every GitHub Action to a full commit SHA with a trailing version comment**
  (`uses: owner/action@<40-hex-sha> # vX.Y.Z`) — tags are mutable and therefore not a pin.
  Enforced in CI by [`scripts/check-workflow-pins.sh`](scripts/check-workflow-pins.sh);
  Dependabot keeps the pins current.
- **New dependencies and toolchains remain ADR-gated** (§3) — this section governs how
  dependencies are referenced, not whether they are added.

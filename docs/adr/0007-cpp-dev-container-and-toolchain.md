# ADR-0007: C++ Dev Container image and CMake toolchain

- **Status:** 🟡 proposed
- **Date:** 2026-08-02
- **Deciders:** Maintainer
- **Note:** Records the state already embodied in
  [`.devcontainer/devcontainer.json`](../../.devcontainer/devcontainer.json) and
  [`CMakeLists.txt`](../../CMakeLists.txt) at the time the agentic-coding template was adopted
  into this existing project.

## Context

[ADR-0002](0002-dev-container-runtime.md) was adopted from the template and decides two things:
the Dev Container has no access to the host Docker daemon, and it is based on the language-free
`mcr.microsoft.com/devcontainers/base:debian` image. libgpkg is a C/C++ project that predates the
template adoption: it already builds with CMake against SQLite, Boost, and optionally ICU, and its
Dev Container is based on the C++ dev container image so the toolchain works out of the box. The
base-image part of ADR-0002 therefore does not fit this project, while its security posture (no
host Docker socket) does.

## Decision

We will base the Dev Container on the `mcr.microsoft.com/devcontainers/cpp` image, pinned by tag
**and digest** (currently `3.0.2-debian13`) in line with the pinning posture of
[ADR-0004](0004-secrets-and-supply-chain.md), with the GitHub CLI Feature on top. The image's own
CMake already satisfies the build requirement (≥ 3.21), so no separate CMake Feature is layered.
We keep the C/C++ toolchain decided by the existing build: CMake (≥ 3.21), C++17, bundled SQLite,
Boost for geometry support, optional ICU. Everything else in
[ADR-0002](0002-dev-container-runtime.md) stays binding — in particular, the host Docker socket is
never mounted into the container. This ADR supersedes ADR-0002 only in the choice of base image
and preinstalled toolchain.

## Alternatives considered

- **Keep `base:debian` and install the C++ toolchain via Features or `postCreateCommand`** —
  slower container builds and more moving parts than the purpose-built C++ image, for no gain.
- **A custom Dockerfile** — more control, but maintenance the prebuilt image already covers.
- **Treat the existing devcontainer as a silent exception to ADR-0002** — violates the ADR rules
  in [`AGENTS.md`](../../AGENTS.md) §3; a deviation from an accepted ADR needs a superseding ADR.

## Sources / Prior art

- Dev Container C++ image — <https://github.com/devcontainers/images/tree/main/src/cpp>.
- Dev Container Features — <https://containers.dev/features>.
- The existing build configuration in [`CMakeLists.txt`](../../CMakeLists.txt) and the CI
  workflows under `.github/workflows/`, which already exercise this toolchain on all target
  platforms.

## Consequences

- Positive: the container is usable for C/C++ work immediately; the security posture of
  [ADR-0002](0002-dev-container-runtime.md) (no host daemon access) is preserved and now
  explicitly reaffirmed for this project.
- Negative / trade-offs: the digest pin makes image updates a manual step (Dependabot does not
  bump Dev Container base images) — the tag and digest must be refreshed deliberately, which is
  exactly the announced, reviewable update the pin exists to force.
- Follow-ups: on acceptance, a human marks the base-image aspect of ADR-0002 as superseded by
  this ADR in the index.

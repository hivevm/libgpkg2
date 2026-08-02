[![Build Geopackage](https://github.com/hivevm/libgpkg2/actions/workflows/build.yaml/badge.svg)](https://github.com/hivevm/libgpkg2/actions/workflows/build.yaml)
[![Checks](https://github.com/hivevm/libgpkg2/actions/workflows/checks.yml/badge.svg)](https://github.com/hivevm/libgpkg2/actions/workflows/checks.yml)

# libgpkg

A SQLite 3 extension that provides a minimal [OGC GeoPackage](https://www.ogc.org/standards/geopackage)
implementation.

GeoPackage is an open, standards-based, application and platform independent, and self-describing
file format for geodata based on SQLite.

The project was originally started by [Luciad](http://www.luciad.com), but hasn't been developed
for years.

> For agent instructions, see [`AGENTS.md`](AGENTS.md) — the single source of truth for all coding
> agents. The work is driven by a written **specification**
> ([`docs/SPECIFICATION.md`](docs/SPECIFICATION.md)) and **Architecture Decision Records**
> ([`docs/adr/`](docs/adr/)), so intent and the reasoning behind every structural choice stay
> explicit and reviewable. The conformance to the OGC GeoPackage standard is tracked in
> [`docs/CONFORMANCE.md`](docs/CONFORMANCE.md).

## Prerequisites

- [VS Code](https://code.visualstudio.com/) with the
  [Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
  extension — or any DevContainer-compatible IDE
- Docker / Podman (rootless) available on the host

## Getting Started

1. Open the repository in VS Code and choose **Reopen in Container** — the Dev Container with the
   C/C++ toolchain and the preconfigured agent extensions builds automatically.
2. Authenticate your coding agent inside the container (for Claude Code: `claude login`).
3. Start working with the agent — drive the work from the specification and the ADRs.

## Build, Test & Run

The build uses [CMake](https://cmake.org) 3.21 or newer and requires the
[Boost libraries](boost.md) for geometry support.

- **Build:** `cmake -S . -B build && cmake --build build`
- **Test:** `cmake -S . -B build -DGPKG_TEST=ON && cmake --build build && ctest --test-dir build`
- **Run:** `./build/shell/geopackage` — a modified SQLite 3 command-line shell that autoloads the
  GeoPackage extension (statically linked with SQLite 3 and the extension)

The same steps are available as CMake presets (`default`, `test`, `debug` — the latter for use
with the LLDB launch configuration in `.vscode/launch.json`), e.g.
`cmake --preset test && cmake --build --preset test && ctest --preset test`.

The build produces two main binaries:

- `build/shell/geopackage`: the standalone shell described above.
- `build/gpkg/shared/libgpkg.so` (or `gpkg.dll` on Windows): a dynamically loadable SQLite 3 extension
  that provides the GeoPackage functionality. This extension library can be used with any
  SQLite 3 that supports extension loading.

## Usage

libgpkg can be loaded into SQLite using the
[sqlite3\_load\_extension](http://sqlite.org/c3ref/load_extension.html) C function or using the
[load\_extension](http://sqlite.org/lang_corefunc.html#load_extension) SQL function. Once loaded,
libgpkg extends SQLite with GeoPackage functions that can be used just like any of the core
functions that SQLite provides.

libgpkg exposes the _init_geopackage_extension_ to load the geopackage implementation through the
boost geometries.

Example session with the bundled shell:

~~~bash
./build/shell/geopackage

.open [PATH_TO_DATABASE]

# Show all tables
.tables

# Query a spatial table
SELECT ST_ASTEXT(geom) FROM [TABLE_NAME]
~~~

## Dependencies

- libgpkg requires SQLite 3.51.0 or higher.
- libgpkg requires Boost 1.90 or higher.
- Optional: ICU for international collation support (`-DGPKG_ICU=ON`).

## Project Layout

```
README.md             # overview & setup for humans
AGENTS.md             # single source of truth for coding agents
docs/SPECIFICATION.md # the specification: problem, goals, vocabulary
docs/CONFORMANCE.md   # conformance to the OGC GeoPackage standard
docs/adr/             # Architecture Decision Records (+ template)
scripts/              # repository consistency checks, enforced in CI
gpkg/                 # the GeoPackage SQLite extension (C/C++)
sqlite/               # bundled SQLite
shell/                # the gpkg command-line shell
boost/                # Boost libraries for geometry support
icu/                  # optional ICU support
test/                 # tests
.github/              # CI workflows, issue & pull request templates
.devcontainer/        # Dev Container definition
.vscode/              # shared editor settings
.claude/CLAUDE.md     # pointer for Claude Code to read AGENTS.md
.claude/settings.json # Claude Code permissions: prompt before git/gh writes
```

## Dev Container

The environment is defined entirely in
[`.devcontainer/devcontainer.json`](.devcontainer/devcontainer.json): it starts from the prebuilt
C++ dev container image and layers Dev Container Features and VS Code extensions on top.

The Dev Container deliberately has **no access to the host Docker daemon** — the socket is not
mounted ([ADR-0002](docs/adr/0002-dev-container-runtime.md)). To manage the host's containers from
VS Code, run the **Container Tools** extension (`ms-azuretools.vscode-containers`) on the **host**
side: install it in your host VS Code. [`.vscode/settings.json`](.vscode/settings.json) already
pins it to run locally via `remote.extensionKind`, so it keeps talking to the host engine even
when this folder is reopened in the container.

## Coding Agents

The Dev Container preinstalls the **Claude Code** VS Code extension (see
[`.devcontainer/devcontainer.json`](.devcontainer/devcontainer.json)); other agents (OpenAI Codex,
Cursor, OpenCode, GitHub Copilot) work too once you add them. Authenticate your agent inside the
container (for Claude Code: `claude login`).

The rules every agent follows live in [`AGENTS.md`](AGENTS.md); how each agent is wired to read
them is recorded in [ADR-0001](docs/adr/0001-agent-governance-model.md).

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the workflow (specification- and ADR-driven, small
reviewable changes) and [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md) for the community standards we
expect of everyone taking part. Security issues: please follow [`SECURITY.md`](SECURITY.md)
instead of opening a public issue.

## License

libgpkg is distributed under the
[Apache Software License](https://www.apache.org/licenses/LICENSE-2.0) version 2.0 — see
[`LICENSE`](LICENSE).

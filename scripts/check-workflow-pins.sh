#!/usr/bin/env bash
#
# GitHub Actions pinning checks for this template.
#
# Enforces — as CI, not just convention — the supply-chain rule from AGENTS.md §7 and ADR-0004:
# every action used in .github/workflows/ is pinned to an immutable reference, because tags are
# mutable and a repointed tag executes attacker code with the repository's CI permissions.
#
# For every non-comment 'uses:' line in .github/workflows/*.yml|*.yaml:
#   - './…' (local actions): exempt — they ship with the repository itself.
#   - 'docker://IMAGE': must carry a '@sha256:' digest.
#   - anything else: the ref must be a full 40-hex commit SHA, and the line must carry a
#     trailing version comment, e.g.  uses: actions/checkout@3d3c42e… # v7.0.1
#
# Pure bash + coreutils/grep/sed only — present in the Dev Container base image.
#
# Usage:
#     scripts/check-workflow-pins.sh        (or: bash scripts/check-workflow-pins.sh)
# Exit code 0 when all checks pass, 1 otherwise.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKFLOW_DIR="$ROOT/.github/workflows"

errors=()
add_error() { errors+=("$1"); }

check_workflow_pins() {
  local f rel lineno line uses
  while IFS= read -r f; do
    rel="${f#"$ROOT"/}"
    while IFS=: read -r lineno line; do
      # Skip commented lines (first non-blank character is '#') — e.g. the inert ci.yml skeleton.
      [[ "$line" =~ ^[[:space:]]*# ]] && continue
      # Extract the value after 'uses:', stripping quotes.
      uses="$(printf '%s' "$line" | sed -E "s/^[[:space:]-]*uses:[[:space:]]*//; s/^[\"']//; s/[\"'][[:space:]]*$//")"
      case "$uses" in
        ./*)
          continue ;;
        docker://*)
          if [[ "$uses" != *@sha256:* ]]; then
            add_error "$rel:$lineno: docker image is not digest-pinned ('@sha256:…'): $uses"
          fi ;;
        *)
          if [[ ! "$uses" =~ @[0-9a-f]{40}([\"\'[:space:]]|$) ]]; then
            add_error "$rel:$lineno: action is not pinned to a full 40-hex commit SHA: $uses"
          elif [[ ! "$line" =~ \#[[:space:]]*v?[0-9] ]]; then
            add_error "$rel:$lineno: SHA-pinned action is missing its trailing version comment (e.g. '# v7.0.1')"
          fi ;;
      esac
    done < <(grep -nE '^[[:space:]#-]*uses:' "$f")
  done < <(find "$WORKFLOW_DIR" -maxdepth 1 -type f \( -name '*.yml' -o -name '*.yaml' \) | sort)
}

check_workflow_pins

if ((${#errors[@]} > 0)); then
  echo "Workflow pin checks FAILED:"
  echo
  for e in "${errors[@]}"; do echo "  - $e"; done
  echo
  echo "${#errors[@]} problem(s) found."
  exit 1
fi

echo "Workflow pin checks passed."
exit 0

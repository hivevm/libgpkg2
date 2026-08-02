#!/usr/bin/env bash
#
# Git convention checks for this template.
#
# Enforces — as CI, not just convention — the git rules from AGENTS.md §6 and ADR-0003:
#
#   1. Branch name: work branches match 'type/short-topic' (kebab-case). 'main' and bot
#      branches (dependabot/*, renovate/*) are exempt.
#   2. Commit subjects: every commit in the checked range is a Conventional Commit subject
#      (type, optional scope, optional '!', ': description'), at most 72 characters, with no
#      trailing period. 'fixup!'/'squash!' commits and git-native 'Revert "…"' subjects are
#      exempt ('fixup!'/'squash!' disappear in the squash merge anyway).
#   3. PR title (optional third argument): same rule as a commit subject — squash merge makes
#      the PR title the permanent commit subject on main.
#
# Pure bash + coreutils/git only — present in the Dev Container base image.
#
# Usage:
#     scripts/check-git-conventions.sh [<branch-name> [<range> [<pr-title>]]]
#
# Without arguments (local use) the current branch and 'origin/main..HEAD' (fallback
# 'main..HEAD') are checked; on 'main' itself there is nothing to compare, so the script
# skips with exit 0. CI passes branch, range, and PR title explicitly.
# Exit code 0 when all checks pass, 1 otherwise.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

TYPES='feat|fix|docs|refactor|perf|test|build|ci|chore'
BRANCH_RE="^($TYPES)/[a-z0-9]+(-[a-z0-9]+)*$"
SUBJECT_RE="^($TYPES|revert)(\([a-z0-9._/-]+\))?!?: [^ ].*$"
EXEMPT_BRANCH_RE='^(dependabot|renovate)/'
MAX_SUBJECT_LEN=72

errors=()
add_error() { errors+=("$1"); }

BRANCH="${1:-$(git -C "$ROOT" rev-parse --abbrev-ref HEAD)}"
RANGE="${2:-}"
PR_TITLE="${3:-}"

if [[ $# -eq 0 && "$BRANCH" == "main" ]]; then
  echo "On 'main' with no explicit range — nothing to check."
  exit 0
fi

if [[ -z "$RANGE" ]]; then
  if git -C "$ROOT" rev-parse --verify --quiet origin/main >/dev/null; then
    RANGE='origin/main..HEAD'
  else
    RANGE='main..HEAD'
  fi
fi

# check_subject <label> <subject> — validate one subject line (commit or PR title).
check_subject() {
  local label="$1" subject="$2"
  if [[ ! "$subject" =~ $SUBJECT_RE ]]; then
    add_error "$label: subject is not a Conventional Commit ('type(scope)?: description'): '$subject'"
    return
  fi
  if (( ${#subject} > MAX_SUBJECT_LEN )); then
    add_error "$label: subject exceeds $MAX_SUBJECT_LEN characters (${#subject}): '$subject'"
  fi
  if [[ "$subject" == *. ]]; then
    add_error "$label: subject ends with a period: '$subject'"
  fi
}

check_branch() {
  [[ "$BRANCH" == "main" || "$BRANCH" == "HEAD" ]] && return
  [[ "$BRANCH" =~ $EXEMPT_BRANCH_RE ]] && return
  if [[ ! "$BRANCH" =~ $BRANCH_RE ]]; then
    add_error "branch '$BRANCH' does not match 'type/short-topic' (types: ${TYPES//|/ }, kebab-case)"
  fi
}

check_commits() {
  local count=0 sha subject
  # An unresolvable range must fail loudly — a silent 'git log' error would pass 0 commits green.
  if ! git -C "$ROOT" rev-list "$RANGE" >/dev/null 2>&1; then
    add_error "commit range '$RANGE' is not resolvable"
    return
  fi
  while IFS=' ' read -r sha subject; do
    [[ -z "$sha" ]] && continue
    count=$((count + 1))
    case "$subject" in
      'fixup! '*|'squash! '*|'Revert "'*) continue ;;
    esac
    check_subject "commit $sha" "$subject"
  done < <(git -C "$ROOT" log --no-merges --format='%h %s' "$RANGE" 2>/dev/null)
  echo "$count commit(s) checked in range '$RANGE'."
}

check_branch
check_commits
[[ -n "$PR_TITLE" ]] && check_subject "PR title" "$PR_TITLE"

if ((${#errors[@]} > 0)); then
  echo "Git convention checks FAILED:"
  echo
  for e in "${errors[@]}"; do echo "  - $e"; done
  echo
  echo "${#errors[@]} problem(s) found."
  exit 1
fi

echo "Git convention checks passed."
exit 0

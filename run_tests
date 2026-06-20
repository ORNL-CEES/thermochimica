#!/usr/bin/env bash
# Run TestThermo* executables, optionally skip some, colorize output, and summarize.

set -Euo pipefail           # NOTE: no -e so we don't exit on first failure
shopt -s nullglob

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)
cd "$script_dir/bin"

# Pattern of tests to run (glob). Example to narrow: TEST_PATTERN="TestThermo7*"
TEST_PATTERN="${TEST_PATTERN:-TestThermo*}"

# --- Optional skip list ---
# Provide via env var, e.g.:
#   SKIP_TESTS="TestThermo40 TestThermo41,TestThermo44
#               TestThermo47"
# If SKIP_TESTS is empty or unset, NOTHING is skipped.
SKIP_TESTS_RAW="${SKIP_TESTS:-}"
# Normalize separators to spaces -> array
if [[ -n "$SKIP_TESTS_RAW" ]]; then
  # shellcheck disable=SC2206
  SKIP_TESTS_ARR=( $(printf '%s' "$SKIP_TESTS_RAW" | tr ',\n' '  ' | xargs -n1 | paste -sd' ' -) )
else
  SKIP_TESTS_ARR=()
fi

is_skipped() {
  local t="$1"
  # No skip list? Never skip.
  ((${#SKIP_TESTS_ARR[@]}==0)) && return 1
  for s in "${SKIP_TESTS_ARR[@]}"; do
    [[ "$t" == "$s" ]] && return 0
  done
  return 1
}

# --- Colors (auto-disable if not a TTY or NO_COLOR is set) ---
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
  if command -v tput >/dev/null 2>&1 && tput colors >/dev/null 2>&1; then
    RED=$(tput setaf 1); GREEN=$(tput setaf 2); YELLOW=$(tput setaf 3); BOLD=$(tput bold); RESET=$(tput sgr0)
  else
    RED=$'\e[31m'; GREEN=$'\e[32m'; YELLOW=$'\e[33m'; BOLD=$'\e[1m'; RESET=$'\e[0m'
  fi
else
  RED=""; GREEN=""; YELLOW=""; BOLD=""; RESET=""
fi

FAILED=()
PASSED=()
TOTAL=0

run_test() {
  local t="$1"
  local log="${t}.log"
  ./"$t" > "$log" 2>&1
  return $?
}

# Build the test list (no mapfile/readarray dependency)
tests_list="$(
  for p in ./$TEST_PATTERN; do
    [[ -x "$p" ]] || continue
    printf '%s\n' "${p#./}"
  done | sort -V
)"

# Run them
while IFS= read -r t; do
  [[ -n "$t" ]] || continue
  if is_skipped "$t"; then
    printf '%s[SKIP]%s %s\n' "$YELLOW" "$RESET" "$t"
    continue
  fi

  ((TOTAL++))
  printf '[RUN ] %s\n' "$t"

  start_ns=$(date +%s%N 2>/dev/null || echo 0)
  if run_test "$t"; then
    end_ns=$(date +%s%N 2>/dev/null || echo 0)
    if [[ $start_ns != 0 && $end_ns != 0 ]]; then
      dur_ms=$(( (end_ns - start_ns)/1000000 ))
      printf '%s[PASS]%s %s (%d ms)\n' "$GREEN" "$RESET" "$t" "$dur_ms"
    else
      printf '%s[PASS]%s %s\n' "$GREEN" "$RESET" "$t"
    fi
    PASSED+=("$t")
  else
    end_ns=$(date +%s%N 2>/dev/null || echo 0)
    if [[ $start_ns != 0 && $end_ns != 0 ]]; then
      dur_ms=$(( (end_ns - start_ns)/1000000 ))
      printf '%s[FAIL]%s %s (%d ms)\n' "$RED" "$RESET" "$t" "$dur_ms"
    else
      printf '%s[FAIL]%s %s\n' "$RED" "$RESET" "$t"
    fi
    FAILED+=("$t")
  fi
done <<< "$tests_list"

echo
if ((${#FAILED[@]})); then
  printf '%sSummary:%s %d run, %s%d passed%s, %s%d failed%s\n' \
    "$BOLD" "$RESET" "$TOTAL" "$GREEN" "${#PASSED[@]}" "$RESET" "$RED" "${#FAILED[@]}" "$RESET"
  printf '%sFailed:%s %s\n' "$RED" "$RESET" "${FAILED[*]}"
  exit 1
else
  printf '%sSummary:%s %d run, %s%d passed%s, %s%d failed%s\n' \
    "$BOLD" "$RESET" "$TOTAL" "$GREEN" "${#PASSED[@]}" "$RESET" "$RED" "0" "$RESET"
fi

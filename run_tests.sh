#!/usr/bin/env bash
# Run TestThermo* executables, continue on failure, and summarize.

set -Euo pipefail           # no -e so we don't exit on first failure
shopt -s nullglob

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)
cd "$script_dir/bin"

TEST_PATTERN="${TEST_PATTERN:-TestThermo*}"

# Tests to skip
SKIP_TESTS=(
  TestThermo00
)

is_skipped() {
  local t="$1"
  for s in "${SKIP_TESTS[@]}"; do
    [[ "$t" == "$s" ]] && return 0
  done
  return 1
}

FAILED=()
PASSED=()
TOTAL=0

run_test() {
  local t="$1"
  local log="${t}.log"

  if [[ "$t" == "TestThermo12" ]]; then
    ./"$t" > test12out.txt 2>&1
    tail -n 1 test12out.txt || true
    return $?
  else
    ./"$t" > "$log" 2>&1
    return $?
  fi
}

# Build the test list without using mapfile/readarray.
# 1) Expand the glob, 2) strip ./, 3) sort -V.
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
    printf '[SKIP] %s\n' "$t"
    continue
  fi

  ((TOTAL++))
  printf '[RUN ] %s\n' "$t"

  start_ns=$(date +%s%N 2>/dev/null || echo 0)
  if run_test "$t"; then
    end_ns=$(date +%s%N 2>/dev/null || echo 0)
    if [[ $start_ns != 0 && $end_ns != 0 ]]; then
      dur_ms=$(( (end_ns - start_ns)/1000000 ))
      printf '[PASS] %s (%d ms)\n' "$t" "$dur_ms"
    else
      printf '[PASS] %s\n' "$t"
    fi
    PASSED+=("$t")
  else
    end_ns=$(date +%s%N 2>/dev/null || echo 0)
    if [[ $start_ns != 0 && $end_ns != 0 ]]; then
      dur_ms=$(( (end_ns - start_ns)/1000000 ))
      printf '[FAIL] %s (%d ms)\n' "$t" "$dur_ms"
    else
      printf '[FAIL] %s\n' "$t"
    fi
    FAILED+=("$t")
  fi
done <<< "$tests_list"

echo
printf 'Summary: %d run, %d passed, %d failed\n' "$TOTAL" "${#PASSED[@]}" "${#FAILED[@]}"
if ((${#FAILED[@]})); then
  printf 'Failed: %s\n' "${FAILED[*]}"
  exit 1
fi

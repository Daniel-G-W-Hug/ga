#!/usr/bin/env bash
#
# compare_algebras.sh
#
# Compare the (LHS, RHS) operand pair sets of one product across all
# same-dimension algebras. Spots cross-algebra coverage gaps — e.g.,
# pga3dp dot was missing entries that sta4d dot had, before we added them.
#
# Pairs same-dim algebras:
#   2D: ega2d
#   3D: ega3d, pga2dp
#   4D: pga3dp, sta4d
#
# Usage:
#   ./compare_algebras.sh <product>
#
# Example:
#   ./compare_algebras.sh dot
#       compares dot pair coverage across ega3d/pga2dp (3D) and
#       pga3dp/sta4d (4D).
#

set -u

if [ $# -ne 1 ]; then
    echo "usage: $0 <product>" >&2
    exit 1
fi

prod=$1
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
algebra_dir="$(cd "$script_dir/.." && pwd)"

# Group same-dim algebras together. Single-algebra dims are skipped (nothing
# to compare). Add new dims here as they're introduced.
groups="3D:ega3d,pga2dp 4D:pga3dp,sta4d"

extract_pairs() {
    local alg=$1
    # Pre-filter out commented-out placeholder entries (lines starting with
    # `//`) and pre-normalise infix product forms ("X * Y -> Z", etc.) into
    # the functional form so either spelling is recognised.
    { grep -v '^[[:space:]]*//' \
        "$algebra_dir/ga_prdxpr_${alg}_config.cpp" 2>/dev/null \
        | sed -E '
            s/\{"([a-z_]+) \* ([a-z_]+) -> /{"gpr(\1,\2) -> /g
            s/\{"([a-z_]+) \^ ([a-z_]+) -> /{"wdg(\1,\2) -> /g
            s/\{"([a-z_]+) << ([a-z_]+) -> /{"l_contract(\1,\2) -> /g
            s/\{"([a-z_]+) >> ([a-z_]+) -> /{"r_contract(\1,\2) -> /g
          ' \
        | grep -hoE "\\{\"${prod}\\([^\"]+\\) -> [^\"]+\"" || true; } \
        | sed -E "s/.*\"${prod}\\(([^,]+), ?([^)]+)\\).*/\\1|\\2/" \
        | sort -u
}

for group in $groups; do
    dim="${group%%:*}"
    algs_csv="${group#*:}"
    algs="$(echo "$algs_csv" | tr ',' ' ')"

    echo "=== ${dim}: ${algs} (product=${prod}) ==="
    union=$(mktemp)
    : > "$union"

    # Build per-algebra pair files and the union.
    pair_files=""
    for alg in $algs; do
        f=$(mktemp)
        extract_pairs "$alg" > "$f"
        cat "$f" >> "$union"
        pair_files="$pair_files $alg=$f"
    done
    sort -u -o "$union" "$union"
    union_n=$(grep -c . "$union" || true)

    for entry in $pair_files; do
        alg="${entry%%=*}"
        f="${entry#*=}"
        n=$(grep -c . "$f" || true)
        missing=$(comm -23 "$union" "$f")
        if [ -z "$missing" ]; then
            printf "  %-8s %2d/%d pairs (complete)\n" "$alg" "$n" "$union_n"
        else
            printf "  %-8s %2d/%d pairs, missing:\n" "$alg" "$n" "$union_n"
            echo "$missing" | sed 's/^/    /'
        fi
        rm -f "$f"
    done
    rm -f "$union"
    echo
done

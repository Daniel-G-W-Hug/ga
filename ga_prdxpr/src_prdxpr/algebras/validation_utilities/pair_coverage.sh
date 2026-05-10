#!/usr/bin/env bash
#
# pair_coverage.sh
#
# Print a coverage matrix: for each algebra, count the unique (LHS, RHS)
# operand pairs declared for each major product. Used to spot completeness
# gaps where one product has fewer pairs than its peers in the same algebra.
#
# Useful follow-up: for any cell that looks low, run
#   compare_products.sh <product> <reference> <algebra>
# to enumerate the specific missing pairs.
#
# Usage:
#   ./pair_coverage.sh                       # all default products
#   ./pair_coverage.sh gpr wdg dot           # custom product list
#

set -u

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
algebra_dir="$(cd "$script_dir/.." && pwd)"

algebras="ega2d ega3d pga2dp pga3dp sta4d"
if [ $# -eq 0 ]; then
    products="gpr wdg rgpr rwdg dot rdot cmt rcmt"
else
    products="$*"
fi

# Header
printf "%-8s" "algebra"
for p in $products; do printf " %7s" "$p"; done
echo

# Body
for alg in $algebras; do
    f="$algebra_dir/ga_prdxpr_${alg}_config.cpp"
    if [ ! -f "$f" ]; then
        printf "%-8s (config not found)\n" "$alg"
        continue
    fi
    printf "%-8s" "$alg"
    for p in $products; do
        # grep returning 1 (no matches) is expected for missing products;
        # use || true to avoid tripping the pipeline. Pre-filter out
        # commented-out placeholder entries (lines starting with `//`)
        # and pre-normalise infix product forms ("X * Y -> Z", "X ^ Y -> Z",
        # "X << Y -> Z", "X >> Y -> Z") into the functional form so a single
        # match handles both equivalent spellings.
        n=$( { grep -v '^[[:space:]]*//' "$f" \
                | sed -E '
                    s/\{"([a-z_]+) \* ([a-z_]+) -> /{"gpr(\1,\2) -> /g
                    s/\{"([a-z_]+) \^ ([a-z_]+) -> /{"wdg(\1,\2) -> /g
                    s/\{"([a-z_]+) << ([a-z_]+) -> /{"l_contract(\1,\2) -> /g
                    s/\{"([a-z_]+) >> ([a-z_]+) -> /{"r_contract(\1,\2) -> /g
                  ' \
                | grep -hoE "\\{\"${p}\\([^\"]+\\) -> [^\"]+\"" || true; } \
            | sed -E "s/.*\"${p}\\(([^,]+), ?([^)]+)\\).*/\\1|\\2/" \
            | sort -u \
            | grep -c . \
            || true )
        printf " %7s" "$n"
    done
    echo
done

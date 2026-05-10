#!/usr/bin/env bash
#
# compare_products.sh
#
# Compare the (LHS, RHS) operand pair sets of two products in one algebra.
# Lists pairs in <product1> missing from <product2> and vice versa. Used
# when adding a new product family or auditing whether a regressive product
# mirrors its non-regressive counterpart.
#
# Usage:
#   ./compare_products.sh <product1> <product2> <algebra>
#
# Example:
#   ./compare_products.sh gpr wdg pga3dp
#       lists pairs gpr declares but wdg doesn't in pga3dp, and vice versa.
#

set -u

if [ $# -ne 3 ]; then
    echo "usage: $0 <product1> <product2> <algebra>" >&2
    exit 1
fi

prod1=$1
prod2=$2
alg=$3

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config="$(cd "$script_dir/.." && pwd)/ga_prdxpr_${alg}_config.cpp"

if [ ! -f "$config" ]; then
    echo "error: config not found: $config" >&2
    exit 1
fi

# Extract sorted unique pair sets for each product (whitespace-normalized).
# grep returning 1 (no matches) is expected for absent products; tolerate it.
# Pre-filter out commented-out placeholder entries (lines starting with `//`)
# and pre-normalise infix product forms ("X * Y -> Z" → "gpr(X,Y) -> Z"
# etc.) so we recognise either spelling.
extract_pairs() {
    local prod=$1
    { grep -v '^[[:space:]]*//' "$config" \
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

tmp1=$(mktemp); tmp2=$(mktemp)
trap 'rm -f "$tmp1" "$tmp2"' EXIT
extract_pairs "$prod1" > "$tmp1"
extract_pairs "$prod2" > "$tmp2"

n1=$(grep -c . "$tmp1" || true)
n2=$(grep -c . "$tmp2" || true)
echo "${alg}: ${prod1}=${n1} pairs, ${prod2}=${n2} pairs"
echo

missing_in_2=$(comm -23 "$tmp1" "$tmp2")
missing_in_1=$(comm -13 "$tmp1" "$tmp2")

if [ -n "$missing_in_2" ]; then
    echo "in ${prod1}, missing from ${prod2}:"
    echo "$missing_in_2" | sed 's/^/  /'
fi

if [ -n "$missing_in_1" ]; then
    echo "in ${prod2}, missing from ${prod1}:"
    echo "$missing_in_1" | sed 's/^/  /'
fi

if [ -z "$missing_in_2" ] && [ -z "$missing_in_1" ]; then
    echo "pair sets match (${n1} pairs each)"
fi

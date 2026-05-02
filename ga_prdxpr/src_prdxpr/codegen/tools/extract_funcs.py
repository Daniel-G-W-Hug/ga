#!/usr/bin/env python3
"""Extract C++ function bodies (keyed by signature) from a header.

Outputs one file per function: signature line followed by body lines. Comments and
blank lines between functions are stripped.
"""
import re
import sys

def extract(path):
    text = open(path).read()
    # crude: split on "template <typename" markers; keep the chunk that follows.
    chunks = re.split(r'(?=^template <typename)', text, flags=re.MULTILINE)
    out = []
    for c in chunks:
        if not c.lstrip().startswith('template <typename'):
            continue
        # Find end of body: track braces from first '{'.
        first_brace = c.find('{')
        if first_brace < 0:
            continue
        depth = 0
        i = first_brace
        while i < len(c):
            if c[i] == '{':
                depth += 1
            elif c[i] == '}':
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
            i += 1
        else:
            continue
        body = c[:end]
        # Strip trailing spaces on each line
        body = '\n'.join(line.rstrip() for line in body.splitlines())
        # Normalize: collapse multi-line signature to single line if short enough
        out.append(body)
    return out

def normalize(body):
    """Collapse multi-line wraps in the signature line so signatures compare cleanly."""
    # Remove trailing newline
    body = body.rstrip()
    return body

if __name__ == '__main__':
    for f in extract(sys.argv[1]):
        sys.stdout.write(normalize(f) + '\n\n')

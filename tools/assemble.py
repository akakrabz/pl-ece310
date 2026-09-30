#!/usr/bin/env python3
"""Inline figures into the site content, in place.

A line `<!-- fig:NAME -->` in content/**/*.md is replaced by the contents of
tools/figs/NAME.html (a single <figure> block with no blank lines), followed by
a blank line: a <figure> is an HTML block in CommonMark and runs until the next
blank line, so without it the following paragraph/callout would be swallowed.
Run from the quartz/ directory: python3 tools/assemble.py
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CONTENT = os.path.join(os.path.dirname(HERE), "content")
FIGS = os.path.join(HERE, "figs")

FIG_RE = re.compile(r"^<!--[ \t]*fig:([a-z0-9-]+)[ \t]*-->[ \t]*\n(?:[ \t]*\n)?", re.M)

def main():
    used, missing = set(), []
    for root, _, files in os.walk(CONTENT):
        for fn in files:
            if not fn.endswith(".md"):
                continue
            sp = os.path.join(root, fn)
            rel = os.path.relpath(sp, CONTENT)
            text = open(sp, encoding="utf-8").read()
            def sub(m):
                name = m.group(1)
                fp = os.path.join(FIGS, name + ".html")
                if not os.path.exists(fp):
                    missing.append((rel, name)); return f"<!-- MISSING FIGURE {name} -->"
                used.add(name)
                return open(fp, encoding="utf-8").read().rstrip("\n") + "\n\n"
            new = FIG_RE.sub(sub, text)
            if new != text:
                open(sp, "w", encoding="utf-8").write(new)
    print(f"figures inlined: {sorted(used) or 'none'}")
    if missing:
        print("MISSING FIGURES:", missing); sys.exit(1)
    unused = {f[:-5] for f in os.listdir(FIGS) if f.endswith('.html')} - used
    if unused:
        print("unused figures:", sorted(unused))

if __name__ == "__main__":
    main()

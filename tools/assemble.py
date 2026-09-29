#!/usr/bin/env python3
"""Assemble the site content: copy content-src -> site/content, inlining figures.

A line `<!-- fig:NAME -->` in a source page is replaced by the contents of
work/figs/NAME.html (a single <figure> block with no blank lines).
"""
import os, re, shutil, sys

SRC = "/home/claude/work/content-src"
DST = "/home/claude/site/content"
FIGS = "/home/claude/work/figs"

FIG_RE = re.compile(r"^<!--\s*fig:([a-z0-9-]+)\s*-->\s*$", re.M)

def main():
    if os.path.exists(DST):
        shutil.rmtree(DST)
    used, missing = set(), []
    for root, _, files in os.walk(SRC):
        for fn in files:
            if not fn.endswith(".md"):
                continue
            sp = os.path.join(root, fn)
            rel = os.path.relpath(sp, SRC)
            dp = os.path.join(DST, rel)
            os.makedirs(os.path.dirname(dp), exist_ok=True)
            text = open(sp, encoding="utf-8").read()
            def sub(m):
                name = m.group(1)
                fp = os.path.join(FIGS, name + ".html")
                if not os.path.exists(fp):
                    missing.append((rel, name)); return f"<!-- MISSING FIGURE {name} -->"
                used.add(name)
                return open(fp, encoding="utf-8").read().rstrip("\n")
            text = FIG_RE.sub(sub, text)
            open(dp, "w", encoding="utf-8").write(text)
    n = sum(len(f) for _, _, f in os.walk(DST))
    print(f"assembled {n} files into {DST}; figures used: {sorted(used)}")
    if missing:
        print("MISSING FIGURES:", missing); sys.exit(1)
    unused = {f[:-5] for f in os.listdir(FIGS) if f.endswith('.html')} - used
    if unused:
        print("unused figures:", sorted(unused))

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Render one question variant to a standalone HTML page (and optionally a PNG) for a visual check.

    python3 tools/preview.py convolution/finite-convolution --seed 3 --png
    python3 tools/preview.py ztransform/inverse-pfe --seed 5 --submit correct --png

Writes tools/.preview/<id>-<seed>.html (+ .png). The page shows the question panel, then (with
--submit) the submission panel for PrairieLearn's "correct"/"incorrect" test answer, then the answer
panel. Styling: Bootstrap 5 (from the harness cache) + the elements' own CSS; math is typeset with a
local KaTeX (PrairieLearn itself uses MathJax 4, which looks slightly different)."""

from __future__ import annotations

import argparse
import pathlib
import sys

TOOLS = pathlib.Path(__file__).resolve().parent
sys.path[:0] = [str(TOOLS), str(TOOLS.parent / "serverFilesCourse")]
import plsim  # noqa: E402

KATEX_DIRS = [pathlib.Path("/opt/npm-tools/node_modules/katex/dist"),
              pathlib.Path.home() / ".npm-global/lib/node_modules/markdownlint-cli2/node_modules/katex/dist",
              TOOLS / ".cache/katex/dist"]


def page(title: str, sections: list[tuple[str, str]]) -> str:
    css = []
    bs = plsim._CACHE / "bootstrap.min.css"
    if bs.exists():
        css.append(bs.read_text())
    for el in ("pl-number-input", "pl-integer-input", "pl-symbolic-input", "pl-multiple-choice", "pl-checkbox",
               "pl-matrix-input", "pl-string-input", "pl-hidden-hints"):
        for f in (plsim.PL_APP / "elements" / el).glob("*.css"):
            css.append(f.read_text())
    katex = next((d for d in KATEX_DIRS if (d / "katex.min.js").exists()), None)
    head = ""
    if katex:
        head = (f'<link rel="stylesheet" href="file://{katex}/katex.min.css">'
                f'<script src="file://{katex}/katex.min.js"></script>'
                f'<script src="file://{katex}/contrib/auto-render.min.js"></script>')
    body = "".join(f'<div class="card mb-4"><div class="card-header"><b>{name}</b></div><div class="card-body">{html}</div></div>'
                   for name, html in sections)
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>
<style>{''.join(css)}</style>{head}</head><body class="p-3" style="max-width:900px">
<h5 class="mb-3">{title}</h5>{body}
<script>
if (window.renderMathInElement) renderMathInElement(document.body, {{delimiters: [
  {{left: "$$", right: "$$", display: true}}, {{left: "$", right: "$", display: false}},
  {{left: "\\\\(", right: "\\\\)", display: false}}, {{left: "\\\\[", right: "\\\\]", display: true}}],
  macros: {{"\\\\lt": "<", "\\\\gt": ">"}}, throwOnError: false}});
</script></body></html>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("qid")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--submit", choices=["correct", "incorrect"], default=None)
    ap.add_argument("--png", action="store_true", help="also save a full-page screenshot (needs playwright)")
    args = ap.parse_args()
    q = plsim.Question(args.qid)
    v = q.generate(args.seed)
    sections = [("Question panel", v.render("question"))]
    if args.submit:
        exp, sub, _ = v.pl_test(args.submit)
        sections.append((f"Submission panel ({args.submit} test answer, score {sub['score']})", v.render("submission", submission=sub)))
    sections.append(("Answer panel", v.render("answer")))
    out = TOOLS / ".preview" / f"{args.qid.replace('/', '__')}-{args.seed}.html"
    out.parent.mkdir(exist_ok=True)
    out.write_text(page(f"{args.qid} — seed {args.seed}", sections), encoding="utf-8")
    print(out)
    if args.png:
        from playwright.sync_api import sync_playwright
        with sync_playwright() as p:
            b = p.chromium.launch()
            pg = b.new_page(viewport={"width": 960, "height": 900})
            pg.goto(out.as_uri())
            pg.wait_for_timeout(400)
            png = out.with_suffix(".png")
            pg.screenshot(path=str(png), full_page=True)
            b.close()
        print(png)
    for w in sorted(set(q.warnings)):
        print("warn:", w)


if __name__ == "__main__":
    main()

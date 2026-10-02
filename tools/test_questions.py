#!/usr/bin/env python3
"""Offline test runner for this PrairieLearn course (no PrairieLearn server needed).

    python3 tools/test_questions.py                    # all questions
    python3 tools/test_questions.py --only convolution # questions whose id contains "convolution"
    python3 tools/test_questions.py --seeds 20 --gen-seeds 200 -v

What it checks
  course   every info*.json validates against PrairieLearn's own JSON schemas; UUIDs are unique;
           question topics/tags are declared in infoCourse.json; assessment question ids exist.
  template question.html Mustache tags resolve; no broken "{{" survives rendering.
  generate --gen-seeds variants: generate() + element prepare() run, data stays JSON-serializable,
           and the question's independent checker (tools/checks/<id>.py, see below) agrees.
           Also reports how many of the variants are distinct.
  pipeline --seeds variants: render (question/answer/submission panels) with the real element code;
           PrairieLearn's "Test" button for test types correct / incorrect / invalid (the parsed and
           graded raw submission must reproduce the expected scores, and "correct" must score 1);
           the checker's extra submissions (equivalent forms must score 1, typical mistakes 0);
           every $...$ in the rendered HTML compiles with KaTeX.

Checker modules: tools/checks/<question id with "/" replaced by "__">.py, optional functions
    check(params, correct)        -> list of problem strings (independent recomputation)
    submissions(params, correct)  -> list of (overrides, expect): overrides = {answers-name: raw
                                     string (or list for pl-checkbox)}, expect = {answers-name:
                                     1 | 0 | "invalid"}; unlisted inputs get the correct raw answer.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import re
import subprocess
import sys
import time
import traceback
import uuid as uuidlib

TOOLS = pathlib.Path(__file__).resolve().parent
ROOT = TOOLS.parent
sys.path[:0] = [str(TOOLS), str(ROOT / "serverFilesCourse")]

import plsim  # noqa: E402

KATEX_CANDIDATES = [
    pathlib.Path("/opt/npm-tools/node_modules/katex/dist/katex.min.js"),
    pathlib.Path.home() / ".npm-global/lib/node_modules/markdownlint-cli2/node_modules/katex/dist/katex.min.js",
    TOOLS / ".cache/katex/dist/katex.min.js",
]


# ============================================================================ course-level checks
def load_json(p: pathlib.Path):
    return json.loads(p.read_text(encoding="utf-8"))


def schema_validator(name: str):
    import jsonschema
    p = plsim.PL_APP / "src" / "schemas" / "schemas" / f"{name}.json"
    if not p.exists():
        return None
    schema = load_json(p)
    cls = jsonschema.validators.validator_for(schema)
    return cls(schema)


def course_checks(qids: list[str]) -> list[str]:
    problems = []
    uuids: dict[str, str] = {}

    def note_uuid(u, where):
        if not u:
            return
        try:
            uuidlib.UUID(u)
        except ValueError:
            problems.append(f"{where}: invalid uuid {u!r}")
        if u in uuids:
            problems.append(f"{where}: duplicate uuid {u} (also in {uuids[u]})")
        uuids[u] = where

    validators = {k: schema_validator(k) for k in ("infoCourse", "infoCourseInstance", "infoAssessment", "infoQuestion")}

    def validate(kind, path):
        try:
            data = load_json(path)
        except Exception as exc:
            problems.append(f"{path.relative_to(ROOT)}: invalid JSON: {exc}")
            return None
        v = validators.get(kind)
        if v is not None:
            for err in sorted(v.iter_errors(data), key=lambda e: list(e.path)):
                problems.append(f"{path.relative_to(ROOT)}: schema: {'/'.join(map(str, err.path))}: {err.message}")
        return data

    course = validate("infoCourse", ROOT / "infoCourse.json") or {}
    note_uuid(course.get("uuid"), "infoCourse.json")
    topics = {t["name"] for t in course.get("topics", [])}
    tags = {t["name"] for t in course.get("tags", [])}
    sets = {s["name"] for s in course.get("assessmentSets", [])}

    all_q = {str(p.parent.relative_to(ROOT / "questions")) for p in (ROOT / "questions").rglob("info.json")}
    for qid in sorted(all_q):
        info = validate("infoQuestion", ROOT / "questions" / qid / "info.json") or {}
        note_uuid(info.get("uuid"), f"questions/{qid}")
        if info.get("topic") not in topics:
            problems.append(f"questions/{qid}: topic {info.get('topic')!r} not declared in infoCourse.json")
        for t in info.get("tags", []):
            if t not in tags:
                problems.append(f"questions/{qid}: tag {t!r} not declared in infoCourse.json")

    for ci in sorted((ROOT / "courseInstances").glob("*/infoCourseInstance.json")):
        d = validate("infoCourseInstance", ci) or {}
        note_uuid(d.get("uuid"), str(ci.relative_to(ROOT)))
        for ap in sorted(ci.parent.glob("assessments/**/infoAssessment.json")):
            a = validate("infoAssessment", ap) or {}
            where = str(ap.relative_to(ROOT))
            note_uuid(a.get("uuid"), where)
            if a.get("set") not in sets:
                problems.append(f"{where}: set {a.get('set')!r} not declared in infoCourse.json assessmentSets")
            for z in a.get("zones", []):
                for q in z.get("questions", []):
                    ids = [q["id"]] if "id" in q else [alt["id"] for alt in q.get("alternatives", [])]
                    for i in ids:
                        if i not in all_q:
                            problems.append(f"{where}: question id {i!r} does not exist")
    return problems


# ============================================================================ template lint
TAG = re.compile(r"\{\{(\{?)\s*([#/^!&]?)\s*([^}]*?)\s*\}?\}\}")


def lint_template(q: plsim.Question, variants: list[plsim.Variant]) -> list[str]:
    """Static checks of question.html. A {{params.x}} tag is reported only when no tested variant sets x
    (tags inside Mustache sections may legitimately be unset in variants that hide the section)."""
    problems = []
    t = q.template
    for m in re.finditer(r"\{\{\{\s*([^}]*?)\s*\}\}\}", t):
        if not m.group(1).endswith("_html"):
            problems.append(f"question.html: triple-brace {{{{{{{m.group(1)}}}}}}} — only params named *_html may be unescaped")
    for m in TAG.finditer(t):
        sigil, name = m.group(2), m.group(3)
        if sigil in ("!", "/"):
            continue
        root = name.split(".")[0]
        if root not in ("params", "correct_answers", "submitted_answers", "feedback", "format_errors", "partial_scores",
                        "options", "score", "editable", "panel", "variant_seed") and sigil == "" and "." in name:
            problems.append(f"question.html: suspicious mustache tag {{{{{name}}}}}")
        if root == "params" and "." in name and sigil not in ("^", "#"):
            if not any(_resolves(v.params, name.split(".")[1:]) for v in variants):
                problems.append(f"question.html: {{{{{name}}}}} is not set by generate() in any of the "
                                f"{len(variants)} variants tested")
    return sorted(set(problems), key=problems.index)


def _resolves(params: dict, parts: list[str]) -> bool:
    cur = params
    for part in parts:
        if isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            return False
    return True


# ============================================================================ KaTeX
MATH = re.compile(r"\$\$(.+?)\$\$|\$(.+?)\$", re.S)


def extract_math(html: str) -> list[str]:
    import lxml.html
    try:
        frag = lxml.html.fragment_fromstring(html, create_parent="div")
    except Exception:
        return []
    for bad in frag.xpath("//script|//style|//code|//pre"):
        bad.drop_tree()
    text = frag.text_content()
    return [(a or b).strip() for a, b in MATH.findall(text) if (a or b).strip()]


def katex_check(exprs: list[str]) -> list[str]:
    katex = next((p for p in KATEX_CANDIDATES if p.exists()), None)
    if katex is None or not exprs:
        return []
    js = r"""
const katex = require(process.argv[1]);
const exprs = JSON.parse(require('fs').readFileSync(0, 'utf8'));
const macros = {"\\lt": "<", "\\gt": ">"};
const bad = [];
for (const e of exprs) {
  try { katex.renderToString(e, {throwOnError: true, strict: 'ignore', macros: {...macros}}); }
  catch (err) { bad.push([e.slice(0, 160), String(err.message).slice(0, 160)]); }
}
console.log(JSON.stringify(bad));
"""
    r = subprocess.run(["node", "-e", js, str(katex)], input=json.dumps(sorted(set(exprs))), capture_output=True, text=True)
    if r.returncode != 0:
        return [f"katex runner failed: {r.stderr[:300]}"]
    return [f"KaTeX cannot render `{e}`: {msg}" for e, msg in json.loads(r.stdout)]


# ============================================================================ per-question run
def load_checker(qid: str):
    p = TOOLS / "checks" / (qid.replace("/", "__") + ".py")
    if not p.exists():
        return None
    spec = importlib.util.spec_from_file_location("check_" + qid.replace("/", "__").replace("-", "_"), p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def strip_element_keys(params: dict, names: set[str]) -> dict:
    return {k: v for k, v in params.items() if k not in names}


def run_question(qid: str, gen_seeds: list[int], pipe_seeds: list[int], verbose: bool, require_checker: bool) -> tuple[list[str], dict]:
    problems: list[str] = []
    stats = {"variants": 0, "distinct": 0, "tests": 0, "submissions": 0, "math": 0}
    try:
        q = plsim.Question(qid)
    except Exception as exc:
        return [f"cannot load: {exc}\n{traceback.format_exc()}"], stats
    checker = load_checker(qid)
    if checker is None:
        msg = "no checker in tools/checks/ (answers are not independently verified)"
        (problems if require_checker else q.warnings).append(msg)

    seen = set()
    variants = {}
    names: set[str] = set()
    for seed in sorted(set(gen_seeds) | set(pipe_seeds)):
        try:
            v = q.generate(seed)
        except Exception as exc:
            problems.append(f"seed {seed}: generate/prepare failed: {exc}")
            continue
        variants[seed] = v
        stats["variants"] += 1
        if not names:
            names = set(plsim.answers_names(v))
        seen.add(json.dumps(strip_element_keys(v.params, names), sort_keys=True))
        if checker is not None and hasattr(checker, "check"):
            try:
                for pr in checker.check(v.params, v.true_answer) or []:
                    problems.append(f"seed {seed}: check: {pr}")
            except Exception as exc:
                problems.append(f"seed {seed}: checker crashed: {exc}\n{traceback.format_exc()}")
        if len(problems) > 40:
            problems.append("... stopping early (too many problems)")
            return problems, stats
    stats["distinct"] = len(seen)

    if variants:
        problems += lint_template(q, list(variants.values()))
    math_exprs: list[str] = []
    for i, seed in enumerate(pipe_seeds):
        v = variants.get(seed)
        if v is None:
            continue
        try:
            html_q = v.render("question")
            html_a = v.render("answer")
            if "{{" in html_q or "{{" in html_a:
                problems.append(f"seed {seed}: '{{{{' survives rendering (broken Mustache tag?)")
            if i < 8:
                math_exprs += extract_math(html_q) + extract_math(html_a)
        except Exception as exc:
            problems.append(f"seed {seed}: render failed: {exc}")
            continue
        correct_raw = None
        for tt in ("correct", "incorrect", "invalid"):
            try:
                exp, sub, mism = v.pl_test(tt)
                stats["tests"] += 1
                for m in mism:
                    problems.append(f"seed {seed}: PL test '{tt}': {m}")
                if tt == "correct":
                    correct_raw = exp["raw_submitted_answers"]
                    if not sub["gradable"] or (sub["score"] or 0) < 1 - 1e-9:
                        problems.append(f"seed {seed}: the correct answers score {sub['score']} (format errors {sub['format_errors']}); partial {json.dumps(sub['partial_scores'])[:400]}")
                    if i < 3:
                        html_s = v.render("submission", submission=sub)
                        if i < 2:
                            math_exprs += extract_math(html_s)
                elif tt == "incorrect" and sub["gradable"] and (sub["score"] or 0) >= 1:
                    problems.append(f"seed {seed}: PL's 'incorrect' test submission scored 1")
            except Exception as exc:
                problems.append(f"seed {seed}: PL test '{tt}' crashed: {exc}")
        if checker is not None and hasattr(checker, "submissions") and correct_raw is not None:
            try:
                cases = checker.submissions(v.params, v.true_answer) or []
            except Exception as exc:
                problems.append(f"seed {seed}: submissions() crashed: {exc}\n{traceback.format_exc()}")
                cases = []
            for overrides, expect in cases:
                raw = dict(correct_raw)
                raw.update(overrides)
                try:
                    sub = v.submit(raw)
                except Exception as exc:
                    problems.append(f"seed {seed}: submit {overrides} crashed: {exc}")
                    continue
                stats["submissions"] += 1
                for name, want in expect.items():
                    if want == "invalid":
                        if name not in (sub["format_errors"] or {}):
                            problems.append(f"seed {seed}: expected a format error for {name}={overrides.get(name)!r}")
                        continue
                    if not sub["gradable"]:
                        problems.append(f"seed {seed}: submission {overrides} not gradable: {sub['format_errors']}")
                        break
                    got = (sub["partial_scores"] or {}).get(name, {}).get("score")
                    if got is None or abs(float(got) - float(want)) > 1e-9:
                        problems.append(f"seed {seed}: {name}={overrides.get(name)!r} scored {got}, expected {want}")
        if len(problems) > 40:
            problems.append("... stopping early (too many problems)")
            break
    katex_problems = katex_check(math_exprs)
    stats["math"] = len(set(math_exprs))
    problems += katex_problems
    if verbose and q.warnings:
        for w in sorted(set(q.warnings)):
            print("   warn:", w)
    return problems, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="", help="substring of question ids to test")
    ap.add_argument("--seeds", type=int, default=25, help="variants pushed through the full PL pipeline")
    ap.add_argument("--gen-seeds", type=int, default=300, help="variants generated and checked")
    ap.add_argument("--start", type=int, default=1, help="first seed")
    ap.add_argument("--no-course", action="store_true", help="skip course-level JSON checks")
    ap.add_argument("--require-checker", action="store_true", help="fail questions without a checker")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    qids = sorted(str(p.parent.relative_to(ROOT / "questions")) for p in (ROOT / "questions").rglob("info.json"))
    if args.only:
        qids = [q for q in qids if args.only in q]
    total_fail = 0
    if not args.no_course:
        cp = course_checks(qids)
        print(f"course: {'OK' if not cp else str(len(cp)) + ' problem(s)'}")
        for p in cp:
            print("   ", p)
        total_fail += len(cp)
    gen_seeds = list(range(args.start, args.start + args.gen_seeds))
    pipe_seeds = list(range(args.start, args.start + args.seeds))
    for qid in qids:
        t0 = time.time()
        probs, st = run_question(qid, gen_seeds, pipe_seeds, args.verbose, args.require_checker)
        status = "OK  " if not probs else "FAIL"
        print(f"{status} {qid:55s} variants {st['variants']:4d} (distinct {st['distinct']:4d})  PL tests {st['tests']:3d}  "
              f"subs {st['submissions']:3d}  math {st['math']:3d}  {time.time() - t0:5.1f}s")
        for p in probs[:25]:
            print("     -", p.rstrip()[:1500])
        if len(probs) > 25:
            print(f"     ... and {len(probs) - 25} more")
        total_fail += bool(probs)
    print("ALL OK" if total_fail == 0 else f"{total_fail} failing item(s)")
    sys.exit(1 if total_fail else 0)


if __name__ == "__main__":
    main()

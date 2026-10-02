"""Local emulation of PrairieLearn's freeform ("v3") question pipeline.

It runs the REAL PrairieLearn element code (``apps/prairielearn/elements``) and the
``prairielearn`` Python library from a PrairieLearn checkout, in the same order as
``apps/prairielearn/src/question-servers/freeform.ts``:

    generate : server.generate(data)
    prepare  : mustache(question.html) -> elements.prepare -> server.prepare
    render   : mustache -> elements.render -> server.render
    test     : mustache -> elements.test   -> score -> server.test
    parse    : mustache -> elements.parse  -> server.parse
    grade    : mustache -> elements.grade  -> score -> server.grade

plus PrairieLearn's "Test" button logic (``src/lib/question-testing.ts``): test() produces a
raw submission and the expected scores, the raw submission is parsed and graded, and the two
results must agree.

Data is round-tripped through JSON (``allow_nan=False``) between phases, exactly like the
JS <-> Python boundary in PrairieLearn, so non-serializable params or NaN/inf fail here too.

Setup: ``bash tools/setup_harness.sh`` fetches the PrairieLearn element code (pinned commit) and the
pure-Python packages it needs (sympy, mpmath, chevron, text_unidecode, coloraide, a pint stand-in)
into ``tools/.cache/``. Override with env ``PL_APP`` (…/apps/prairielearn of any PrairieLearn
checkout) and ``PL_PYLIB``. numpy, lxml, pandas, networkx, jsonschema must be importable.

Typical use (see tools/test_questions.py for the full runner):

    import plsim
    q = plsim.Question("convolution/finite-convolution")   # id under questions/
    v = q.generate(seed=7)                                   # generate + element prepare
    html = v.render("question")                              # real element HTML
    expected, submission, mismatches = v.pl_test("correct")  # PrairieLearn's "Test" button
    sub = v.submit({"ny": "-1", "y": "[1, 2, 3]"})           # parse + grade a raw submission
    sub["score"], sub["partial_scores"], sub["format_errors"]
"""

from __future__ import annotations

import builtins
import contextlib
import copy
import io
import json
import math
import os
import pathlib
import random
import sys
import traceback
from typing import Any

_CACHE = pathlib.Path(os.environ.get("PL_HARNESS_CACHE", pathlib.Path(__file__).resolve().parent / ".cache"))
PL_APP = pathlib.Path(os.environ.get("PL_APP", _CACHE / "PrairieLearn" / "apps" / "prairielearn"))
PL_PYLIB = os.environ.get("PL_PYLIB", str(_CACHE / "pylib"))
if not (PL_APP / "elements").is_dir():
    raise SystemExit(f"PrairieLearn element code not found at {PL_APP}.\n"
                     f"Run:  bash {pathlib.Path(__file__).resolve().parent / 'setup_harness.sh'}   (or set PL_APP)")
for _p in (PL_PYLIB, str(PL_APP / "python")):
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)

import chevron  # noqa: E402
import chevron.renderer as _chevron_renderer  # noqa: E402
import numpy as np  # noqa: E402
from prairielearn.internal import question_phases as _qp  # noqa: E402
from prairielearn.internal.check_data import check_data as _check_data  # noqa: E402

COURSE_DIR = pathlib.Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------- mustache (as mustache.js)
def _js_str(v: Any) -> str:
    """How mustache.js stringifies a JSON value."""
    if v is True:
        return "true"
    if v is False:
        return "false"
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isfinite(v) and v == int(v) and abs(v) < 1e21:
            return builtins.str(int(v))
        return repr(v)
    return builtins.str(v)


def _mustache_js_escape(s: str) -> str:
    table = {"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
             "/": "&#x2F;", "`": "&#x60;", "=": "&#x3D;"}
    return "".join(table.get(ch, ch) for ch in s)


_chevron_renderer.str = _js_str  # type: ignore[attr-defined]  (module-level lookup shadows the builtin)
_chevron_renderer._html_escape = _mustache_js_escape  # type: ignore[attr-defined]


def mustache_render(template: str, data: dict) -> str:
    return chevron.render(template, data)


# --------------------------------------------------------------------------- elements
def _load_core_elements() -> dict[str, dict]:
    elems = {}
    for d in sorted((PL_APP / "elements").iterdir()):
        info = d / "info.json"
        if info.exists():
            j = json.loads(info.read_text())
            elems[d.name] = {"name": d.name, "controller": j["controller"], "type": "core"}
    return elems


CORE_ELEMENTS = _load_core_elements()


def _load_course_elements(course_dir: pathlib.Path) -> dict[str, dict]:
    """Course-specific elements (<course>/elements/<name>/info.json), as PrairieLearn loads them."""
    elems = {}
    root = course_dir / "elements"
    if root.is_dir():
        for d in sorted(root.iterdir()):
            info = d / "info.json"
            if info.exists():
                j = json.loads(info.read_text())
                elems[d.name] = {"name": d.name, "controller": j["controller"], "type": "course"}
    return elems


def _jsonify(data: dict, allow_nan: bool = False) -> dict:
    return json.loads(json.dumps(data, allow_nan=allow_nan))


class PLError(Exception):
    """A failure PrairieLearn would report as a (fatal) course issue."""


@contextlib.contextmanager
def _sandbox(cwd: pathlib.Path, extra_paths: list[str]):
    saved_path, saved_cwd = list(sys.path), os.getcwd()
    out, err = io.StringIO(), io.StringIO()
    try:
        os.chdir(cwd)
        sys.path[:0] = extra_paths
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            yield out, err
    finally:
        sys.path[:] = saved_path
        os.chdir(saved_cwd)


class Question:
    def __init__(self, qid: str, course_dir: pathlib.Path = COURSE_DIR):
        self.qid = qid
        self.course_dir = pathlib.Path(course_dir)
        self.dir = self.course_dir / "questions" / qid
        self.info = json.loads((self.dir / "info.json").read_text())
        self.partial_credit = self.info.get("partialCredit", True)
        self.template = (self.dir / "question.html").read_text(encoding="utf-8")
        self.warnings: list[str] = []
        self._server: dict | None = None
        if (self.dir / "server.py").exists():
            self._server = {"__file__": str(self.dir / "server.py")}
            code = compile((self.dir / "server.py").read_text(encoding="utf-8"), str(self.dir / "server.py"), "exec")
            with _sandbox(self.dir, self._server_paths()) as (out, err):
                exec(code, self._server)
            self._note_output("import server.py", out, err)

    # ---------------------------------------------------------------- helpers
    def _server_paths(self) -> list[str]:
        return [str(self.dir), str(self.course_dir / "serverFilesCourse"), str(PL_APP / "python")]

    def _note_output(self, what: str, out: io.StringIO, err: io.StringIO) -> None:
        text = (out.getvalue() + err.getvalue()).strip()
        if text:
            self.warnings.append(f"{what}: output logged on console: {text[:300]}")

    def options(self, render: bool = False) -> dict:
        o = {
            "question_path": str(self.dir),
            "client_files_question_path": str(self.dir / "clientFilesQuestion"),
            "client_files_course_path": str(self.course_dir / "clientFilesCourse"),
            "server_files_course_path": str(self.course_dir / "serverFilesCourse"),
            "course_extensions_path": str(self.course_dir / "elementExtensions"),
            "user": {"uid": "student@example.com", "name": "Test Student", "uin": None, "user_id": "1"},
            "group": None,
        }
        if render:
            o.update({
                "client_files_question_url": "/pl/clientFilesQuestion",
                "client_files_course_url": "/pl/clientFilesCourse",
                "client_files_question_dynamic_url": "/pl/generatedFilesQuestion",
                "course_element_files_url": "/pl/course/elements",
                "course_element_extension_files_url": "/pl/course/elementExtensions",
                "submission_files_url": None,
                "external_image_capture_url": None,
                "base_url": "/pl",
                "workspace_url": None,
            })
        return o

    def _call_server(self, phase: str, data: dict, html: str | None = None):
        if self._server is None or phase not in self._server or not callable(self._server[phase]):
            return html if phase == "render" else data
        orig = copy.deepcopy(data)
        args = [data] + ([html] if phase == "render" else [])
        with _sandbox(self.dir, self._server_paths()) as (out, err):
            try:
                ret = self._server[phase](*args)
            except Exception as exc:
                raise PLError(f"server.py {phase}(): {type(exc).__name__}: {exc}\n{traceback.format_exc()}") from exc
        self._note_output(f"server.py {phase}()", out, err)
        if phase == "render":
            return ret
        if ret is not None and ret is not data:
            self.warnings.append(f"server.py {phase}() returned a data object other than the one passed in")
            data = ret
        if phase not in ("render", "file"):
            try:
                _check_data(orig, data, phase)
            except ValueError as exc:
                raise PLError(f"server.py: invalid state after {phase}(): {exc}") from exc
        try:
            return _jsonify(data)
        except (TypeError, ValueError) as exc:
            raise PLError(f"server.py {phase}(): data is not JSON-serializable (allow_nan=False): {exc}") from exc

    def _process_html(self, phase: str, data: dict) -> tuple[dict, str | None]:
        try:
            html = mustache_render(self.template, data)
        except Exception as exc:
            raise PLError(f"question.html mustache error: {exc}") from exc
        if "<markdown" in html:
            self.warnings.append("question.html uses <markdown> (not emulated)")
        context = {"html": html, "elements": {**CORE_ELEMENTS, **_load_course_elements(self.course_dir)},
                   "element_extensions": {}, "course_path": str(self.course_dir)}
        saved_path, saved_cwd = list(sys.path), os.getcwd()
        out, err = io.StringIO(), io.StringIO()
        try:
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                result, _processed = _qp.process(phase, data, context)
        except Exception as exc:
            raise PLError(f"element {phase}(): {type(exc).__name__}: {exc}\n{traceback.format_exc()}") from exc
        finally:
            sys.path[:] = saved_path
            os.chdir(saved_cwd)
        self._note_output(f"elements {phase}()", out, err)
        if phase in ("grade", "test"):
            ps = data.get("partial_scores") or {}
            if self.partial_credit:
                tw = sum(float(v.get("weight", 1)) if isinstance(v, dict) and isinstance(v.get("weight", 1), (int, float)) else 1.0 for v in ps.values())
                tws = sum((float(v.get("weight", 1)) if isinstance(v.get("weight", 1), (int, float)) else 1.0) * (v.get("score", 0) if isinstance(v.get("score"), (int, float)) else 0) for v in ps.values() if isinstance(v, dict))
                data["score"] = tws / (tw if tw else 1)
            else:
                data["score"] = 1 if ps and all(isinstance(v, dict) and (v.get("score") or 0) >= 1 for v in ps.values()) else 0
        try:
            data = _jsonify(data)
        except (TypeError, ValueError) as exc:
            raise PLError(f"elements {phase}(): data is not JSON-serializable: {exc}") from exc
        return data, result

    # ---------------------------------------------------------------- phases
    def generate(self, seed: int) -> "Variant":
        random.seed(seed)
        np.random.seed(seed)
        data = {"params": {}, "correct_answers": {}, "variant_seed": seed, "options": self.options(), "preferences": {}}
        data = self._call_server("generate", _jsonify(data))
        v = Variant(self, seed, data["params"], data["correct_answers"])
        # prepare (same worker in PL: seeded once with the variant seed)
        random.seed(seed)
        np.random.seed(seed)
        pdata = {"params": copy.deepcopy(v.params), "correct_answers": copy.deepcopy(v.true_answer),
                 "variant_seed": seed, "options": self.options(), "preferences": {}, "answers_names": {}}
        pdata, _ = self._process_html("prepare", pdata)
        pdata = self._call_server("prepare", pdata)
        v.params, v.true_answer = pdata["params"], pdata["correct_answers"]
        return v


class Variant:
    def __init__(self, q: Question, seed: int, params: dict, true_answer: dict):
        self.q, self.seed, self.params, self.true_answer = q, seed, params, true_answer

    # ----------------------------------------------------------------- render
    def render(self, panel: str = "question", submission: dict | None = None) -> str:
        s = submission or {}
        data = {
            "params": copy.deepcopy(s.get("params", self.params)),
            "correct_answers": copy.deepcopy(s.get("true_answer", self.true_answer)),
            "submitted_answers": copy.deepcopy(s.get("submitted_answers", {})),
            "format_errors": copy.deepcopy(s.get("format_errors", {})),
            "partial_scores": copy.deepcopy(s.get("partial_scores", {})),
            "score": s.get("score", 0),
            "feedback": copy.deepcopy(s.get("feedback", {})),
            "variant_seed": self.seed,
            "options": self.q.options(render=True),
            "preferences": {},
            "raw_submitted_answers": copy.deepcopy(s.get("raw_submitted_answers", {})),
            "editable": submission is None,
            "manual_grading": False,
            "ai_grading": False,
            "panel": panel,
            "num_valid_submissions": 1 if submission else 0,
            "correct_answer_shown": True,
        }
        data, html = self.q._process_html("render", data)
        html = self.q._call_server("render", data, html)
        return html or ""

    # ----------------------------------------------------------------- submit
    def submit(self, raw: dict) -> dict:
        """Parse + grade a raw submission (dict of answers-name -> string) like PrairieLearn."""
        raw = {k: v for k, v in raw.items()}
        data = {
            "params": copy.deepcopy(self.params),
            "correct_answers": copy.deepcopy(self.true_answer),
            "submitted_answers": copy.deepcopy(raw),
            "feedback": {},
            "format_errors": {},
            "variant_seed": self.seed,
            "options": self.q.options(),
            "preferences": {},
            "raw_submitted_answers": copy.deepcopy(raw),
            "gradable": True,
        }
        data, _ = self.q._process_html("parse", data)
        data = self.q._call_server("parse", data)
        if data["format_errors"]:
            data["gradable"] = False
        sub = {"params": data["params"], "true_answer": data["correct_answers"],
               "submitted_answers": data["submitted_answers"], "format_errors": data["format_errors"],
               "raw_submitted_answers": data["raw_submitted_answers"], "feedback": data["feedback"],
               "gradable": data["gradable"], "partial_scores": None, "score": None}
        if not data["gradable"]:
            return sub
        gdata = {
            "params": copy.deepcopy(sub["params"]),
            "correct_answers": copy.deepcopy(sub["true_answer"]),
            "submitted_answers": copy.deepcopy(sub["submitted_answers"]),
            "format_errors": copy.deepcopy(sub["format_errors"]),
            "partial_scores": {},
            "score": 0,
            "feedback": copy.deepcopy(sub["feedback"]),
            "variant_seed": self.seed,
            "options": self.q.options(),
            "preferences": {},
            "raw_submitted_answers": copy.deepcopy(sub["raw_submitted_answers"]),
            "gradable": True,
        }
        gdata, _ = self.q._process_html("grade", gdata)
        gdata = self.q._call_server("grade", gdata)
        if gdata["format_errors"]:
            gdata["gradable"] = False
        sub.update({"params": gdata["params"], "true_answer": gdata["correct_answers"],
                    "submitted_answers": gdata["submitted_answers"], "format_errors": gdata["format_errors"],
                    "partial_scores": gdata["partial_scores"], "score": gdata["score"],
                    "feedback": gdata["feedback"], "gradable": gdata["gradable"]})
        return sub

    # ----------------------------------------------------------------- PL "Test" button
    def pl_test(self, test_type: str) -> tuple[dict, dict, list[str]]:
        """Return (expected, submission, mismatches) as PrairieLearn's question tester does."""
        data = {
            "params": copy.deepcopy(self.params),
            "correct_answers": copy.deepcopy(self.true_answer),
            "format_errors": {},
            "partial_scores": {},
            "score": 0,
            "feedback": {},
            "variant_seed": self.seed,
            "options": self.q.options(),
            "preferences": {},
            "raw_submitted_answers": {},
            "gradable": True,
            "test_type": test_type,
        }
        data, _ = self.q._process_html("test", data)
        data = self.q._call_server("test", data)
        if data["format_errors"]:
            data["gradable"] = False
        expected = data
        sub = self.submit(expected["raw_submitted_answers"])
        problems = []

        def norm(x):  # JSON.stringify has no int/float distinction
            if isinstance(x, bool) or x is None or isinstance(x, str):
                return x
            if isinstance(x, (int, float)):
                return float(x)
            if isinstance(x, dict):
                return {k: norm(v) for k, v in x.items()}
            if isinstance(x, list):
                return [norm(v) for v in x]
            return x

        def chk(name, a, b):
            if json.dumps(norm(a), sort_keys=True) != json.dumps(norm(b), sort_keys=True):
                problems.append(f'"{name}" mismatch: expected {json.dumps(a)[:300]} but got {json.dumps(b)[:300]}')

        chk("gradable", expected["gradable"], sub["gradable"])
        chk("format_errors keys", sorted(expected["format_errors"]), sorted(sub["format_errors"]))
        if expected["gradable"] and sub["gradable"]:
            chk("partial_scores", expected["partial_scores"], sub["partial_scores"])
            chk("score", expected["score"], sub["score"])
        return expected, sub, problems


def answers_names(variant: Variant) -> list[str]:
    """answers-name attributes of all input elements in the rendered question panel template."""
    import lxml.html
    html = mustache_render(variant.q.template, {"params": variant.params, "correct_answers": variant.true_answer})
    frag = lxml.html.fragment_fromstring(html, create_parent="div")
    return [el.get("answers-name") for el in frag.iter() if isinstance(el.tag, str) and el.get("answers-name")]

# Writing and testing questions for this course

This is the working guide for adding randomized questions to the ECE 310 PrairieLearn course —
for people and for Claude sessions alike. Everything here is tested offline with PrairieLearn's
**own** element code, so a question that passes `tools/test_questions.py` behaves the same in
PrairieLearn.

## 1. One-time setup (≈15 s, needs git + python3 ≥ 3.12)

```bash
bash tools/setup_harness.sh        # fetches PrairieLearn element code (pinned) + sympy etc. into tools/.cache/
python3 tools/test_lib.py          # self-test of the shared math library (serverFilesCourse/ece310)
python3 tools/test_questions.py    # every question: schemas, generate, render, PL "Test" button, checkers
```

`numpy`, `scipy`, `lxml`, `pandas`, `networkx`, `jsonschema` must be importable (`pip install …`).
Previews also use `playwright` (optional).

The whole suite takes ~15 min on one core (some checkers are slow on purpose: exact sums, root
grids). Run the topic folders in parallel when you need everything:

```bash
for t in signals systems convolution lccde ztransform stability dtft; do
  python3 tools/test_questions.py --only "$t/" --no-course > "/tmp/pl-$t.log" 2>&1 &
done; wait; tail -n 2 /tmp/pl-*.log; python3 tools/test_questions.py --only none   # course checks
```

## 2. The tools

| file | what it does |
|---|---|
| `tools/setup_harness.sh` | downloads PrairieLearn (sparse, pinned commit; `PL_REF=master` for the newest) and the pure-Python packages it needs into `tools/.cache/` |
| `tools/plsim.py` | runs one question through PrairieLearn's pipeline: `generate → prepare → render → test → parse → grade`, with the real element code, Mustache and JSON round-trips (see its docstring for the API) |
| `tools/test_questions.py` | the runner: course JSON vs PrairieLearn's schemas, template lint, `--gen-seeds` variants through `generate` + checker, `--seeds` variants through render + PrairieLearn's **Test** button (correct / incorrect / invalid) + extra submissions + KaTeX compile of every `$…$` |
| `tools/checks/<id with / → __>.py` | independent answer check for one question (`check`, `submissions`, see §5) |
| `tools/checklib.py` | helpers for checkers: SymPy parsing, `zsum` (overflow-safe truncated z-transform sums), `impulse_response` (lfilter) |
| `tools/preview.py` | renders a variant (question, submission, answer panels) to HTML/PNG for a visual check |
| `tools/test_lib.py` | self-test of `serverFilesCourse/ece310` against scipy / truncated sums |

```bash
python3 tools/test_questions.py --only ztransform/inverse-pfe --seeds 25 --gen-seeds 300 -v
python3 tools/preview.py ztransform/inverse-pfe --seed 4 --submit correct --png
```

The harness is course-independent: copy `tools/plsim.py`, `tools/test_questions.py`,
`tools/checklib.py`, `tools/preview.py`, `tools/setup_harness.sh` into any PrairieLearn course.

## 3. Anatomy of a question

```
questions/<topic-folder>/<slug>/
  info.json      uuid (fresh uuid4), title, topic (from infoCourse.json), tags (only declared ones), "type": "v3"
  server.py      generate(data) — all randomness here; optional grade(data) for custom grading
  question.html  Mustache template + PrairieLearn elements
tools/checks/<topic-folder>__<slug>.py   independent checker
```

Topic folders: `signals`, `systems`, `convolution`, `lccde`, `ztransform`, `stability`, `dtft`.
Topics and tags are declared in `infoCourse.json` (tags: `L1`…`L14`, `hw1`…`hw5`, `midterm1`,
`fam-*` exam families, input types `numeric integer symbolic MC checkbox matrix`, `easy medium hard`).

### server.py rules

- PrairieLearn seeds `random` and `numpy.random` with the variant seed before `generate`; use them
  directly. Never seed yourself, never print (console output is a course issue).
- Write LaTeX in raw strings (`r"\frac{1}{2}"`): `\right`, `\frac`, `\tfrac`, `\beta` start with the
  escapes `\r`, `\f`, `\t`, `\b` and silently become control characters otherwise (this also bites
  Python scripts that edit question files).
- A custom `grade(data)` that compares numbers must treat NaN/inf as wrong — write
  `if not (abs(err) <= tol): …wrong…`, never `if abs(err) > tol` (NaN makes that False, and `0/0` or
  `zoo` scored 100 % once) — then recompute the score with `pl.set_weighted_score_data(data)`.
- `from ece310 import fmt, poly, zt, seq, mc` — `serverFilesCourse/` is on `sys.path`. Do exact
  arithmetic with `fractions.Fraction`, then convert: `params` and `correct_answers` must be
  **plain JSON** (str, int, float, bool, list, dict) — no `Fraction`, no numpy scalars, no NaN/inf.
- Build every piece of LaTeX in `server.py` (`fmt.tex_*`, `zt.ROC.tex()`, `zt.term_*`…) and put
  it in `params["..._tex"]`. Generated LaTeX never contains raw `<`, `>`, `&` (use `\lt \gt \le \ge`).
- Correct answers, per element:
  - `pl-number-input` → `float`. Exact (rational) answers: `comparison="relabs" rtol="1e-3" atol="1e-6"`
    and, next to the inputs, the sentence *"Enter exact values: fractions such as -7/3 are accepted;
    decimals need 4 significant digits."* (with `rtol="1e-6"` a student's `1.667` for 5/3 scores 0).
    Wrong-answer probes in the checker must then differ from the answer by more than ~1 %.
  - `pl-integer-input` → `int`
  - `pl-symbolic-input` → a SymPy string from `fmt` (`**`, rationals as `(3/4)`, never floats), and the element's `variables="..."` must list every symbol (`z`, `n`, `w`…).
    **Answers that are functions of `n`: store `answers.symbolic(sym_string)`** (n declared an integer),
    otherwise correct forms such as `(-3)^(-n)` for `(-1/3)^n` or `(-1)^n 2^n` for `(-2)^n` score 0.
    The element refuses decimals ("enter fractions, not decimals" belongs in the prompt) and students
    cannot type `u[n]`: ask for the expression "for $n \ge n_0$" (and "for $n \le -1$") instead.
  - Equivalences SymPy cannot prove (exponential vs. cosine form of a DTFT, long sums of shifted
    exponentials) or that time out: grade numerically in `server.py` — either leave
    `correct_answers[name]` unset and implement `grade(data)` + `test(data)` (see
    `dtft/dtft-of-finite-sequence`), or keep it and override the score in `grade(data)` after
    evaluating at sample points (see `lccde/system-response`). Recompute `data["score"]` afterwards
    (`pl.set_weighted_score_data(data)`).
  - `pl-matrix-input` → `pl.to_json(np.array([[...]]))` (a 1×N row for sequences). Entries must be
    numbers — `1/2` is a format error — so keep vector answers integer (or short decimals, and say so).
  - `pl-number-input` with `allow-blank="true" blank-value=""` and a correct answer of `""` gives a
    "leave empty if …" input (see `stability/bibo-from-impulse-response`).
  - `pl-multiple-choice` / `pl-checkbox` → **don't** set `correct_answers`; put the options in
    `params` as `[{"text": ..., "correct": true/false}, …]` (`mc.choices`, `mc.yes_no`) and loop over
    them in the template; the element computes the correct answer in its `prepare`. A `pl-checkbox`
    needs at least one correct option, even with `allow-blank="true"`.
- No giveaways: an option list must not reveal a count or a magnitude asked in another part (add
  invalid options, write options in terms of $p_1, p_2$), an input must not appear only in the variants
  where it applies, and a format example must not have the shape of an answer.
- Design variants so the answer is unambiguous (e.g. nonzero end samples in convolutions, distinct
  pole magnitudes when a question needs distinct ROCs) and so the numbers stay "exam nice".

### question.html rules

- Inside math write `{{params.x_tex}}` (double braces). Triple braces only for params named
  `*_html` that hold real HTML (tables). Never put `{{` inside static LaTeX (`x^{{2}}` breaks
  Mustache) — build that LaTeX in `server.py` instead.
- Math delimiters: `$…$` and `$$…$$` (PrairieLearn uses MathJax 4).
- Input `label`/`suffix`: all math or all text. The label box is a flex container, so text next to
  math loses its spaces ("for $n$" shows as "forn"); write `label="$A_1\;(\text{pole } p_1 = \tfrac12)$"`.
- Several templates in one question: wrap each in a Mustache section (`{{#params.kind_a}}…{{/params.kind_a}}`,
  elements included). A `{{params.x}}` used only inside a section needs to be set only in the
  variants that show it (the runner's lint looks at the union over all tested variants).
- Elements are validated against their schemas — an unknown attribute is a fatal error. Read the
  attribute table in `tools/.cache/PrairieLearn/docs/elements/<element>.md` before using one
  (e.g. `pl-checkbox` uses `partial-credit="each-answer"` and `allow-blank`, not the old
  `partial-credit="true"`/`allow-none-correct`).
- Multiple choice from params:
  ```html
  <pl-multiple-choice answers-name="roc">
    {{#params.roc_choices}}<pl-answer correct="{{correct}}">{{text}}</pl-answer>{{/params.roc_choices}}
  </pl-multiple-choice>
  ```
  (`order="fixed"` for Yes/No, True/False or naturally ordered lists; `hide-letter-keys="true"`
  whenever the question's parts are labelled (a), (b), …, otherwise the submission panel reads
  "(b) Is it causal? (b) No".)
- `pl-checkbox`: prefer `partial-credit="net-correct"` or `"coverage"`; with `"each-answer"` a blank
  submission earns credit for every correctly unticked option (74 % on average in one question).
- Two or more input elements in a row: wrap each in `<div class="mb-1">…</div>` (outside a `<p>`), or
  the answer panel shows their answers run together on one line.
- Every question gets `<pl-hidden-hints>` (one free hint, one after a submission) and a
  `<pl-answer-panel>` with the full worked solution, written from generated LaTeX.

### Course notation (match the lectures and exams)

`x[n]`, `h[n]`, `y[n]`, `δ[n]`, `u[n]`; sequences with an arrow under the `n = 0` sample
(`fmt.tex_seq`); transforms in powers of $z^{-1}$, always with the ROC written `|z| > a`;
LCCDE in the Lecture 9 form $y[n] + \sum_k a_k y[n-k] = \sum_k b_k x[n-k]$ (same `b, a` as
`scipy.signal.lfilter`); "LTI" (older exams say "LSI"); DTFT written $X_d(\omega)$.

## 4. The shared library `serverFilesCourse/ece310`

| module | main functions |
|---|---|
| `fmt` | `Q`, `tex_num`, `sym_num`, `plain_num`, `tex_sum`, `sym_sum`, `tex_pow`, `sym_pow`, `zpow_tex`, `tex_poly_zinv`, `sym_poly_zinv`, `tex_factor`, `sym_factor`, `tex_frac`, `tex_seq`, `plain_seq` |
| `poly` | polynomials in $z^{-1}$: `pmul`, `padd`, `pscale`, `from_roots`, `pdivmod` (Lecture 10 long division), `peval`, `trim` |
| `zt` | `ROC` (`right/left/finite`, `intersect`, `contains_unit_circle`, `is_causal_system`, `kind`, `tex`, `to_json`), `rocs_for_poles`, `pfe` (exact cover-up incl. improper part), `inverse_pfe` → `InverseZ` (`value(n)`, `sym_nonneg`, `sym_neg`, `tex`), signal terms `term_right_exp`, `term_left_exp`, `term_n_exp`, `term_finite`, `term_cos` (each with `x_tex`, `X_tex`, `X_sym`, `roc`, `poles`, `spec`), `sum_terms`, `x_value` |
| `seq` | `conv` (with start indices), `trim_seq`, `value` |
| `mc` | `choices` (correct + distinct distractors, shuffled), `yes_no`, `true_false` |
| `answers` | `symbolic(expr, integer=("n",))` — correct answer for `pl-symbolic-input` in `n` (see §3) |
| `convlib`, `zinv`, `lccde`, `systems_bank` | helpers written for particular questions (convolution closed forms; inverse-z formatting; LCCDE/H(z) formatting and exact recursion; the vetted system-property bank with LaTeX, flags and reasons) — reuse them, read their docstrings |

Treat the library as read-only from a question; if a helper is missing, write it in your
`server.py` (or add a new module, never change the behaviour of an existing function).

## 5. Checkers and the definition of done

`tools/checks/<id>.py` recomputes the answer **independently** (numpy / scipy / SymPy, not the
code path `server.py` used) and probes the grading:

```python
def check(params, correct):          # -> list of problem strings ([] = OK)
    ...
def submissions(params, correct):    # -> [(overrides, expect), ...]
    # overrides: {answers-name: raw string (list for pl-checkbox keys)}; others get the correct answer
    # expect:    {answers-name: 1 | 0 | "invalid" | exact partial score, e.g. 0.8 for pl-checkbox each-answer}
    return [({"X": "z/(z-1/2)"}, {"X": 1}),          # equivalent form must score 1
            ({"X": "1/(1+(1/2)*z^(-1))"}, {"X": 0})]  # classic mistake must score 0
```

For `pl-multiple-choice`/`pl-checkbox`, `params[name]` is the element's option list
(`[{"key": "a", "html": …}, …]`) and `correct[name]` is `{"key", "html"}` (a list for checkbox):
submit a key (or a list of keys) and find options by their `html`.

A question is done when
1. `python3 tools/test_questions.py --only <id> --seeds 25 --gen-seeds 300 -v` prints `OK` with no warnings;
2. at least ~90 % of the variants are distinct (the runner prints the count);
3. the checker covers every answer and at least two equivalent forms and two classic mistakes;
4. you looked at `tools/preview.py <id> --seed … --submit correct --png` for two seeds.

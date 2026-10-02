# ECE 310 — PrairieLearn course (randomized practice)

29 randomized PrairieLearn "v3" questions for ECE 310 (Digital Signal Processing): everything on
Midterm 1 (Lectures 1–11, Homework 1–4) plus a short DTFT unit (Lectures 13–14, Homework 5).
Every variant is a new random problem — new signals, poles, coefficients — with "exam-nice"
numbers, an autograded answer and a worked solution in the answer panel that names the lecture,
homework or past-exam problem it is modelled on. Each question has an independent checker that
recomputes its answers with numpy/scipy/SymPy, and the whole course is tested offline with
PrairieLearn's own element code (see [Testing](#testing-without-a-prairielearn-server)).

```
infoCourse.json                       topics, tags (lectures L1–L14, hw1–hw5, exam families fam-*), assessment sets
courseInstances/Fa26/
  infoCourseInstance.json             publishing dates
  assessments/
    PR1-signals-systems-convolution/  Practice 1  (L1–L5)
    PR2-z-transform/                  Practice 2  (L6–L8)
    PR3-transfer-functions-stability/ Practice 3  (L9–L11)
    PR4-dtft/                         Practice 4  (L13–L14, not on Midterm 1)
    ME1-mock-midterm1/                Mock Midterm 1: 10 problems in past-exam shape, 2 h, retakeable
questions/<topic>/<slug>/             info.json, server.py, question.html
serverFilesCourse/ece310/             shared exact-arithmetic library (fractions, polynomials in z^-1, ROCs, PFE …)
tools/                                offline test harness (reusable in any PrairieLearn course)
AUTHORING.md                          how to write and test new questions — start here
```

## Questions

| id | what a variant asks | inputs |
|---|---|---|
| `signals/index-transformations` | y[n] = x[an + b] for a finite x (shift, reversal, decimation): first index and samples | integer, row vector |
| `signals/complex-numbers` | magnitude and principal angle of products/quotients/powers; roots of z^N = c; Re/Im of sums of complex exponentials | numbers |
| `systems/system-properties` | linear? time-invariant? causal? BIBO stable? — 47 vetted exam/HW systems plus a generator of combinations | 4 × Yes/No |
| `systems/true-false` | six statements from the Midterm 1 true/false bank (parameterized ones recomputed) | 6 × True/False |
| `convolution/finite-convolution` | finite x * h with start indices (shift-and-add table in the solution) | integer, row vector |
| `convolution/infinite-convolution` | aⁿu[n−k₁] * bⁿu[n−k₂] (incl. a = b) or finite * exponential: first index, closed form, one sample | integer, symbolic, number |
| `convolution/find-h-deconvolution` | find h[n] from one input–output pair (finite deconvolution or h[n] = y[n] − a·y[n−1]) | integer, row vector |
| `convolution/step-response` | step response from h (running sum) or h from a step response (first difference) | row vector, numbers, symbolic |
| `convolution/lti-superposition` | response to αx₁[n−k] + βx₁[n−m] given x₁ ↦ y₁ (sometimes x₂ only as a sequence) | integer, row vector |
| `lccde/recursion` | h[0..3] of a first/second-order difference equation at rest; FIR or IIR | numbers, MC |
| `lccde/lccde-to-transfer-function` | difference equation → H(z), or H(z) → coefficient vectors b, a | symbolic or 2 × row vector |
| `lccde/system-response` | y[n] for x = aⁿu[n] or u[n] via partial fractions, incl. pole–zero cancellations | symbolic, numbers |
| `ztransform/z-transform-pair` | X(z) and ROC of a 1–2 term signal (shifted/left-sided exponentials, n aⁿu[n], finite, damped cosine) | symbolic, MC |
| `ztransform/roc-of-signal` | ROC of a 2–3 term signal, including empty ROCs and 0/∞ subtleties | MC |
| `ztransform/z-transform-properties` | transform of x[n−k], n x[n], cⁿx[n], x[−n], x * x, x[n] − x[n−1] from a known pair, with ROC | symbolic, MC |
| `ztransform/inverse-pfe` | partial fractions, ROC from "causal / stable / left-sided", h[n] for n ≥ 0 and n ≤ −1 | numbers, MC, 2 × symbolic |
| `ztransform/all-possible-rocs` | number of ROCs, the stable ROC, h[−2], h[0], h[2] for a given ROC (cancellation traps) | integer, MC, numbers |
| `ztransform/improper-long-division` | improper H(z): long division in z⁻¹, residues, h[0..3] | numbers |
| `ztransform/complex-poles` | conjugate poles r e^{±jθ}: h[0..2] and h[n] = rⁿ(A cos θn + B sin θn) | numbers, symbolic |
| `stability/series-parallel` | overall H(z) of a series/parallel connection; stability of each part and of the whole | symbolic, 3 × Yes/No |
| `stability/stability-from-roc` | causal? stable? for a given ROC; which ROC would be stable | 3 × MC |
| `stability/bibo-from-impulse-response` | is h absolutely summable? Σ\|h[n]\| when it is | Yes/No, number |
| `stability/unbounded-outputs` | marginally stable system: which bounded inputs give unbounded outputs (2/3 rad vs 2π/3 trap, cancellations) | checkbox |
| `lccde/two-sided-recursion` | stable two-sided system as a forward + a backward recursion: which direction for each pole, the backward recursion's coefficients, y[−1], y[0], y[1] for a two-sample input | MC, numbers |
| `stability/parameters-for-stability` | values of a parameter K for BIBO stability (interval endpoints, single value or MC) | numbers or MC |
| `dtft/dtft-quantities` | X_d(0), X_d(π), ∫X_d dω, ∫\|X_d\|² dω without computing the DTFT | numbers |
| `dtft/dtft-of-finite-sequence` | X_d(ω) of a short sequence (exp or cosine form accepted) and \|X_d(ω₀)\| | symbolic, number |
| `dtft/dtft-properties` | DTFT of a shifted, modulated, reversed or self-convolved signal from a known pair (MC of concrete expressions), plus one value the property gives directly | MC, number |
| `dtft/frequency-response` | \|H_d(0)\|, \|H_d(π)\|, the output amplitude for a cosine at π/2, and the filter type (lowpass, highpass, …) | numbers, MC |

Tags let you filter in PrairieLearn: `L1`–`L14`, `hw1`–`hw5`, `midterm1`, the past-exam problem
families `fam-*` (with how many of the seven past Midterm 1 exams contain each), input types and
`easy`/`medium`/`hard`.

## Assessments

- **Practice 1–4** (Homework type, open all semester): each question is worth 6 points — 2 for the
  first correct variant, more for consecutive correct variants (PrairieLearn's mastery scoring) —
  with 3 tries per variant and a new random variant whenever a student wants one.
- **Mock Midterm 1** (Exam type, 2 hours, `multipleInstance`): ten problems in the shape of a past
  Midterm 1 — true/false, property table, finite convolution, infinite convolution or
  superposition, finding h, z-transform with ROC, inverse transform / all ROCs / two-sided recursion, transfer function
  or system response, unbounded outputs, parameters for stability — drawn from alternative pools.
  Three attempts per problem (10/7/4 points); questions and solutions are visible after finishing;
  students can start a fresh instance any time.

Dates are in `courseInstances/Fa26/` (course instance published 2026-08-24 → 2026-12-20).

## Run it in PrairieLearn

Push this repository to the GitHub repo connected to your PrairieLearn course and sync, or run a
local PrairieLearn on it:

```bash
docker run -it --rm -p 3000:3000 -v "$PWD":/course prairielearn/prairielearn
# open http://localhost:3000 → Load from disk → ECE 310
```

## Testing without a PrairieLearn server

`tools/` contains an offline harness that runs each question through PrairieLearn's pipeline
(`generate → prepare → render → parse → grade`, plus the instructor **Test** button) with
PrairieLearn's **real element code**, fetched from GitHub at a pinned commit.

```bash
bash tools/setup_harness.sh               # once, ~15 s: PrairieLearn elements + sympy etc. into tools/.cache/
python3 tools/test_lib.py                 # self-test of serverFilesCourse/ece310 against scipy
python3 tools/test_questions.py           # everything (~15 min; see AUTHORING.md to run topics in parallel)
python3 tools/test_questions.py --only ztransform/inverse-pfe --seeds 25 --gen-seeds 300 -v
python3 tools/preview.py ztransform/inverse-pfe --seed 4 --submit correct --png   # look at a variant
```

For every question the runner checks the course JSON against PrairieLearn's schemas, generates 300
variants and compares them with the question's independent checker (`tools/checks/`), pushes 25
variants through rendering, the Test button (correct / incorrect / invalid) and extra
submissions (equivalent forms must score 1, classic mistakes 0), and compiles every formula.
Last full run: 29 questions, `ALL OK`.

The harness does not depend on this course: copy `tools/plsim.py`, `tools/test_questions.py`,
`tools/checklib.py`, `tools/preview.py` and `tools/setup_harness.sh` into another PrairieLearn
course to test its questions the same way.

## Adding questions

Read [AUTHORING.md](AUTHORING.md): question anatomy, rules for `server.py`/`question.html`,
course notation, the shared library, checker format and the definition of done. In short: write
the question, write `tools/checks/<topic>__<slug>.py`, run the runner until it prints `OK`, look
at two previews.

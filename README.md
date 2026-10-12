# ECE 310 — PrairieLearn course (randomized practice)

60 randomized PrairieLearn "v3" questions for ECE 310 (Digital Signal Processing), covering every
lecture so far (Lectures 1–17, Homework 1–6) with conceptual, symbolic and numerical questions, plus
the Midterm 2 topics that come after Lecture 17 (ideal D/A, digital filtering of analog signals, the
DFT; tag `ahead`), learned from the seven past Midterm 2 exams and the 2019 course notes.
Every variant is a new random problem — new signals, poles, coefficients — with "exam-nice"
numbers, an autograded answer and a worked solution in the answer panel that names the lecture,
homework or past-exam problem it is modelled on. Each question has an independent checker that
recomputes its answers with numpy/scipy/SymPy, and the whole course is tested offline with
PrairieLearn's own element code (see [Testing](#testing-without-a-prairielearn-server)).

```
infoCourse.json                       topics, tags (lectures L1–L17, ahead, hw1–hw6, exam families fam-*/fam2-*), assessment sets
courseInstances/Fa26/
  infoCourseInstance.json             publishing dates
  assessments/
    PR1-signals-systems-convolution/  Practice 1  (L1–L5)
    PR2-z-transform/                  Practice 2  (L6–L8)
    PR3-transfer-functions-stability/ Practice 3  (L9–L11)
    PR4-dtft/                         Practice 4  (L12–L14: template matching, the DTFT)
    PR5-frequency-response/           Practice 5  (L15–L16: frequency response, magnitude and phase)
    PR6-sampling/                     Practice 6  (L17: ideal A/D conversion)
    PR7-ahead-reconstruction-dft/     Practice 7  (after L17: D/A, A/D → H_d → D/A, the DFT)
    ME1-mock-midterm1/                Mock Midterm 1: 10 problems in past-exam shape, 2 h, retakeable
    ME2-mock-midterm2/                Mock Midterm 2: 9 problems in past-exam shape, 2 h, retakeable
questions/<topic>/<slug>/             info.json, server.py, question.html
serverFilesCourse/ece310/             shared library (exact arithmetic, polynomials in z^-1, ROCs, PFE, SVG figures …)
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
| `signals/sampling-and-sequences` | a sequence built from u and δ → its samples; sampling a sinusoid x_a(nT) → ω₀ and one sample; classify a signal (analog / discrete-time / digital …) | integer, row vector, numbers, MC |
| `signals/geometric-sums` | finite, infinite, left-sided, two-sided and roots-of-unity sums; a closed form in n; where a series converges | number, symbolic, MC |
| `signals/complex-exponentials` | Euler split c e^{jω₀n} + c* e^{−jω₀n}; phasor sum of two sinusoids → A cos(ω₀n + θ); \|x[n]\| and x[k] of B aⁿ | numbers, symbolic |
| `systems/lti-from-io-pairs` | what three input–output pairs prove (consistent with LTI / not linear / not time-invariant / cannot be LTI); the LTI output for a new input | MC, integer, row vector |
| `convolution/convolution-facts` | first and last index, length, Σy, end samples and one interior sample of x * h without the full convolution; a property MC | integers, MC |
| `lccde/block-diagram` | a drawn direct-form-I diagram → b, a (term by term), FIR or IIR, h[0], h[1], h[2] | 2 × row vector, MC, numbers |
| `convolution/template-matching` | matched filter h[n] = t[−n], where y peaks and its value, the normalized score (Lecture 12; not tested on exams) | integer, row vector, numbers |
| `ztransform/polynomial-transforms` | finite sequences ↔ polynomials in z and z⁻¹, the ROC at 0 and ∞, h[n] from H(z), an output sample from X(z)H(z), Y(±1) | symbolic, MC, integer, row vector |
| `stability/pole-zero-plot` | a drawn pole–zero plot and a gain condition → H(z); causal ROC and stability; the stable system's ROC and sidedness | symbolic, 4 × MC |
| `convolution/system-interconnection` | a drawn series/parallel network of FIR blocks with a delay → overall h[n], H(z) and one output | integer, row vector, symbolic |
| `dtft/periodicity` | fundamental period of sums/products of sinusoids (or "not periodic"); which sequences equal cos(ω₀n); the alias of ω₁ > π in [0, π] | MC, integer, checkbox, number |
| `dtft/inverse-dtft-ideal` | inverse DTFT of a drawn ideal low-/high-/band-pass or two-level spectrum, possibly delayed: x[n] for n ≠ n₀, x[n₀], one more sample | symbolic (graded numerically), numbers |
| `dtft/symmetry-and-existence` | X_d at negative or shifted frequencies from Hermitian symmetry and periodicity; is h real? from H_d; does the DTFT exist, and in what sense | numbers, MC |
| `dtft/modulation` | Y_d of x[n]cos(ω₀n), x[n]cos²(ω₀n) or (−1)ⁿx[n] for a drawn X_d: values (overlaps, wrap-around) and band edges | numbers |
| `dtft/true-false` | six Midterm 2 true/false statements on Lectures 13–16 from a 34-statement bank (19 parameterized) | 6 × True/False |
| `freqresp/sinusoidal-response` | output for a constant + exponential + sinusoid through H_d given in closed form, as h[n], H(z) or an LCCDE; is h real?; the e^{±jω₂n} coefficients | numbers, MC |
| `freqresp/magnitude-phase-fir` | symmetric or antisymmetric FIR: H_d = e^{−jαω}A(ω), magnitudes, principal phase (π jumps), group delay, first zero | numbers, symbolic |
| `freqresp/phase-and-group-delay` | closed-form H_d (e^{−jMω}sin Kω, (b + cos ω)e^{±j\|ω\|}, first-order IIR): magnitude, principal phase, group delay, where the phase jumps by π | numbers, checkbox |
| `freqresp/pole-zero-and-magnitude` | match four drawn pole–zero plots to low-/high-/band-pass/notch shapes (enrichment: the distance rule); \|H_d(0)\| or \|H_d(π)\| | matching, number |
| `sampling/nyquist-rate` | bandlimit, Nyquist rate and Nyquist period of x_a(t) built from sums, sincs, products, x(ct), modulation, convolution; which f_s avoid aliasing | numbers, checkbox |
| `sampling/sampled-spectrum` | X_d(ω) of a sampled triangle/trapezoid/rectangle/V (drawn): band edge, X_d(0), X_d(π), one more value; guard band, touching or overlapping copies | numbers, MC |
| `sampling/sketch-sampled-spectrum` | the same problem answered by drawing: mark the corners of X_d(ω) on a `pl-sketch` canvas | sketch |
| `sampling/sinusoid-aliasing` | a sampled sinusoid → ω and φ in [0, π] (folding flips the phase), the lowest consistent analog frequency, which frequencies or periods give stated samples | numbers, checkbox |
| `sampling/ideal-reconstruction` | *ahead:* ideal D/A of a sampled cosine; zero-order-hold and sinc-interpolation values; which inputs give a stated output | numbers, checkbox |
| `sampling/ad-filter-da` | *ahead:* A/D → ideal digital filter → D/A: output amplitudes at listed frequencies (aliasing included); largest T and ω_c for a design | numbers |
| `sampling/true-false` | *ahead:* six Midterm 2 true/false statements on sampling, D/A, the DFT, windows and the FFT (37-statement bank, 24 parameterized) | 6 × True/False |
| `dft/dft-values` | *ahead:* X[0], X[N/2], X[k], Parseval and inverse values for N = 4, 6, 8 | numbers, row vector |
| `dft/dft-properties` | *ahead:* y[n] from Y[k] defined by circular shifts, reversals, modulation, products; the formula y = c e^{j2πℓn/N} x[⟨n − n₀⟩]; which list is the DFT of a real signal | row vector (complex), integers, MC |
| `dft/dft-and-dtft` | *ahead:* DFT bins as DTFT samples, the smallest N containing two frequencies, zero-padding relations X[a] = Y[b] or Y*[b] | integer, numbers, checkbox, MC |
| `dft/spectral-analysis` | *ahead:* bins ↔ hertz, bin spacing, aliases; amplitudes and frequencies of two cosines read off a \|X[k]\| plot | numbers |
| `dft/circular-convolution` | *ahead:* N-point circular vs linear convolution, the smallest N, which samples differ, zero padding for an FFT | row vector, integers, checkbox |

Tags let you filter in PrairieLearn: `L1`–`L17`, `ahead` (taught after Lecture 17), `hw1`–`hw6`,
`midterm1`, `midterm2`, the past-exam problem families `fam-*` (Midterm 1) and `fam2-*` (Midterm 2),
each with how many of the seven past exams contain it, input types (`numeric integer symbolic MC
checkbox matrix matching sketch`) and `easy`/`medium`/`hard`.

## Assessments

- **Practice 1–7** (Homework type, open all semester): each question is worth 6 points — 2 for the
  first correct variant, more for consecutive correct variants (PrairieLearn's mastery scoring) —
  with 3 tries per variant and a new random variant whenever a student wants one.
- **Mock Midterm 1** (Exam type, 2 hours, `multipleInstance`): ten problems in the shape of a past
  Midterm 1 — true/false, property table, finite convolution, infinite convolution or
  superposition, finding h, z-transform with ROC, inverse transform / all ROCs / two-sided recursion, transfer function
  or system response, unbounded outputs, parameters for stability — drawn from alternative pools.
  Three attempts per problem (10/7/4 points); questions and solutions are visible after finishing;
  students can start a fresh instance any time.
- **Mock Midterm 2** (Exam type, 2 hours, `multipleInstance`): nine problems in the shape of a past
  Midterm 2 — true/false; a DTFT or inverse DTFT; magnitude and phase; response to sinusoids; the
  sampled spectrum (numbers or a sketch); a sampling puzzle; A/D → H_d → D/A; DFT properties; DFT
  analysis — drawn from alternative pools, scored like Mock Midterm 1. The last three are `ahead`.

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
bash tools/setup_harness.sh               # once, ~20 s: PrairieLearn elements, sympy, sketchresponse … into tools/.cache/
python3 tools/test_lib.py                 # self-test of serverFilesCourse/ece310 against scipy
python3 tools/test_questions.py           # everything (~45 min on one core; see AUTHORING.md to run topics in parallel)
python3 tools/test_questions.py --only ztransform/inverse-pfe --seeds 25 --gen-seeds 300 -v
python3 tools/preview.py ztransform/inverse-pfe --seed 4 --submit correct --png   # look at a variant
```

For every question the runner checks the course JSON against PrairieLearn's schemas, generates 300
variants and compares them with the question's independent checker (`tools/checks/`), pushes 25
variants through rendering, the Test button (correct / incorrect / invalid) and extra
submissions (equivalent forms must score 1, classic mistakes 0), compiles every formula and flags
typesetting slips ("display lint": `e^{j-…}`, double signs, a written-out `1\cos`). `pl-sketch`
questions are graded by the real element with the `sketchresponse` package PrairieLearn uses, and
checkers submit drawings exactly as a browser sends them (`checklib.sketch_raw`).
Last full run (Oct 12, 2026): 60 questions, `ALL OK`.

The harness does not depend on this course: copy `tools/plsim.py`, `tools/test_questions.py`,
`tools/checklib.py`, `tools/preview.py` and `tools/setup_harness.sh` into another PrairieLearn
course to test its questions the same way.

## Adding questions

Read [AUTHORING.md](AUTHORING.md): question anatomy, rules for `server.py`/`question.html`,
course notation, the shared library, checker format and the definition of done. In short: write
the question, write `tools/checks/<topic>__<slug>.py`, run the runner until it prints `OK`, look
at two previews.

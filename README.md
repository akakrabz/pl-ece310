# ECE 329 — PrairieLearn course (sample questions)

A minimal PrairieLearn course directory: thirteen randomized "v3" questions — five modelled on FA26 Exam 1
(and the Lecture 3 flux challenge), four on HW4 / Lecture 10 material (layered dielectrics, lossy capacitors)
and four on magnetostatics (Lectures 12–15) — one course instance, and three practice assessments.
Formats were checked against docs.prairielearn.com (September 2026).

```
prairielearn/
├── infoCourse.json                      name, title, topics, tags, the "Practice" assessment set
├── courseInstances/Fa26/
│   ├── infoCourseInstance.json          uuid, longName, publishing dates
│   └── assessments/                     Homework-type, autoPoints/maxAutoPoints
│       ├── exam1-practice/              Practice 1: the five Exam 1 families
│       ├── dielectrics-practice/        Practice 2: HW4-style (Lectures 8-10)
│       └── magnetostatics-practice/     Practice 3: Lectures 12-15
├── questions/electrostatics/
│   ├── conservative-field-analysis/     Exam 1 #1 family: curl (MC), V (symbolic), rho (symbolic), line integral (number)
│   ├── nonuniform-slab/                 Exam 1 #2 family: power-law slab, Gauss's law + potential (numeric)
│   ├── dielectric-interface/            Exam 1 #3 family: boundary conditions, D, P, bound charge (numeric)
│   ├── coax-two-layer/                  Exam 1 #4 family: D, E, outer charge, V(a)-V(b), C' (numeric)
│   ├── flux-through-plane/              Lecture 3 challenge: half-space flux argument (MC)
│   ├── two-layer-parallel-plates/       HW4 #1/#3 family: E2, V(d), Vp, plate charge, C/A (numeric)
│   ├── sheet-between-grounded-plates/   HW4 #2 family: V0, E above/below, induced plate charges (numeric + MC)
│   ├── dielectric-insertion/            HW4 #5 family: fixed-Q vs fixed-V fill, energy ratio, conducting fluid (numeric + MC)
│   └── lossy-capacitor-relaxation/      HW4 #4 family: C, G = (sigma/eps)C, tau = eps/sigma, Q(t), I(0); random geometry
├── questions/magnetostatics/
│   ├── ampere-coax-fields/              H_phi in the four regions of a coax (numeric + MC)
│   ├── current-sheets-superposition/    two current sheets: H by region, force per area, direction (numeric + MC)
│   ├── solenoid-inductance-energy/      H, flux per turn, L, W, L/R (numeric)
│   └── faraday-loop-emf/                flux, signed emf and current, Lenz direction (numeric + MC)
└── tools/test_questions.py              offline checker (see below)
```

Each question is `info.json` + `question.html` + `server.py`. All randomization happens in `generate(data)`;
symbolic correct answers are given as sympy-parsable strings (accepted by `pl-symbolic-input`), so the
questions have no dependencies beyond the PrairieLearn runtime.

## Install

Copy the contents of this folder into your PrairieLearn course repository (or point a local PrairieLearn at it):

```bash
# local docker PrairieLearn, from https://docs.prairielearn.com/installing/
docker run -it --rm -p 3000:3000 -v "$PWD":/course prairielearn/prairielearn
```

Then open http://localhost:3000, load the course, and open each question's preview. Sync errors to expect if you
merge into an existing course: duplicate UUIDs (regenerate with `python3 -c "import uuid;print(uuid.uuid4())"`)
and undeclared topics/tags (add them to your `infoCourse.json` or change the questions' `topic`/`tags`).

## Test without PrairieLearn

```bash
python3 tools/test_questions.py 30     # ~1 min; the slab check is the slow one (numerical quadrature per seed)
```

runs every `generate()` on N seeds (default 300) and checks: JSON-serializability, that every `{{params.x}}` used by the
template is set, that every input has a correct answer, that multiple-choice questions have exactly one correct
and distinct choices, and — most importantly — that the closed-form answers agree with independent numerical
computation (finite-difference gradients/divergence for the symbolic field question, quadrature for the slab
and coax; midpoint-rule integrals for the layered and coax fields; shell resistances for G). The delivered version passes all checks; the harness is the place to add a test when you add a question.

## Notes on the design

- **Problem 1 family** generates the field *from* a potential, so it is conservative by construction; a field with
  random polynomial components almost never is, and then "find V" has no answer.
- **Dielectric interface** picks $E_{1z}$ so that $E_{2z}$ is an integer (via the gcd of the two $\epsilon_r$).
- **Flux MC** builds its distractors from the actual mistakes (both charges counted +, closed-surface answer,
  flipped normal, swapped above/below) and de-duplicates them against the correct value.
- Points: Homework-type with `autoPoints`/`maxAutoPoints` (e.g. 5/15) so a student can bank three correct variants;
  `triesPerVariant` 2–3 so the hidden hints (`show-after-submission`) have a chance to appear.
- Hints use `pl-hidden-hints` (the only place `pl-hint` is allowed); `order="fixed"` replaces the deprecated
  `fixed-order`; `accessControl`/`publishing` are the current (non-deprecated) access mechanisms.

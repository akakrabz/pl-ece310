---
title: "ECE 329 · Fields and Waves"
description: "Teaching notes for ECE 329 (Fields and Waves I): lectures rewritten as narratives, a concept glossary that links everything, worked exam-style problems, and interactive demos."
---

Notes for **ECE 329 — Fields and Waves I** (University of Illinois), written to *teach* the material rather than to summarize it. The source material is the course itself — Prof. Kudeki's lecture notes, the lecture slides with their in-class annotations, the homework and exams — digested and rewritten in one consistent notation, with the connections made explicit.

> [!tip] How this site is organized — three layers, one graph
> - **Lectures** (the spine, in course order) tell the story: motivation → definition → derivation → worked example → what goes wrong → what the exam does with it.
> - **Concepts** are the glossary: one short page per idea (Gauss's law, curl, flux, …) with the key equation, when it applies, the traps, and links to every lecture and problem that uses it.
> - **Problems** and **demos** are where the ideas get exercised — exam-style problems worked in full (with the numbers changed), and interactive pages you can drag things around in.
>
> Every page links to its neighbours; hover a link for a preview, and use the **graph view** at the top right of any page to see what connects to what. Fields and Waves is a web of ideas, not a list — the site is built the same way.

## Course map

### [[0-toolkit/index|Toolkit]] — the mathematics assumed
[[0-toolkit/01-coordinates-and-differential-elements|Coordinates and differential elements]] · [[0-toolkit/02-vector-calculus-cheatsheet|Vector calculus cheat sheet]] · [[0-toolkit/03-units-and-constants|Units and constants]] · [[0-toolkit/04-errata-in-the-course-materials|Errata in the course materials]] · [[0-toolkit/04-errata-in-the-course-materials|Errata in the course materials]]

### [[1-electrostatics/index|Unit 1 · Electrostatics]] — Lectures 1–11 (Exam 1: L1–10)
1. [[1-electrostatics/01-fields-forces-and-the-maxwell-roadmap|Fields, forces, and the Maxwell roadmap]]
2. [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Coulomb's law, superposition, and the birth of Gauss's law]]
3. [[1-electrostatics/03-gauss-law-at-work|Gauss's law at work: flux, symmetry, and charge densities]]
4. [[1-electrostatics/04-divergence-and-curl|Divergence and curl: Maxwell's equations at a point]]
5. [[1-electrostatics/05-curl-free-fields-and-the-electrostatic-potential|Curl-free fields and the electrostatic potential]]
6. [[1-electrostatics/06-circulation-and-boundary-conditions|Circulation and boundary conditions]]
7. [[1-electrostatics/07-poisson-and-laplace|Poisson's and Laplace's equations]]
8. [[1-electrostatics/08-conductors-dielectrics-and-polarization|Conductors, dielectrics, and polarization]]
9. [[1-electrostatics/09-static-fields-in-dielectric-media|Static fields in dielectric media]]
10. [[1-electrostatics/10-capacitance-and-conductance|Capacitance and conductance]]
11. Lorentz–Drude models for conductivity and susceptibility — *coming (not on Exam 1)*

### [[2-magnetostatics/index|Unit 2 · Magnetostatics]] — Lectures 12–15
Magnetic force and Ampère's law · current sheets, solenoids, vector potential · Faraday's law · inductance — *planned*

### [[3-maxwell-and-waves/index|Unit 3 · Maxwell's equations and waves]] — Lectures 16–26
Displacement current and the full Maxwell equations · plane TEM waves · Poynting · phasors · lossy media · polarization · reflection and standing waves — *planned*

### [[4-transmission-lines/index|Unit 4 · Transmission lines]] — Lectures 27–38
Guided TEM waves · bounce diagrams · sinusoidal steady state and input impedance · quarter-wave transformers · the Smith chart · impedance matching — *planned*

### Cross-cutting
[[concepts/index|Concept glossary]] · [[problems/index|Worked problems]] · [[demos/index|Interactive demos]]

## Conventions used throughout

| | this site writes | you may also see |
|---|---|---|
| vectors | bold: $\mathbf{E}$, $\mathbf{r}$ | $\vec{E}$ on the slides |
| unit vectors | $\hat{x},\hat{y},\hat{z}$, $\hat{r},\hat{\phi},\hat{\theta}$ | $\hat{a}_x$ on the slides |
| permittivity of free space | $\epsilon_0$ | $\epsilon_o$ in the course notes |
| charge densities | $\rho$ [C/m³], $\rho_s$ [C/m²], $\rho_l$ [C/m] | $\lambda$ for line charge |
| electric flux density | $\mathbf{D} = \epsilon_0\mathbf{E}$ (vacuum), $\epsilon_0\mathbf{E}+\mathbf{P}$ (matter) | "displacement flux density" |
| differential elements | $d\mathbf{l}$, $d\mathbf{S}$ (outward on closed surfaces), $dV$ | $d\mathbf{L}$, $d\mathbf{s}$ |
| Gauss's law | $\oint_S\mathbf{D}\cdot d\mathbf{S} = Q_{\text{enc}}$ | $\nabla\cdot\mathbf{E} = \rho/\epsilon_0$ (same thing in vacuum) |
| the constant | $\dfrac{1}{4\pi\epsilon_0}$, never "$k$" | $9\times10^9$ |

This follows the Midterm 1 formula sheet (the exam papers themselves write $\lambda$ for line charge and overbar vectors). Units go in square brackets after every result.

## Callouts

The notes use a few recurring boxes, so you can skim for what you need:

> [!key] Key
> The result to remember — a law, a boxed formula, a final answer.

> [!recipe] Recipe
> A procedure: how to actually do the calculation, step by step.

> [!trap] Trap
> A mistake that costs points, and how to avoid it. Several are taken from real exam keys.

> [!exam] Exam
> How the idea shows up on exams, with the FA26 Exam 1 mapping.

> [!intuition] Intuition
> The picture or analogy behind the mathematics.

> [!derivation]- Derivation (click to expand)
> Longer derivations are folded so the narrative stays readable. Expand when you want the details.

*Status: Unit 1 through the Exam 1 scope (Toolkit, Lectures 1–10, 25 concept pages, one worked problem per Exam 1 problem) is written; Lecture 11 and the remaining units are outlined. Sources: E. Kudeki, ECE 329 Lecture Notes (2026); lecture slides (Shao, adapted from Goddard and Cunningham); FA26 homework and Exam 1; N. N. Rao, Fundamentals of Electromagnetics for Electrical and Computer Engineering.*

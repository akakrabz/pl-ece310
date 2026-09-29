---
title: "The 5-step recipe for fields of charge distributions"
description: "Draw, chop, write dE of one piece, kill components by symmetry, integrate — the procedure for every superposition integral, and when to abandon it for Gauss's law."
tags: [concept, electrostatics, exam-1]
aliases: ["5-step program", "superposition integral recipe"]
---

> [!recipe] The patented 5-step program
> 1. **Draw** — large and clear; a cross-section if the problem is 3-D; a coordinate system that matches the symmetry.
> 2. **Chop** the distribution into pieces $dQ$: $\rho_l\,dl$, $\rho_s\,dS$, $\rho\,dV$, written in your coordinates (e.g. $dQ = \rho_s\,r\,dr\,d\phi$ on a disk).
> 3. **One piece**: $d\mathbf{E} = \dfrac{dQ}{4\pi\epsilon_0R^2}\hat{R}$, with $R$ and $\hat{R}$ (source → field point) expressed in the coordinates.
> 4. **Symmetry**: pair each piece with its mirror image; drop the components that cancel; keep the surviving component along a direction that is *fixed at the field point*.
> 5. **Integrate** the scalar that is left.

**Why step 4 matters.** It converts a vector integral with a direction that changes from piece to piece into a scalar integral times a fixed unit vector. Without it you must resolve $\hat{R}$ into Cartesian components and integrate each — legal, but three times the work. Never pull a cylindrical or spherical unit vector through an integral unless it refers to the field point ([[0-toolkit/01-coordinates-and-differential-elements]]).

**Worked instances.**
- Finite line charge → $E_r = \dfrac{\rho_l}{4\pi\epsilon_0r}\dfrac{2a}{\sqrt{a^2+r^2}}$; infinite limit $\dfrac{\rho_l}{2\pi\epsilon_0r}$ ([[1-electrostatics/02-coulombs-law-superposition-and-gauss#4-continuous-distributions--the-5-step-program|Lecture 2 §4]]).
- Charged rod's force on a point charge beyond its end → $\dfrac{qQ}{4\pi\epsilon_0\,a(L+a)}\hat{x}$ (Cartesian $\hat{x}$ comes out of the integral).
- Two equal and opposite charges, field point on their perpendicular bisector → double one contribution's surviving component (Lecture 2 §3). On the dipole *axis* the contributions are collinear and must be subtracted instead.

**When to skip it.** If the distribution is spherically, cylindrically or planarly symmetric *and infinite/complete* (whole sphere, infinite line, infinite sheet), [[concepts/gauss-law|Gauss's law]] gives the same answer in four lines. Use the recipe when the symmetry is partial — a finite rod, a ring, a disk, a dipole — or when you need the field somewhere the Gaussian trick doesn't reach.

**Where it appears.** [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]]; homework 1–2; [[1-electrostatics/05-curl-free-fields-and-the-electrostatic-potential|Lecture 5]] repeats it for the potential (a scalar integral from the start — even easier).

Related: [[concepts/superposition]] · [[concepts/coulombs-law]].

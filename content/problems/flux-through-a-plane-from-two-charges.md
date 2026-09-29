---
title: "Worked problem — Flux through a plane from two point charges"
description: "The half-space argument: a point charge sends exactly half its flux through any plane that doesn't contain it. Superpose, watch the sign of the normal, and the Lecture 3 challenge question is a one-line count."
tags: [problem, electrostatics, exam-1]
---

*Problem · style of the Lecture 3 in-class challenge · uses [[1-electrostatics/03-gauss-law-at-work|Lecture 3]] · concepts: [[concepts/flux]], [[concepts/gauss-law]], [[concepts/superposition]]*

> [!question] Problem
> A charge $+3Q$ sits at $(0,0,+h)$ and a charge $-Q$ at $(0,0,-h)$. Find the electric flux $\psi_E = \int\mathbf{D}\cdot d\mathbf{S}$ through the entire $z=0$ plane, taking the plane's normal to be $\hat{n} = +\hat{z}$.
>
> Then: what changes if $\hat{n} = -\hat{z}$? And what is the flux through a closed sphere of radius $2h$ centred at the origin?

## The one idea: half of a point charge's flux crosses any plane

The flux of $\mathbf{D}$ out of a *closed* surface around $Q$ is $Q$ (Gauss). A point charge radiates its field lines uniformly in all directions, so an infinite plane that does not pass through the charge is crossed by exactly **half** of them: $|\psi| = Q/2$. (Rigorous version: the plane subtends a solid angle $2\pi$ out of $4\pi$ as seen from the charge; or, close the plane with a huge hemisphere on the charge's side — the hemisphere catches $Q/2$ by symmetry, and the closed combination must catch $Q$.)

The **sign** is fixed by whether the lines cross the plane *along* or *against* the chosen $\hat{n}$.

## Count the contributions

Normal $\hat{n} = +\hat{z}$.

- **$+3Q$ above the plane.** Its lines point *away* from it; the half that reach the plane cross it going **downward**, against $\hat{n}$: contribution $-\dfrac{3Q}{2}$.
- **$-Q$ below the plane.** Its lines point *toward* it; the ones crossing the plane come from above and go **downward** into it, again against $\hat{n}$: contribution $-\dfrac{Q}{2}$.

$$
\psi_E = -\frac{3Q}{2} - \frac{Q}{2} = -2Q .
$$

> [!key] Answer
> $\psi_E = -2Q$ [C] through the $z=0$ plane with $\hat{n}=+\hat{z}$.

The two charges' flux lines *both* cross the plane downward — a positive charge above pushes lines down through it, a negative charge below pulls lines down through it — which is why the magnitudes add rather than cancel.

## The follow-ups

**Flip the normal** to $\hat{n} = -\hat{z}$: every sign flips, $\psi_E = +2Q$. The physics did not change; the bookkeeping convention did. This is why the problem statement must specify $\hat{n}$ for an open surface, and why a closed surface (outward normal, no choice) is the natural home of Gauss's law.

**Closed sphere of radius $2h$**: it encloses both charges, so $\oint\mathbf{D}\cdot d\mathbf{S} = Q_{\text{enc}} = 3Q - Q = 2Q$. No half-space argument needed — only what is *inside* counts, and the geometry of the surface is irrelevant.

**A sphere of radius $h/2$ centred at the origin** encloses neither charge: flux zero, even though $\mathbf{D}\ne0$ on its surface (what enters leaves).

> [!trap] Three ways to get this wrong
> - Assigning $+3Q/2$ to the positive charge "because it is positive". Sign comes from the crossing direction relative to $\hat{n}$, not from the sign of the charge alone.
> - Adding the charges first ($3Q - Q = 2Q$) and then halving gives $\pm Q$, which is wrong: through the plane each charge's half-flux carries its *own* sign. (The unhalved $2Q$ is the flux through a *closed* surface around both charges.)
> - Forgetting that lines from $-Q$ arrive from *both* sides in general — but the ones that cross the plane must come from the far side, downward.

## Variants for practice

| charges | $\hat{n}$ | $\psi$ through $z=0$ |
|---|---|---|
| $+Q$ at $+h$, $-Q$ at $-h$ (dipole) | $+\hat{z}$ | $-Q/2 - Q/2 = -Q$ |
| $+2Q$ at $+h$, $-Q$ at $-h$ | $+\hat{z}$ | $-Q - Q/2 = -3Q/2$ |
| $+Q$ at $+h$, $+Q$ at $-h$ | $+\hat{z}$ | $-Q/2 + Q/2 = 0$ |
| $+Q$ at $+h$ only | $-\hat{z}$ | $+Q/2$ |

Related demo: drag charges and a *closed* Gaussian loop in [[demos/point-charges-and-gauss]] and watch the flux track the enclosed charge — the closed-surface half of this problem.

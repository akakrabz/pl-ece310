---
title: "Flux"
description: "The surface integral of the normal component of a field — 'how many field lines pierce the surface' — and the two cases where it is trivial to evaluate."
tags: [concept, electrostatics, exam-1]
aliases: ["electric flux", "surface integral"]
---

> [!key] Definition
> $$
> \psi = \int_S \mathbf{F}\cdot d\mathbf{S} = \int_S F_n\,dS,\qquad d\mathbf{S} = \hat{n}\,dS
> $$
> Through a small flat patch: $\Delta\psi = \mathbf{F}\cdot\Delta\mathbf{S} = F\cos\alpha\,\Delta S$, with $\alpha$ the angle between $\mathbf{F}$ and the normal. Electric flux $\psi_E = \oint\mathbf{D}\cdot d\mathbf{S}$ is in coulombs; magnetic flux $\psi_B = \int\mathbf{B}\cdot d\mathbf{S}$ is in webers.

**Picture.** Count the field lines that pass through the surface, with sign: along $\hat{n}$ counts $+$, against counts $-$. Flux depends on the field strength, the area, and the tilt.

> [!tip] The two cases that make Gauss's law usable
> 1. $\mathbf{F}$ uniform in magnitude and parallel to $d\mathbf{S}$ over a piece of surface → that piece contributes $F\cdot A$.
> 2. $\mathbf{F}$ parallel to the surface (perpendicular to $d\mathbf{S}$) → that piece contributes $0$.
>
> Every Gaussian surface in the course (sphere, coaxial cylinder, pillbox) is assembled from pieces of these two kinds.

**Orientation.** Closed surface: $d\mathbf{S}$ outward, no choice. Open surface: *you* choose $\hat{n}$; the sign of the flux follows your choice. A point charge's flux through an infinite plane that doesn't contain it is $\pm Q/2$ ([[1-electrostatics/03-gauss-law-at-work#2-gausss-law-restated-for-use|half-space argument]]).

**Superposition.** Flux is linear in the field, so fluxes of separate sources add.

**Two theorems about it.** [[concepts/gauss-law|Gauss's law]] (flux of $\mathbf{D}$ out of a closed surface = enclosed charge) is physics; the [[concepts/divergence-theorem|divergence theorem]] (flux out = volume integral of the divergence) is mathematics.

**Where it appears.** [[1-electrostatics/03-gauss-law-at-work#1-flux-counting-arrows-through-a-surface|Lecture 3 §1]] (definition and intuition), Lecture 14 (Faraday's law: EMF = −dψ_B/dt), [[problems/flux-through-a-plane-from-two-charges]].

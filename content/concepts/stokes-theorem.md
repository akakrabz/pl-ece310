---
title: "Stokes' theorem"
description: "Circulation around a closed curve equals the flux of the curl through any surface it bounds — the bridge between the integral and differential forms of Faraday's and Ampère's laws."
tags: [concept, toolkit, exam-1]
aliases: ["Stokes theorem", "curl theorem"]
---

> [!key] Statement
> $$
> \oint_C\mathbf{A}\cdot d\mathbf{l} = \int_S(\nabla\times\mathbf{A})\cdot d\mathbf{S},\qquad C = \partial S,\ \text{right-hand rule between } d\mathbf{l} \text{ and } d\mathbf{S}
> $$
> Any surface bounded by $C$ gives the same right-hand side.

**Why it is true.** Tile $S$ with tiny loops. Each interior edge is traversed in opposite directions by its two neighbouring loops and cancels; only the outer boundary survives. Each tiny loop's circulation is (curl)·(area) — that is the *definition* of curl as circulation per unit area.

**What it is for.**
- Converts the curl equations between forms: $\oint_C\mathbf{E}\cdot d\mathbf{l} = -\frac{\partial}{\partial t}\int_S\mathbf{B}\cdot d\mathbf{S}$ for every $S$ ⟺ $\nabla\times\mathbf{E} = -\partial\mathbf{B}/\partial t$; likewise Ampère.
- Proves the equivalences for [[concepts/conservative-field|conservative fields]]: $\nabla\times\mathbf{E} = 0$ everywhere ⟹ $\oint_C\mathbf{E}\cdot d\mathbf{l} = 0$ for every loop ⟹ path independence ⟹ a potential.
- Interprets an EMF: $\oint_C\mathbf{E}\cdot d\mathbf{l}$ around a loop equals the total curl threading it.

> [!example] A curl trap, checked with Stokes
> $\mathbf{E} = E_0(2xy\hat{x}+\kappa x^2\hat{y})$ has $\nabla\times\mathbf{E} = 2(\kappa-1)E_0\,x\,\hat{z}$. Around the unit square in the $xy$-plane (counter-clockwise), Stokes predicts $\oint\mathbf{E}\cdot d\mathbf{l} = \int_0^1\!\!\int_0^1 2(\kappa-1)E_0x\,dx\,dy = (\kappa-1)E_0$; direct integration of the four legs gives the same. For $\kappa = 1$ the field is conservative and the loop integral vanishes.

**Where it appears.** [[1-electrostatics/04-divergence-and-curl#3-curl-circulation-per-unit-area|Lecture 4 §3]], [[1-electrostatics/06-circulation-and-boundary-conditions|Lecture 6]] (circulation, boundary conditions), Lecture 14 (Faraday). Its divergence-side twin is the [[concepts/divergence-theorem]].

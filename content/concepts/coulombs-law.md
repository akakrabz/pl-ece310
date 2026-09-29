---
title: "Coulomb's law"
description: "The 1/r² radial field of a stationary point charge, with the source-to-field-point convention for the unit vector."
tags: [concept, electrostatics, exam-1]
aliases: ["Coulomb field", "Coulomb's law"]
---

> [!key] Statement
> Field of a stationary charge $Q$ at $\mathbf{r}'$, observed at $\mathbf{r}$:
> $$
> \mathbf{E}(\mathbf{r}) = \frac{Q}{4\pi\epsilon_0 R^2}\,\hat{R},\qquad
> \mathbf{R} = \mathbf{r}-\mathbf{r}',\ R = |\mathbf{R}|,\ \hat{R} = \mathbf{R}/R .
> $$
> Force on a second charge $q$ there: $\mathbf{F} = q\mathbf{E}$. Like charges repel (force along $+\hat{R}$), unlike attract.

**The convention that prevents sign errors.** Number the *source* 1 and the *field point* 2. $\hat{R}$ points **from the source to the field point** ("final minus initial"). The field of a positive charge points along $+\hat{R}$ (away from it); of a negative charge along $-\hat{R}$ (toward it). Newton's third law: $\mathbf{F}_1 = -\mathbf{F}_2$.

**Off-origin form** (used in every superposition integral):
$$
\mathbf{E}(\mathbf{r}) = \frac{Q}{4\pi\epsilon_0}\frac{\mathbf{r}-\mathbf{r}'}{|\mathbf{r}-\mathbf{r}'|^3}.
$$
The cube in the denominator is not a typo: one power normalizes the direction vector.

**Scope.** Exact for charges at rest; for slowly moving charges it is the leading (quasi-static) term. For continuous distributions, replace $Q$ by $dQ = \rho_l\,dl,\ \rho_s\,dS,\ \rho\,dV$ and integrate ([[concepts/five-step-recipe]]). [[concepts/gauss-law]] is derived from it and outlives it: Gauss holds even when Coulomb fails (moving charges, time-varying fields).

**Numbers.** $\dfrac{1}{4\pi\epsilon_0}\approx 9\times10^{9}$ m/F — but the course keeps $4\pi\epsilon_0$ explicit because $\epsilon_0\to\epsilon$ inside dielectrics (Lectures 8–9).

**Where it appears.** [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]] (statement, superposition, line charge), [[1-electrostatics/04-divergence-and-curl#6-static-electric-fields-are-curl-free|Lecture 4]] (its curl is zero).

---
title: "Superposition"
description: "Fields of separate sources add as vectors — the linearity that turns Coulomb's law into integrals and lets you assemble complicated distributions from known building blocks."
tags: [concept, electrostatics, exam-1]
---

> [!key] Principle
> Maxwell's equations are linear in the fields and sources. Therefore the field of several sources is the **vector sum** of the fields each would produce alone:
> $$
> \mathbf{E} = \sum_n\mathbf{E}_n = \sum_n\frac{Q_n}{4\pi\epsilon_0 R_n^2}\hat{R}_n
> \quad\longrightarrow\quad
> \mathbf{E} = \frac{1}{4\pi\epsilon_0}\int\frac{dQ'}{R^2}\hat{R}.
> $$
> Fluxes, potentials and charges add too.

**Two ways to use it.**
1. **Integrate.** Chop a distribution into $dQ$, write the Coulomb field of one piece, use symmetry, integrate — the [[concepts/five-step-recipe]].
2. **Assemble.** Recognize a distribution as a sum of things whose fields you already know (sheet, slab, line, point, cylinder), shift each to its own location, and add. This is how [[1-electrostatics/03-gauss-law-at-work#4-superposition-of-the-building-blocks|Lecture 3]] gets the parallel-plate field and the pn-junction field with no new integrals.

**Superposition of fluxes.** In [[concepts/gauss-law|Gauss's law]] $Q_{\text{enc}}$ is the algebraic sum of enclosed charges; for an *open* surface each charge contributes its own $\pm$ share (the half-space argument in [[problems/flux-through-a-plane-from-two-charges]]).

**Symmetry is superposition in disguise.** "By symmetry, the $z$-components cancel" means: pair each element with its mirror image and add their two contributions first (the dipole on its bisector in Lecture 2 §3 is the two-charge version).

> [!trap] Superpose *vectors*, not magnitudes
> $|\mathbf{E}_1+\mathbf{E}_2|\ne|\mathbf{E}_1|+|\mathbf{E}_2|$ unless the two point the same way. Resolve into components (Cartesian, or a fixed direction at the field point) before adding.

**Where it appears.** [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]] (dipole on its bisector, continuous distributions), [[1-electrostatics/03-gauss-law-at-work|Lecture 3]] (sheets and slabs), [[1-electrostatics/05-curl-free-fields-and-the-electrostatic-potential|Lecture 5]] (potentials add as scalars — even easier).

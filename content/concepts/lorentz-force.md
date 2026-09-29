---
title: "Lorentz force"
description: "F = q(E + v×B): the definition of the electric and magnetic fields through the force on a charge."
tags: [concept, electrostatics, magnetostatics]
aliases: ["Lorentz force law"]
---

> [!key] Definition
> $$
> \mathbf{F} = q\left(\mathbf{E}+\mathbf{v}\times\mathbf{B}\right)
> $$
> $\mathbf{E}$ = force per unit *stationary* charge (V/m). $\mathbf{B}$ is defined by the velocity-dependent part of the force, $\mathbf{F}_M = q\mathbf{v}\times\mathbf{B}$: the magnetic force per unit charge is $\mathbf{v}\times\mathbf{B}$, always perpendicular to both $\mathbf{v}$ and $\mathbf{B}$ (units of $\mathbf{B}$: T = Wb/m² = N·s/(C·m)).

**What it is for.** It is not derived from anything — it is how $\mathbf{E}$ and $\mathbf{B}$ are *defined*, operationally, by what they do to a test charge. Newton's law $m\,d\mathbf{v}/dt = \mathbf{F}$ then predicts the motion (valid for $v\ll c$).

**Consequences worth knowing.**
- The magnetic force does no work ($\mathbf{F}_M\perp\mathbf{v}$); it changes direction, not speed → for $\mathbf{v}\perp\mathbf{B}$, circular motion of radius $R = mv/(|q|B)$ in a uniform $\mathbf{B}$ (a helix in general).
- Crossed fields select a velocity: $\mathbf{F}=0$ when $v = E/B$.
- For a current element, $q\mathbf{v}\to I\,d\mathbf{l}$, so $d\mathbf{F} = I\,d\mathbf{l}\times\mathbf{B}$ (motors, Lecture 12).
- $\mathbf{E}$ and $\mathbf{B}$ mix between reference frames: $\mathbf{E}'\approx\mathbf{E}+\mathbf{v}\times\mathbf{B}$. A charge at rest for you is a current for me — which is why Maxwell's equations couple the two fields.

**Where it appears.** [[1-electrostatics/01-fields-forces-and-the-maxwell-roadmap|Lecture 1]] (definition), [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]] (velocity and mass selectors), Lecture 12 (magnetic force on currents), Lecture 14 (motional EMF).

Related: [[concepts/electric-field]] · [[concepts/coulombs-law]].

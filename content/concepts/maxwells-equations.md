---
title: "Maxwell's equations"
description: "The four equations in integral and differential form, what each says physically, the static special case, and where in the course each one is derived and used."
tags: [concept, electrostatics, magnetostatics, waves]
aliases: ["Maxwell equations", "field equations"]
---

> [!key] Free-space form (course notation: $\mathbf{D} = \epsilon_0\mathbf{E}$, $\mathbf{B} = \mu_0\mathbf{H}$)
> | law | integral form | differential form | says |
> |---|---|---|---|
> | Gauss | $\oint_S\mathbf{D}\cdot d\mathbf{S} = \int_V\rho\,dV$ | $\nabla\cdot\mathbf{D} = \rho$ | E-lines start/end on charge |
> | Gauss (magnetic) | $\oint_S\mathbf{B}\cdot d\mathbf{S} = 0$ | $\nabla\cdot\mathbf{B} = 0$ | B-lines never start or end |
> | Faraday | $\oint_C\mathbf{E}\cdot d\mathbf{l} = -\dfrac{\partial}{\partial t}\int_S\mathbf{B}\cdot d\mathbf{S}$ | $\nabla\times\mathbf{E} = -\dfrac{\partial\mathbf{B}}{\partial t}$ | changing B makes E circulate |
> | Ampère–Maxwell | $\oint_C\mathbf{H}\cdot d\mathbf{l} = \int_S\Big(\mathbf{J}+\dfrac{\partial\mathbf{D}}{\partial t}\Big)\cdot d\mathbf{S}$ | $\nabla\times\mathbf{H} = \mathbf{J}+\dfrac{\partial\mathbf{D}}{\partial t}$ | current and changing D make H circulate |
>
> Plus the [[concepts/lorentz-force|Lorentz force]] $\mathbf{F} = q(\mathbf{E}+\mathbf{v}\times\mathbf{B})$, which says what the fields *do*.

**How to read them.** Left column: what the fields emanate from (divergence = sources). Right column: what makes them circulate (curl). Sources on the right-hand sides, fields on the left. The integral and differential forms are equivalent through the [[concepts/divergence-theorem]] and [[concepts/stokes-theorem]] plus "true for every surface/volume ⟹ true at every point".

**Built in: charge conservation.** $\nabla\cdot(\nabla\times\mathbf{H}) = 0$ gives $\nabla\cdot\mathbf{J} = -\partial\rho/\partial t$.

**The static split (Lectures 1–13).** With $\partial/\partial t = 0$ the system decouples:
- *Electrostatics*: $\nabla\cdot\mathbf{D} = \rho$, $\nabla\times\mathbf{E} = 0$ — curl-free, so $\mathbf{E} = -\nabla V$ and $\nabla^2V = -\rho/\epsilon_0$.
- *Magnetostatics*: $\nabla\cdot\mathbf{B} = 0$, $\nabla\times\mathbf{H} = \mathbf{J}$ — divergence-free, so $\mathbf{B} = \nabla\times\mathbf{A}$.

By the Helmholtz theorem, divergence + curl (with decay at infinity) determine a field uniquely, so each pair is a complete description.

**The dynamic coupling (Lectures 14–38).** Lecture 14 (Faraday's law) puts $\partial\mathbf{B}/\partial t$ back into the curl-$\mathbf{E}$ equation — the first coupling, and the end of "voltage is path-independent"; Lecture 16 adds the displacement current, and from then on Faraday and Ampère–Maxwell feed each other: a changing $\mathbf{B}$ makes $\mathbf{E}$, a changing $\mathbf{E}$ makes $\mathbf{B}$ (the displacement current $\partial\mathbf{D}/\partial t$ is Maxwell's contribution). The result is the wave equation with speed $c = 1/\sqrt{\mu_0\epsilon_0}$ — light, radio, and signals on transmission lines.

**In matter.** $\mathbf{D} = \epsilon_0\mathbf{E}+\mathbf{P}$, $\mathbf{H} = \mathbf{B}/\mu_0-\mathbf{M}$; $\rho$ and $\mathbf{J}$ then mean *free* charge and current (Lectures 8–9, 17).

**Where each is derived in the course.** Gauss: [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]] (from Coulomb) and [[1-electrostatics/04-divergence-and-curl|Lecture 4]] (differential). Magnetic Gauss: [[1-electrostatics/03-gauss-law-at-work#6-the-magnetic-counterpart--bds--0|Lecture 3]]. Faraday: Lecture 14. Ampère: Lectures 12–13, completed with displacement current in Lecture 16. The roadmap is in [[1-electrostatics/01-fields-forces-and-the-maxwell-roadmap#2-the-roadmap-maxwells-equations-on-day-one|Lecture 1 §2]].

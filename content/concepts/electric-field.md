---
title: "Electric field E"
description: "Force per unit stationary charge; a vector field with units V/m, sourced by charge, drawn as arrows or field lines."
tags: [concept, electrostatics]
aliases: ["E field", "electric field"]
---

> [!key] Definition
> $$
> \mathbf{E}(\mathbf{r}) = \lim_{q\to 0}\frac{\mathbf{F}}{q}\qquad [\text{V/m} = \text{N/C}]
> $$
> The force a small stationary test charge $q$ would feel at $\mathbf{r}$, divided by $q$. It exists whether or not a test charge is there.

**How it is produced.** By charge: a point charge gives the Coulomb field $\dfrac{Q}{4\pi\epsilon_0 r^2}\hat{r}$ ([[concepts/coulombs-law]]); many charges give the vector sum ([[concepts/superposition]]); charge densities give integrals. Time-varying magnetic fields *also* produce $\mathbf{E}$ (Faraday, Lecture 14) — that part of $\mathbf{E}$ is not curl-free.

**How to picture it.**
- *Arrow plots*: length ∝ $|\mathbf{E}|$ at grid points.
- *Field lines*: tangent to $\mathbf{E}$ everywhere; start on + charge, end on − charge or at infinity; density ∝ $|\mathbf{E}|$; never cross (one direction per point), except at null points where $\mathbf{E}=0$. Twice the charge, twice the lines.

**Its two local properties (the laws of electrostatics).**
$$
\nabla\cdot\mathbf{E} = \rho/\epsilon_0 \quad(\text{sources are charges}),\qquad \nabla\times\mathbf{E} = 0 \quad(\text{static: no whirlpools}).
$$
The second says a static $\mathbf{E}$ is [[concepts/conservative-field|conservative]], so $\mathbf{E} = -\nabla V$ — voltage exists ([[1-electrostatics/05-curl-free-fields-and-the-electrostatic-potential|Lecture 5]], [[concepts/electrostatic-potential]]).

**Relatives.** $\mathbf{D} = \epsilon_0\mathbf{E}$ in vacuum, $\epsilon_0\mathbf{E}+\mathbf{P}$ in matter ([[concepts/electric-flux-density]]). Inside a perfect conductor in statics, $\mathbf{E}=0$ ([[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]], [[concepts/conductors]]).

**Where it appears.** [[1-electrostatics/01-fields-forces-and-the-maxwell-roadmap|Lecture 1]] · [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]] · every lecture after.

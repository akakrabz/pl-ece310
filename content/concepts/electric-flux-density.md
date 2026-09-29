---
title: "Electric flux density D"
description: "D = ε₀E in vacuum, D = ε₀E + P in matter: the field whose flux counts free charge in coulombs."
tags: [concept, electrostatics, exam-1]
aliases: ["D field", "displacement field", "electric displacement"]
---

> [!key] Definition
> $$
> \mathbf{D} = \epsilon_0\mathbf{E}\ \ (\text{vacuum}),\qquad
> \mathbf{D} = \epsilon_0\mathbf{E}+\mathbf{P} = \epsilon\mathbf{E}\ \ (\text{linear dielectric, } \epsilon = \epsilon_0\epsilon_r)
> \qquad[\text{C/m}^2]
> $$
> Defined so that **Gauss's law counts charge in coulombs**: $\oint_S\mathbf{D}\cdot d\mathbf{S} = Q_{\text{enc}}$. Maxwell called $\mathbf{D}$ the *electric displacement*; the slides say "displacement flux density"; this site says "electric flux density".

**Why two fields?** In vacuum $\mathbf{D}$ is just $\epsilon_0\mathbf{E}$ and carries no new information. Inside matter ([[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]], [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]]) the medium polarizes, producing *bound* charge that also sources $\mathbf{E}$. $\mathbf{D}$ is built so that its sources are the **free** charge only:
$$
\nabla\cdot\mathbf{D} = \rho_{\text{free}},\qquad \nabla\cdot(\epsilon_0\mathbf{E}) = \rho_{\text{free}}+\rho_{\text{bound}},\qquad \rho_{\text{bound}} = -\nabla\cdot\mathbf{P}.
$$
So a Gauss's-law calculation gives $\mathbf{D}$ from the free charge *regardless of the material*, and then $\mathbf{E} = \mathbf{D}/\epsilon$ region by region. That is exactly the structure of the FA26 coax problem (problem 4a–b): one $\mathbf{D}$, two different $\mathbf{E}$'s.

**Boundary behaviour** ([[1-electrostatics/06-circulation-and-boundary-conditions|Lecture 6]], [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]]; [[concepts/boundary-conditions]]): the normal component of $\mathbf{D}$ is continuous across an interface carrying no free surface charge; the tangential component of $\mathbf{E}$ is continuous. (Not the other way round — FA26 problem 3.)

**Units check.** $\epsilon_0$ [F/m] × $\mathbf{E}$ [V/m] = C/(V·m)·(V/m) = C/m² ✓.

**Where it appears.** [[1-electrostatics/02-coulombs-law-superposition-and-gauss#5-from-coulomb-to-gauss|Lecture 2 §5]] (introduced), [[1-electrostatics/03-gauss-law-at-work|Lecture 3]] (every Gauss's-law example), [[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]], [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]], [[1-electrostatics/10-capacitance-and-conductance|Lecture 10]] (materials, capacitance), [[problems/two-layer-coaxial-capacitor]].

Related: [[concepts/gauss-law]] · [[concepts/electric-field]] · [[concepts/flux]] · [[concepts/polarization]] · [[concepts/permittivity]].

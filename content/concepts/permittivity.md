---
title: "Permittivity, susceptibility, and D = εE"
description: "For linear dielectrics P = ε₀χeE, so D = ε₀(1 + χe)E = εE with ε = εrε₀. One number replaces the whole polarization story; typical εr values; what changes and what does not when a dielectric is inserted."
tags: [concept, electrostatics, exam-1]
aliases: ["dielectric constant", "relative permittivity", "electric susceptibility", "constitutive relation", "linear dielectric"]
---

> [!key] Constitutive relation of a linear, isotropic dielectric
> $$
> \mathbf{P} = \epsilon_0\chi_e\mathbf{E}\quad\Longrightarrow\quad \mathbf{D} = \epsilon_0(1+\chi_e)\mathbf{E} = \epsilon\mathbf{E},\qquad
> \epsilon = \epsilon_r\epsilon_0,\quad \epsilon_r = 1+\chi_e .
> $$
> $\chi_e\ge0$ (electric susceptibility) and $\epsilon_r$ (relative permittivity, "dielectric constant") are dimensionless; $\epsilon$ is in F/m. $\epsilon_0 = 8.854\times10^{-12}$ F/m $\approx 10^{-9}/36\pi$. Typical $\epsilon_r$: vacuum 1, air 1.0006, glass 4–10, soil 5–10, silicon 11–12, distilled water 81.

**Where the ε goes.** Everywhere $\epsilon_0$ appeared in a vacuum formula whose field region is filled with the dielectric: Coulomb's law and the point-charge potential ($Q/4\pi\epsilon r$), Poisson's equation ($-\rho/\epsilon$), capacitance ($\epsilon A/d$, $2\pi\epsilon/\ln(b/a)$), energy density ($\tfrac12\epsilon E^2$). Gauss's law is best used in its $\mathbf{D}$ form, $\oint\mathbf{D}\cdot d\mathbf{S} = Q_{\text{free}}$, which needs no $\epsilon$ at all.

**What changes when a dielectric is inserted.** With the *free charge held fixed*, $\mathbf{D}$ is unchanged and $\mathbf{E}$ drops by $\epsilon_r$ (so does $V$; $C$ rises by $\epsilon_r$). With the *voltage held fixed*, $\mathbf{E}$ is unchanged and $\mathbf{D}$, $\rho_s$ and $Q$ rise by $\epsilon_r$ ($C$ rises by $\epsilon_r$ either way). $\mathbf{P} = \mathbf{D}-\epsilon_0\mathbf{E} = (1-1/\epsilon_r)\mathbf{D}$ makes up the difference. Limits: $\chi_e = 0$ is vacuum; $\chi_e\to\infty$ makes the material behave like a conductor ($\mathbf{E}\to0$).

**At an interface** between $\epsilon_1$ and $\epsilon_2$ (no free charge): $E_t$ continuous, $D_n$ continuous, so $E_{2n} = (\epsilon_1/\epsilon_2)E_{1n}$ and $D_{2t} = (\epsilon_2/\epsilon_1)D_{1t}$; $\tan\theta_1/\tan\theta_2 = \epsilon_1/\epsilon_2$.

**Homogeneous vs inhomogeneous.** Laplace's and Poisson's equations hold only where $\epsilon_r$ is independent of position. Layered media: solve each layer and match $V$ and $D_n$. Graded $\epsilon(\mathbf{r})$: use $\nabla\cdot\mathbf{D} = \rho$ directly ($D$ constant in 1-D with no free charge, then $E = D/\epsilon(z)$). Anisotropic media ($\epsilon$ a tensor) are not treated in this course; plasmas ($\chi_e<0$) are deferred to ECE 350.

> [!trap]
> - $\mathbf{D} = \epsilon\mathbf{E}$ uses the permittivity *of the region where you evaluate it*; the FA26 key's struck-out $-10\epsilon_0$ (instead of $2\epsilon_0\cdot(-10)$) is this slip.
> - $\epsilon_r$ and $\chi_e$ differ by 1; $\epsilon_r = 6$ means $\chi_e = 5$.
> - "$\epsilon\oint\mathbf{E}\cdot d\mathbf{S} = Q$" holds only when one uniform $\epsilon$ fills the whole surface; $\oint\mathbf{D}\cdot d\mathbf{S} = Q_{\text{free}}$ always does.

**Where it appears.** [[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]], [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]], [[1-electrostatics/10-capacitance-and-conductance|Lecture 10]], [[problems/fields-across-a-dielectric-interface]], [[problems/two-layer-coaxial-capacitor]].

Related: [[concepts/polarization]] · [[concepts/electric-flux-density]] · [[concepts/boundary-conditions]] · [[concepts/capacitance]].

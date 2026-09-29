---
title: "Worked problem — Fields across a dielectric interface"
description: "Given E just above a plane between two dielectrics, find E just below, then D and P on both sides, and the bound charge on the interface. Modelled on FA26 Exam 1, problem 3, with different permittivities and a three-component field."
tags: [problem, electrostatics, exam-1]
---

*Problem · style of FA26 Exam 1 #3 (24 pts) · uses [[1-electrostatics/06-circulation-and-boundary-conditions|Lecture 6]], [[1-electrostatics/08-conductors-dielectrics-and-polarization|Lecture 8]], [[1-electrostatics/09-static-fields-in-dielectric-media|Lecture 9]] · concepts: [[concepts/boundary-conditions]], [[concepts/permittivity]], [[concepts/polarization]]*

> [!question] Problem
> Two perfect dielectrics meet at the plane $z=0$: $\epsilon_1 = 4\epsilon_0$ for $z>0$ and $\epsilon_2 = 6\epsilon_0$ for $z<0$. There is no free charge on the interface. Immediately above it the field is measured to be
> $$
> \mathbf{E}_1 = 3\hat{x} + 2\hat{y} - 6\hat{z}\quad[\text{V/m}].
> $$
> (a) Find $\mathbf{E}_2$, the field immediately below the interface. (b) Find $\mathbf{D}_1$ and $\mathbf{D}_2$. (c) Find $\mathbf{P}_1$ and $\mathbf{P}_2$. (d) Find the bound surface charge density on the interface, and the angle each field makes with the normal.

## Before writing anything: which component is which

The interface is the plane $z=0$, so the **normal** direction is $\hat{z}$ and the **tangential** directions are $\hat{x}$ and $\hat{y}$. Take medium 1 to be the upper half-space; then $\hat{n} = +\hat{z}$ points from medium 2 into medium 1, as the boundary conditions require. Split the given field:

$$
\mathbf{E}_{1t} = 3\hat{x}+2\hat{y},\qquad E_{1n} = -6 .
$$

Two rules govern everything ([[concepts/boundary-conditions]]): tangential $\mathbf{E}$ is continuous, and — because the interface carries no free charge — normal $\mathbf{D}$ is continuous.

## (a) The field below

**Tangential:** $\mathbf{E}_{2t} = \mathbf{E}_{1t} = 3\hat{x}+2\hat{y}$. Copied, no arithmetic.

**Normal:** $D_{1n} = D_{2n}$ with $D_n = \epsilon E_n$ in each medium:

$$
\epsilon_1E_{1n} = \epsilon_2E_{2n}\quad\Longrightarrow\quad E_{2n} = \frac{\epsilon_1}{\epsilon_2}E_{1n} = \frac{4\epsilon_0}{6\epsilon_0}(-6) = -4 .
$$

> [!key] Answer (a)
> $\mathbf{E}_2 = 3\hat{x}+2\hat{y}-4\hat{z}$ [V/m]. Entering the higher-permittivity medium, the normal component shrinks by $\epsilon_1/\epsilon_2 = 2/3$; the tangential part is untouched.

> [!trap] The four wrong answers, and the rule each one breaks
> - $\mathbf{E}_2 = \mathbf{E}_1$: normal $\mathbf{E}$ treated as continuous. It is normal $\mathbf{D}$ that is continuous.
> - $\mathbf{E}_2 = \tfrac23\mathbf{E}_1 = 2\hat{x}+\tfrac43\hat{y}-4\hat{z}$: the whole vector scaled, i.e. tangential $\mathbf{D}$ treated as continuous.
> - $E_{2n} = \tfrac32(-6) = -9$: the ratio inverted. Sanity check: the field is *weaker* in the higher-$\epsilon$ medium (its bound charge opposes the field), so $\lvert E_{2n}\rvert<\lvert E_{1n}\rvert$.
> - $D_{1n} = \epsilon_0(-6)$ instead of $4\epsilon_0(-6)$: $\mathbf{D} = \epsilon\mathbf{E}$ must use the permittivity of the medium you are in. The FA26 key itself contains this slip, struck out.

## (b) The displacement fields

$\mathbf{D} = \epsilon\mathbf{E}$ in each medium, each with its own $\epsilon$:

$$
\mathbf{D}_1 = 4\epsilon_0(3\hat{x}+2\hat{y}-6\hat{z}) = \epsilon_0\,(12\hat{x}+8\hat{y}-24\hat{z}),\qquad
\mathbf{D}_2 = 6\epsilon_0(3\hat{x}+2\hat{y}-4\hat{z}) = \epsilon_0\,(18\hat{x}+12\hat{y}-24\hat{z}).
$$

> [!key] Answer (b)
> $\mathbf{D}_1 = \epsilon_0(12\hat{x}+8\hat{y}-24\hat{z})$ and $\mathbf{D}_2 = \epsilon_0(18\hat{x}+12\hat{y}-24\hat{z})$ [C/m²]; numerically $\epsilon_0 = 8.854\times10^{-12}$ F/m, so $D_{1z} = D_{2z}\approx-2.1\times10^{-10}$ C/m².

**Built-in check.** The $z$-components agree ($-24\epsilon_0$ on both sides) — they must, since that was the condition used in (a). The tangential components of $\mathbf{D}$ *do* jump, from $(12,8)\epsilon_0$ to $(18,12)\epsilon_0$, by the factor $\epsilon_2/\epsilon_1 = 3/2$. If your two $D_z$ values differ, part (a) is wrong; fix it before going on.

## (c) The polarization

$\mathbf{P} = \mathbf{D}-\epsilon_0\mathbf{E}$, medium by medium, pairing each $\mathbf{D}$ with the $\mathbf{E}$ of the *same* medium:

$$
\mathbf{P}_1 = \epsilon_0(12\hat{x}+8\hat{y}-24\hat{z}) - \epsilon_0(3\hat{x}+2\hat{y}-6\hat{z}) = \epsilon_0(9\hat{x}+6\hat{y}-18\hat{z}),
$$
$$
\mathbf{P}_2 = \epsilon_0(18\hat{x}+12\hat{y}-24\hat{z}) - \epsilon_0(3\hat{x}+2\hat{y}-4\hat{z}) = \epsilon_0(15\hat{x}+10\hat{y}-20\hat{z}).
$$

> [!key] Answer (c)
> $\mathbf{P}_1 = \epsilon_0(9\hat{x}+6\hat{y}-18\hat{z})$, $\mathbf{P}_2 = \epsilon_0(15\hat{x}+10\hat{y}-20\hat{z})$ [C/m²].

**Shortcut and check.** $\mathbf{P} = (\epsilon_r-1)\epsilon_0\mathbf{E} = \epsilon_0\chi_e\mathbf{E}$: with $\chi_{e1} = 3$, $\mathbf{P}_1 = 3\epsilon_0\mathbf{E}_1 = \epsilon_0(9,6,-18)$ ✓; with $\chi_{e2} = 5$, $\mathbf{P}_2 = 5\epsilon_0\mathbf{E}_2 = \epsilon_0(15,10,-20)$ ✓. Each $\mathbf{P}$ is parallel to its own $\mathbf{E}$, as it must be in an isotropic medium.

> [!trap] Two lines that are always wrong
> - $\mathbf{P} = \mathbf{D}-\epsilon\mathbf{E}$ is identically zero. The vacuum part $\epsilon_0\mathbf{E}$ is what gets subtracted.
> - $\mathbf{P} = \chi_e\mathbf{E}$ has the wrong units (V/m, not C/m²). The $\epsilon_0$ belongs there.
> - And a subtler one: $\mathbf{D}_2 - \epsilon_0\mathbf{E}_1$, pairing fields from different media, gives a meaningless vector.

## (d) Bound charge on the interface, and the bending of the field

There is no *free* charge on $z=0$, but the normal component of $\mathbf{P}$ jumps there, so the interface carries **bound** surface charge. With $\hat{n} = \hat{z}$ from 2 into 1, $\hat{n}\cdot(\mathbf{P}_1-\mathbf{P}_2) = -\rho_{sb}$:

$$
\rho_{sb} = P_{2z}-P_{1z} = (-20+18)\epsilon_0 = -2\epsilon_0\ \text{C/m}^2\approx-1.8\times10^{-11}\ \text{C/m}^2 .
$$

Cross-check with Gauss's law for *all* charge: $\epsilon_0(E_{1z}-E_{2z}) = \rho_s+\rho_{sb} = 0 + \rho_{sb}$, and indeed $\epsilon_0(-6-(-4)) = -2\epsilon_0$ ✓. The negative bound layer is what takes the $z$-component of $\mathbf{E}$ from $-6$ to $-4$ (weaker) across the interface.

The **angles** from the normal: $\tan\theta_1 = \lvert\mathbf{E}_{1t}\rvert/\lvert E_{1n}\rvert = \sqrt{13}/6$, $\tan\theta_2 = \sqrt{13}/4$, so $\theta_1\approx31.0°$, $\theta_2\approx42.0°$, and $\tan\theta_1/\tan\theta_2 = 4/6 = \epsilon_1/\epsilon_2$ ✓ — the refraction law of [[1-electrostatics/09-static-fields-in-dielectric-media#2-refraction-of-field-lines-at-an-interface|Lecture 9 §2]]. The line bends *away* from the normal on entering the higher-$\epsilon$ side.

## Variant: with free charge on the interface

Suppose the same interface also carried free charge $\rho_s = 6\epsilon_0$ C/m². The tangential rule is unchanged. The normal rule becomes $\hat{n}\cdot(\mathbf{D}_1-\mathbf{D}_2) = \rho_s$, i.e. $D_{2z} = D_{1z}-\rho_s = -24\epsilon_0 - 6\epsilon_0 = -30\epsilon_0$, so $E_{2z} = -30\epsilon_0/6\epsilon_0 = -5$ and $\mathbf{E}_2 = 3\hat{x}+2\hat{y}-5\hat{z}$ V/m. In words: crossing the sheet *downward* (against $\hat{n}$), $D_z$ drops by $\rho_s$ — a positive sheet pushes field lines away from itself on both sides, making $E_z$ more negative below and less negative above.

## The general pattern

For $\epsilon_1 = \epsilon_{r1}\epsilon_0$ above and $\epsilon_2 = \epsilon_{r2}\epsilon_0$ below, with free charge $\rho_s$ on the plane and $\hat{n} = \hat{z}$:

$$
\mathbf{E}_{2t} = \mathbf{E}_{1t},\qquad
E_{2n} = \frac{\epsilon_{r1}E_{1n} - \rho_s/\epsilon_0}{\epsilon_{r2}},\qquad
\mathbf{P}_i = (\epsilon_{ri}-1)\epsilon_0\mathbf{E}_i,\qquad
\rho_{sb} = P_{2n}-P_{1n} = \epsilon_0(E_{1n}-E_{2n}) - \rho_s .
$$

The exam instance was $\epsilon_{r1} = 2$, $\epsilon_{r2} = 5$, $\mathbf{E}_1 = 4\hat{x}-10\hat{z}$, giving $\mathbf{E}_2 = 4\hat{x}-4\hat{z}$ and $\rho_{sb} = -6\epsilon_0$. Whatever the numbers, the structure is: split into tangential and normal, copy the tangential part, scale the normal part by $\epsilon_1/\epsilon_2$ (and subtract $\rho_s/\epsilon_2$ if there is free charge), then $\mathbf{D} = \epsilon\mathbf{E}$ and $\mathbf{P} = \mathbf{D}-\epsilon_0\mathbf{E}$ medium by medium. If the interface is $x=0$ instead, the normal component is $E_x$ — the geometry decides which component gets scaled, not the letter.

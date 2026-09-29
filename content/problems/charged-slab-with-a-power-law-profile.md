---
title: "Worked problem — Charged slab with a non-uniform profile"
description: "Total charge per area, E inside and outside by Gauss's law with a symmetric pillbox, and the potential at a point outside. Modelled on FA26 Exam 1, problem 2, with a quadratic profile instead of the exam's |x| profile."
tags: [problem, electrostatics, exam-1]
---

*Problem · style of FA26 Exam 1 #2 (21 pts) · uses [[1-electrostatics/03-gauss-law-at-work|Lecture 3]] and Lecture 5 · concepts: [[concepts/gauss-law]], [[concepts/charge-density]], [[concepts/conservative-field]]*

> [!question] Problem
> An infinite slab occupies $-a<x<a$ and extends without limit in $y$ and $z$. Its charge density is
> $$
> \rho(x) = \begin{cases}\rho_0\,\dfrac{x^2}{a^2}, & |x|<a\\[6pt] 0, & |x|>a\end{cases}\qquad[\text{C/m}^3],\quad \rho_0>0 .
> $$
> (a) Find the total charge per unit area (per unit of $yz$-plane area). (b) Using Gauss's law, find $\mathbf{E}(x)$ for $|x|<a$ and for $|x|>a$. (c) Taking $V(0)=0$, find $V(b)$ for a point $b>a$.

## (a) Charge per unit area

Integrate through the thickness. The profile is even in $x$, so integrate half and double:

$$
\rho_s^{\text{tot}} = \int_{-a}^{a}\rho_0\frac{x^2}{a^2}\,dx = \frac{2\rho_0}{a^2}\int_0^a x^2\,dx = \frac{2\rho_0}{a^2}\cdot\frac{a^3}{3} = \frac{2\rho_0a}{3}\quad[\text{C/m}^2].
$$

> [!key] Answer (a)
> $\rho_s^{\text{tot}} = \dfrac{2}{3}\rho_0 a$ [C/m²]. (Compare: a uniform slab of the same width and peak density would hold $2\rho_0a$; the parabolic profile holds a third of that.)

## (b) The field, by Gauss's law

**Symmetry.** $\rho$ depends only on $x$ and is even, so $\mathbf{E} = E_x(x)\,\hat{x}$ with $E_x$ *odd*: $E_x(-x) = -E_x(x)$, hence $E_x(0) = 0$. Field points away from the slab on both sides.

**Surface.** A pillbox with caps of area $A$ at $\pm x$, straddling the symmetry plane. Sides: $\mathbf{D}\parallel$ surface, zero flux. Both caps: outward normal and field both point away from the centre, each contributing $D_x(x)\,A$.

**Inside** ($0<x<a$):

$$
2\,D_x(x)\,A = Q_{\text{enc}} = A\int_{-x}^{x}\rho_0\frac{x'^2}{a^2}\,dx' = A\,\frac{2\rho_0x^3}{3a^2}
\;\Rightarrow\; E_x = \frac{\rho_0\,x^3}{3\epsilon_0a^2}.
$$

**Outside** ($x>a$): the pillbox encloses the whole slab, $Q_{\text{enc}} = A\cdot\frac{2\rho_0a}{3}$ from part (a):

$$
2\,D_x\,A = \frac{2\rho_0a}{3}A \;\Rightarrow\; E_x = \frac{\rho_0a}{3\epsilon_0}.
$$

> [!key] Answer (b)
> $$
> \mathbf{E}(x) = \begin{cases}
> \hat{x}\,\dfrac{\rho_0\,x^3}{3\epsilon_0 a^2}, & |x|\le a\\[8pt]
> \hat{x}\,\text{sgn}(x)\,\dfrac{\rho_0 a}{3\epsilon_0}, & |x|\ge a
> \end{cases}\qquad[\text{V/m}]
> $$
> (The inside expression is already odd in $x$, so it needs no sgn; the outside one does.)

**Three checks that cost ten seconds.**
1. Continuity at $x=a$: inside gives $\rho_0a/(3\epsilon_0)$, outside gives the same ✓ (no surface charge at the faces, so $E$ must be continuous).
2. Differential Gauss: $\epsilon_0\,dE_x/dx = \rho_0x^2/a^2 = \rho(x)$ ✓ inside; $0$ outside ✓.
3. Far away the slab is a sheet with $\rho_s = \frac{2}{3}\rho_0a$, whose field is $\rho_s/(2\epsilon_0) = \rho_0a/(3\epsilon_0)$ ✓ — the outside field is exactly the sheet result.

> [!trap] Where the points go missing
> - Copying the *uniform*-slab formula $\rho x/\epsilon_0$: here $\rho$ varies, so $Q_{\text{enc}}$ is an integral and $E_x\propto x^3$.
> - One cap instead of two → a factor 2 wrong outside and inside. (A one-cap pillbox is fine only if the other cap sits at $x=0$ where $E=0$ — say so.)
> - No direction for $x<0$.

## (c) The potential outside

$V(b) - V(0) = -\int_0^b E_x\,dx$, and the integrand changes form at $x=a$, so split there:

$$
V(a) = -\int_0^a\frac{\rho_0x^3}{3\epsilon_0a^2}\,dx = -\frac{\rho_0}{3\epsilon_0a^2}\cdot\frac{a^4}{4} = -\frac{\rho_0a^2}{12\epsilon_0},
\qquad
V(b)-V(a) = -\int_a^b\frac{\rho_0a}{3\epsilon_0}\,dx = -\frac{\rho_0a\,(b-a)}{3\epsilon_0}.
$$

Add:

$$
V(b) = -\frac{\rho_0a^2}{12\epsilon_0} - \frac{\rho_0ab}{3\epsilon_0} + \frac{\rho_0a^2}{3\epsilon_0}
= \frac{\rho_0a^2}{4\epsilon_0} - \frac{\rho_0ab}{3\epsilon_0}
= \frac{\rho_0a\,(3a - 4b)}{12\epsilon_0}.
$$

> [!key] Answer (c)
> $V(b) = \dfrac{\rho_0a^2}{4\epsilon_0} - \dfrac{\rho_0ab}{3\epsilon_0}$ [V] for $b\ge a$ — negative (the potential falls as you move away from positive charge). $V(x)$ itself is even, so $V(-b) = V(b)$: replace $b$ by $|b|$ in the formula.

Sanity: at $b=a$ this gives $-\rho_0a^2/(12\epsilon_0) = V(a)$ ✓. For $b = 2a$: $V = -\frac{5}{12}\rho_0a^2/\epsilon_0$.

> [!trap]
> - Using the outside field all the way from 0 gives $-\rho_0ab/(3\epsilon_0)$ — wrong by the constant $\rho_0a^2/(4\epsilon_0)$.
> - Dropping the minus sign in $V = -\int\mathbf{E}\cdot d\mathbf{l}$ gives a positive potential outside a positive slab.
> - Trying to reference $V$ at infinity: for an *infinite* slab $V\to-\infty$ there, which is why the problem supplies $V(0)=0$.

## The general pattern

For $\rho = \rho_0|x/a|^k$ on $|x|<a$: total $\dfrac{2\rho_0a}{k+1}$; inside $E_x = \dfrac{\rho_0\,x|x|^{k}}{(k+1)\epsilon_0a^k}$; outside $\pm\dfrac{\rho_0a}{(k+1)\epsilon_0}$; $V(b) = \dfrac{\rho_0a^2}{(k+2)\epsilon_0} - \dfrac{\rho_0ab}{(k+1)\epsilon_0}$ for $b\ge a$. The exam used $k=1$; this page used $k=2$; $k=0$ is the uniform slab of [[1-electrostatics/03-gauss-law-at-work#3c-uniform-slab-planar-with-interior|Lecture 3 §3c]]. A slab with an *odd* profile (no absolute value) is a different animal: net charge zero, $\mathbf{E}=0$ outside, and the symmetric pillbox no longer works because $E_x$ is then even — integrate $dE_x/dx = \rho/\epsilon_0$ from a face instead.

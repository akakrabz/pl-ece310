---
title: "Errata in the course materials"
description: "Slips found in the lecture slides, the course notes, the FA26 Exam 1 key and the homework solutions while writing these notes — each with the correct statement. None of them is propagated into this site; several are the exact mistakes students make."
tags: [toolkit, exam-1]
aliases: ["errata", "known errors in the slides", "errors in the exam key"]
---

Course materials are written fast and by hand, and every set has a few slips. This page lists the ones found while digesting the sources for Lectures 1–10, so that you (a) do not copy them onto your notecard and (b) recognize them for what they are when a formula on a slide looks wrong. Only substantive items are listed — things that would change an answer or a sign — not spelling. Each entry gives the source, what it says, and what it should say.

> [!tip] How to read this page
> "Slide" means the annotated lecture decks (Shao, after Goddard and Cunningham); "notes" means Prof. Kudeki's lecture notes; "key" means the FA26 Exam 1 solutions; "HW sol." means the FA26 homework solution sets. Nothing here is a criticism of the materials — most of these are momentary hand slips that the authors corrected verbally or on the next line. They are listed because they are copied verbatim into notecards every semester.

## Lecture slides

| where | as written | should be | why it matters |
|---|---|---|---|
| L1 slides 25, 28 | $\vec A\times\vec B = \lvert\vec A\rvert\lvert\vec B\rvert\sin\theta$ | $\lvert\vec A\times\vec B\rvert = \lvert\vec A\rvert\lvert\vec B\rvert\sin\theta$ | the cross product is a vector; only its magnitude is that expression |
| L3 slide 17 | line charge as $\lambda\,\delta(x-x_0)$ | $\lambda\,\delta(x-x_0)\,\delta(y-y_0)$ | a line needs two δ's to have units C/m³ (slide 16's point charge was corrected by hand to three δ's) |
| L3 slide 20 | $\vec D = \frac{\rho_s}{2}\,\text{sgn}(z)$ | $\vec D = \hat z\,\frac{\rho_s}{2}\,\text{sgn}(z)$ | the unit vector is missing; slide 22 then places the same sheet at $x=0$ |
| L4 slides 13–14, pp. 26–28 | $d\vec B/dt$, $d\vec D/dt$, $d\rho/dt$ in the differential Maxwell equations | $\partial\vec B/\partial t$, etc. | fields depend on $(x,y,z,t)$; the time derivative is partial |
| L5 slide 3 | "Gradient of $E$: $\nabla V\equiv\ldots$" | gradient of $V$ | one takes the gradient of the scalar |
| L5 slide 18 | "$\vec E = -\nabla\Phi$ satisfies Gauss's law by the identity $\nabla\cdot\nabla\Phi = 0$" | the identity is $\nabla\times\nabla\Phi = 0$, which makes $-\nabla\Phi$ satisfy $\nabla\times\vec E = 0$ | $\nabla\cdot\nabla\Phi = \nabla^2\Phi$ is the Laplacian, and Gauss's law turns it into Poisson's equation, not zero ([[1-electrostatics/07-poisson-and-laplace\|Lecture 7]]) |
| L6–7 slide 5 | $D_1 - D_2 = \rho$ | $D_{1n} - D_{2n} = \rho_s$ | the jump is set by *surface* charge [C/m²]; the boxed equation on the same slide is right |
| L6–7 slide 12 (ink) | $J_s = I/S$, $I = \int_S J_s\,d\vec S$ | these are the *volume* current relations ($\vec J$ in A/m²); a surface current density $\vec J_s$ is current per unit *width* [A/m], $I = \int\vec J_s\cdot\hat u\,dl$ across a line | units of the Ampère boundary condition $\hat n\times(\vec H_1-\vec H_2) = \vec J_s$ |
| L6–7 slide 13 (ink) | "$D_{2n} = 5\hat y$" | $D_{1n} = 5$ | the unknown was $D_{1n}$; the circled answer (d) is right |
| L8 slide 5 vs 14 | $\rho$ used for resistivity [Ω·m] and then for charge density | resistivity $= 1/\sigma$ | the same letter for two quantities on one deck |
| L8 slide 15 | header "$H_{\text{cond}} = 0$" presented as a boundary fact | a course *assumption* (slide 17 says so) | a static $\vec H$ can exist inside a perfect conductor; the course chooses to take it as zero |
| L9 slide 7, 13 | $\vec P$ "due to the total $\vec E_a$ field" / "under applied $\vec E$" | $\vec P = \epsilon_0\chi_e\vec E_{\text{tot}}$, the **total** field (corrected in ink on slide 7) | using the applied field alone overestimates $\vec P$ by $(1+\chi_e)$ |
| L9 slide 17 (ink) | $V = -\nabla E$ | $\vec E = -\nabla V$ | inverted; the circled answer (a) is still right |
| L10 slide 3 | $V(z) = V_0x/d$; $\rho_s = \vec a_n\cdot(D_{n1}-D_{n2})$ | $V(x) = V_0x/d$; $\hat n\cdot(\vec D_1-\vec D_2) = \rho_s$ | mismatched coordinate; a vector dotted into scalar components |
| L10 slide 4 | $V = \int_0^d\vec E\cdot d\vec l$ | $V(d)-V(0) = -\int_0^d\vec E\cdot d\vec l$, or "the drop from + to −" | as written it gives a negative plate voltage for the stated field; the pen "= $E\,d$" silently uses magnitudes |
| L10 slide 7 | Gauss on a cylinder "of radius $r>a$" | $a<r<b$ | for $r>b$ the enclosed charge is zero |
| L10 slide 10 | "work by the battery … $dW = V_c\,dq$" | the battery supplies $V\,dq$; $V_c\,dq$ is what the *capacitor* stores | over a full charge the battery supplies $CV^2$ and half is dissipated ([[1-electrostatics/10-capacitance-and-conductance\|Lecture 10 §5]]) |
| L10 slide 11 (ink) | energy density "[W/m³]" | J/m³ | watts are power |
| Exam 1 review, slide 21 | "$\vec E = \nabla V$ for some scalar $V$" | $\vec E = -\nabla V$ | the course convention (the same deck's slide 22 and equation sheet have the minus sign) |
| Exam 1 review, slide 24 (ink) | $\hat n\cdot(\vec P_1-\vec P_2) = \rho_{sb}$ | $\hat n\cdot(\vec P_1-\vec P_2) = -\rho_{sb}$ | the minus sign is what makes $\epsilon_0E_n$ jump by $\rho_s+\rho_{sb}$ while $D_n$ jumps by $\rho_s$ ([[concepts/polarization]]) |
| Exam 1 review, slides 13–15 and equation sheet | $\epsilon\oint\vec E\cdot d\vec S = Q_{\text{enc}}$ | valid only when one uniform $\epsilon$ fills the region; in general $\oint\vec D\cdot d\vec S = Q_{\text{free}}$ | with two dielectrics inside the surface there is no single $\epsilon$ to pull out |

## Course notes (Kudeki)

| where | as written | should be |
|---|---|---|
| L5 p. 7 | "$V(x,y,x) = -x^2-3(y+1)z$" | $V(x,y,z)$ |
| L5 p. 1 | determinant expansion written $\hat x0-\hat y0-\hat z0$ | the cofactor signs are $+,-,+$ (harmless here, every term is zero) |
| L5 | Examples are numbered 1, 2, 3, 5, 6 | there is no Example 4; nothing is missing |
| L6 p. 9 | "$D_x$ is zero in $z<0$" | in $x<0$ |
| L7 p. 6 | "in terms of $V_{12}$" | $V_{21}$, as defined two lines earlier |
| L8 p. 12 | boundary condition "relevant for $\vec D = \epsilon\vec E+\vec P$" | $\vec D = \epsilon_0\vec E+\vec P$ |
| L9 pp. 4, 7 (margin graphs) | $V(z)$ for the two-layer plates, and $\epsilon(z) = 4\epsilon_0/(4-z)$, drawn as straight lines | qualitative sketches: the actual rise is 2 V then 1 V (not the other way round), and $\epsilon(z)$ is convex |
| L10 p. 7 | coax length written "$l$" once | $\ell$, as everywhere else |

## FA26 Exam 1 solution key

None of these changes a boxed answer; all of them are exactly the slips graders take points for.

| part | as written | should be |
|---|---|---|
| 1(a) | $\hat z$-component of the curl written $\big(\partial_yE_x - \partial_xE_y\big)$ | $\big(\partial_xE_y - \partial_yE_x\big)$ — harmless only because both partials are equal for this field; copied into a field with nonzero curl it flips the sign |
| 3(a) | $D_{1z} = -10\epsilon_0$ (struck out by the grader) | $D_{1z} = \epsilon_1E_{1z} = 2\epsilon_0\cdot(-10)$: use the permittivity of the medium you are in |
| 3(a) | a stray "$4\epsilon$" in $E_{2z}$ (struck out) | $E_{2z} = D_{2z}/\epsilon_2 = -4$, no $\epsilon$ left over |
| 4(e) | $\mathcal C = Q/V$ evaluated with $V = V(b)-V(a)<0$, giving a negative $\mathcal C$ in the intermediate line; the box adds the minus sign silently | $\mathcal C = \lambda/[V(a)-V(b)]$: the drop from the positive conductor to the negative one ([[problems/two-layer-coaxial-capacitor]]) |

## Homework solution sets (FA26)

| where | as written | should be |
|---|---|---|
| HW3 #1(d) | "$\rho_s = \vec D\cdot\hat n = D_x\big\rvert_{x=d}$" | $\hat n = -\hat x$ at the anode, so $\rho_s = -D_x(d) = +8\epsilon_0/3d$; the final value is right, the intermediate sign is not |
| HW3 #5 (hint) | $d\vec l = (-\hat x-\hat y)\,dx$ *and* limits $\int_1^{-1}$ | one or the other: reversing the direction in $d\vec l$ and in the limits counts it twice (harmless here because the slant-edge integrands vanish) |
| HW3 #6(e) | "$\int 4\,dz$" | $\int 4\,dy$ |
| HW3 #7(f) | bound surface charges left as functions of $r$ | evaluate at $r=a$ (inner face) and $r=b$ (outer face) |
| HW4 #4(e) | general solution "$a_1e^{-Gt/C}+a_2$" | $a_2 = 0$ is forced by the ODE itself, not by a boundary condition |
| HW4 #6 | a stray label $\vec P_4$; only *free* surface charges reported | the four surfaces also carry bound charge (listed in [[1-electrostatics/09-static-fields-in-dielectric-media\|Lecture 9]]'s pattern: $\rho_{sb} = \vec P\cdot\hat n$) |

## Notation that varies between sources (not errors, but easy to trip on)

- Line charge density: $\lambda$, $\lambda_0$, $\rho_l$, $\rho_\ell$. This site uses $\rho_l$, and $\lambda$ only where the exam does.
- Unit normal at an interface: $\hat a_n$ (slides, "points into medium 1"), $\hat n$ (ink), and Kudeki's $\pm$ superscripts ($\hat n$ from $-$ to $+$). Same content; see [[concepts/boundary-conditions]].
- Component subscripts: $D_{n1}$ (typed slides) vs $D_{1n}$ (ink, notes). This site writes $D_{1n}$.
- Potential difference: $V_{ab}$ means $V_a - V_b$ in the notes and HW, but the slides also write $V_{ab} = -\int_a^b\vec E\cdot d\vec l$, which is $V_b - V_a$. This site always writes the difference out: $V(b) - V(a)$.
- Permittivity of free space: $\epsilon_0$ (slides), $\epsilon_o$ (notes); relative permittivity is given as "$2\epsilon_0$" on the slides and as $\epsilon_r$ in the notes.
- Volume: the notes use $dV$ and $V$ for a volume on the same page where $V$ is the potential; this site writes $d^3\mathbf r'$ or $d\mathcal V$ where it could be confused.
- Unit vectors: $\hat a_x$ (slides) vs $\hat x$ (notes, this site); the spherical radius is $r$ in the notes and $R$ on the review's differential-element table.
- The review deck's "Stoke's theorem" is Stokes' theorem, after G. G. Stokes.

*If you find another one, it belongs here — the pages that quote a slide slip also say so inline, but this is the list to check against before an exam.*

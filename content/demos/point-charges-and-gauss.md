---
title: "Demo — Point charges, field lines, and a movable Gaussian loop"
description: "Drag charges around, watch the field lines rearrange, and move a Gaussian loop to see that the flux through it depends only on the enclosed charge."
tags: [demo, electrostatics, exam-1]
---

*Demo · pairs with [[1-electrostatics/02-coulombs-law-superposition-and-gauss|Lecture 2]] and [[1-electrostatics/03-gauss-law-at-work|Lecture 3]] · concepts: [[concepts/superposition]], [[concepts/flux]], [[concepts/gauss-law]]*

<div class="ece-demo">
<iframe src="/static/demos/point-charges/" title="Point charges and Gauss's law — interactive demo" loading="lazy" style="height:720px"></iframe>
</div>

[Open the demo in its own tab](/static/demos/point-charges/) if the frame is cramped on your screen.

## What to try

1. **Superposition.** Add a second positive charge and drag it toward the first. Watch the field lines between them bend away from each other and a null point appear on the line joining them. Add a negative charge: lines that used to run off to infinity now terminate on it.
2. **Field-line rules.** Lines never cross (except at a null point, where $\mathbf{E}=0$ has no direction); they leave positive charges and end on negative ones; a $+2Q$ charge sprouts twice as many lines as a $+1Q$ charge — this is the "line density ∝ field strength" convention from Lecture 2.
3. **Gauss's law.** Put the loop around a single charge and read the flux: it equals the enclosed charge. Now **drag the loop** so the charge is off-centre, or enlarge it, or move a second charge just *outside* it. The flux does not change. Only what is inside counts, however lopsided the field on the loop looks.
4. **Zero flux, non-zero field.** Put the loop in empty space between two charges. Field lines pass straight through it — in one side, out the other — and the flux reads 0. Compare with the [[problems/flux-through-a-plane-from-two-charges|flux-through-a-plane problem]].
5. **The sign of the normal.** The loop's normals point outward (small ticks). A negative charge inside gives negative flux: lines flow *in* through the loop.

## What the demo is (and isn't)

It is two-dimensional, so the "point" charges are really **line charges seen end-on**: $|\mathbf{E}|\propto 1/r$ rather than $1/r^2$, and Gauss's law reads $\oint\mathbf{E}\cdot\hat{n}\,dl = 2\pi\sum q_{\text{enc}}$ (per unit length, in the demo's units). The *topology* — which lines start and end where, and what the loop encloses — is identical to the 3-D case, which is the point.

The readout integrates $\mathbf{E}\cdot\hat{n}$ numerically around the loop (720 samples of the exact $1/r$ field) and is accurate to the displayed digits unless a charge sits within a pixel or two of the loop itself, where the sampled field is enormous and the sum is unreliable. (The *drawing* uses a slightly softened field near each charge so the arrows stay finite; the readout does not.)

## Related

[[concepts/electric-field]] · [[concepts/coulombs-law]] · [[1-electrostatics/03-gauss-law-at-work#2-gausss-law-restated-for-use|the symmetry ladder in Lecture 3]]

I now have derived a sligthly different understanding:

The dot-product between vectos depends on the vector metric given directly and is defined
only for vectors in the classical approach. What I mixed-up (or conflated as you said) for
the metric for sta4ds is that I tried to put characteristics (spacelike, timelike)
directly into the extended metric. This is simply wrong. When we look for example into the
extended metric of pga3dp the diagonal elements are either one (resulting from the wedge
exomorphism) or zero (splitting between bulk and weight), despite the fact that the
bivectors for example square to -1. So the metric values and the result of the product of
basis blades are completely independent. I understand that now. The the connection between
the metric and the regressive metric is given by the formula "G rG = det(g) I" (G:
extended metric, rG regressive metric (both multiplied with the regular matrix product), g
is the vector metric (defined by the signature of the algebra), det(g) is its determinant,
and I is the identity matrix). This completely describes how extended metric and extended
regressive metric are connected. This formula is valid for every algebra independent of
the metric.

Looking at the EGA and PGA implementation, what Lengyel does is the following:

An extended metric is DEFINED from an exomorphism, extending the given vector metric of
the algebra, which only depends on the signature (e.g. {+,+,+} for ega3d or {-,-,-,+} for
sta4ds). It is derived purely from the wedge product and is defined by G(wdg(a,b)) =
wdg(G(a),G(b)). The antiexomorphism is defined by rG(rwdg(a,b)) = rwdg(rG(a),rG(b)). Both
exomorphisms are completely ignorant of the vector metric as such (they only depend on the
rules for the wedge product or regressive wedge product of the blades), and the resulting
extended metric only reflects the mechanics/rules of the (regressive) wedge product. The
resulting extended metric matrix is then used to DEFINE the dot products and thus the
nrm_sq for the corresponding algebra. This enables dot products from scalars, vectors
(which already existed through the given signature), bivectors, ... to pseudoscalars. Two
scalars always are 1 (s^s==1). The rest is defined by the exomorphism. There is only one
unique solution.

Based on this understanding, I assume that the current implementation of everything that
is metric-related in the current sta4ds implementation is still wrong and needs rework.

Following points need to be covered:

- Does the extended metric calculation in ga_prdxpr reflect the process of generating
  the extended metric from the given vector metric as described above at all? If not, it
  needs to be changed accordingly. We first need to prove that the new version reproduces
  the extended metrics of ega and the (regressive) extended metrics of pga.
- For sta4ds we need to get rid of the is_minkowski part completely. This is just wrong.
  The generation must only depend on the wedge product (or the rules defined for the wedge
  product and the corresponding basis vectors of the algebra). We need to regenerate the
  extended metric and regressive metric. After that is done the tests should show that the
  (anti-)/exomorphism is fulfilled.
- Then we need to regenerate the dot product and the nrm_sq fitting to the extended
  metric.
- Then we need to regenerate all geometric products. They are wrong, since they are
  currently based on the wrong extended metric.

All sufficiently clear or do you still have questions? Can you generate an implementation
plan that reworks ga_prdxpr in terms of extended matrix generation and regeneration of the
geometric products for sta4ds? Then we need to rebuild ga_std4ds_ops*.hpp in a stepwise
procedure to reflect the changes. Please create a stepwise implementation plan.

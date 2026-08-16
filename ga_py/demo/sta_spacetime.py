"""STA in Python — causal character, observer splits, boosts and rotations.

Spacetime algebra is G(1,3,0): the time direction g4 squares to +1, the three
spatial directions to -1. That one signature carries the physics — the sign of a
vector's square IS its causal character, and a Lorentz boost is the same rotor
sandwich that performs a rotation, just with a time-space generator instead of a
purely spatial one.

Run from the project root (see ga_py/README.md §3.3 for PYTHONPATH):

    PYTHONPATH="$PWD/build/ga_py:$PWD/ga_py/python" \\
        ga_py/.venv/bin/python ga_py/demo/sta_spacetime.py

The numbers printed here are pinned in the C++ suite (ga_sta_test); the closed
forms they reproduce (gamma = cosh phi, beta*gamma = sinh phi, additive
rapidity) are stated next to each line.
"""
import math

from ga_py import sta


print("--- 1. the signature is the physics --------------------------------")
u_rest = sta.vec4ds(0.0, 0.0, 0.0, 1.0)  # an observer at rest: 4-velocity g4
x_axis = sta.vec4ds(1.0, 0.0, 0.0, 0.0)  # a purely spatial direction
photon = sta.vec4ds(1.0, 0.0, 0.0, 1.0)  # a light ray: |dx| = |dt|
for name, v in (("u_rest", u_rest), ("x_axis", x_axis), ("photon", photon)):
    print(f"  {name:8s} nrm_sq = {float(sta.nrm_sq(v)):+.1f}   "
          f"timelike={sta.is_timelike(v)!s:5s} spacelike={sta.is_spacelike(v)!s:5s} "
          f"lightlike={sta.is_lightlike(v)}")
print("  (the causal test reads the GEOMETRIC square gr0(X*X), not nrm_sq --- the two")
print("   differ in sign at grades 2 and 3 under the exomorphism metric)")


print("\n--- 2. an observer splits a vector into time and space --------------")
v = sta.vec4ds(1.0, 2.0, 3.0, 4.0)
print("  event v          =", v)
print("  time_split(v, u) =", sta.time_split(v, u_rest), "  (the observer's clock)")
print("  space_split(v,u) =", sta.space_split(v, u_rest),
      "  (the relative 3-vector, as a bivector)")
print("  the split is observer-dependent: it is a contraction with the 4-velocity,")
print("  not a choice of coordinates.")


print("\n--- 3. a boost is a rotor sandwich ---------------------------------")
phi = 0.5  # rapidity
B = sta.get_boost(sta.g14_4ds, phi)  # generator: the time-space plane g14
u = sta.transform(u_rest, B)
print(f"  boost by rapidity phi = {phi} in the g14 (t,x) plane")
print("  boosted 4-velocity =", u)
print(f"    gamma      = {u.w:.13f}   vs cosh(phi) = {math.cosh(phi):.13f}")
print(f"    beta*gamma = {u.x:.13f}   vs sinh(phi) = {math.sinh(phi):.13f}")
print(f"    beta       = {u.x / u.w:.13f}   vs tanh(phi) = {math.tanh(phi):.13f}")
print(f"  rapidity(u_rest, u) = {sta.rapidity(u_rest, u):.13f}  (recovers phi)")
print(f"  nrm_sq stays +1: {float(sta.nrm_sq(u)):.13f}  (a boost preserves the metric)")

# collinear boosts compose by ADDING rapidity --- which is why rapidity, not
# velocity, is the natural parameter (velocities add by the relativistic formula)
B2 = sta.get_boost(sta.g14_4ds, 0.3)
u12 = sta.transform(u_rest, B * B2)
print(f"  0.5 then 0.3 -> rapidity {sta.rapidity(u_rest, u12):.13f} (= 0.8, additive)")
beta1, beta2 = math.tanh(0.5), math.tanh(0.3)
print(f"    the same in velocities: (b1+b2)/(1+b1*b2) = "
      f"{(beta1 + beta2) / (1 + beta1 * beta2):.13f} = tanh(0.8) = {math.tanh(0.8):.13f}")


print("\n--- 4. a rotation is the same sandwich, spatial generator -----------")
theta = math.pi / 2
R = sta.get_rotor(sta.g12_4ds, theta)  # generator: the spatial plane g12
g1, g2 = sta.vec4ds(1.0, 0.0, 0.0, 0.0), sta.vec4ds(0.0, 1.0, 0.0, 0.0)
print("  rotor for theta = pi/2 in the g12 plane")
print("    g1 ->", sta.transform(g1, R))
print("    g2 ->", sta.transform(g2, R))
print("  i.e. g1 -> cos(theta) g1 + sin(theta) g2 --- the same sense ega3d gives for")
print("  e1 in the e12 plane. Both builders here use the POSITIVE half-argument,")
print("  exp(+p/2 B): ega3d needs an explicit -theta/2 for that sense, STA does not,")
print("  because the spatial vectors squaring to -1 already supply the flip.")
print("  A rotation leaves the causal character alone:",
      sta.is_spacelike(sta.transform(g1, R)))

# boost and rotation compose by multiplying the rotors --- one group, one sandwich
M = R * B
print("  R * B applied to u_rest =", sta.transform(u_rest, M),
      "(boost, then rotate)")


print("\n--- 5. exp / log / sqrt and the double cover ------------------------")
print("  log(B)      =", sta.log(B), " --- the generator carries phi/2 = 0.25")
print("  exp(log(B)) == B :", sta.is_same_transform(sta.exp(sta.log(B)), B, 1e-12))
half = sta.sqrt(B)
print("  sqrt(B)     =", half)
print("  sqrt(B)^2   == B :", sta.is_same_transform(half * half, B, 1e-12))
print(f"    half the rapidity: {sta.rapidity(u_rest, sta.transform(u_rest, half)):.13f}"
      f"  (= {phi / 2})")
print("  is_same_transform(B, -B) :", sta.is_same_transform(B, -B, 1e-12),
      "(rotors double cover the transformations)")
print("  B == -B component-wise   :", B == -B,
      "--- same motion, different representation")

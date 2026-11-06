import math
from scipy.optimize import minimize

T = 4955
L = 1000e-3
shearModulus = 95e9
maxTwist = math.radians(1.5)
maxStress = 140e6
rHole = 15e-3


def stress(radii):
    ri = radii[0]
    ro = radii[1]
    j = J(ri, ro)
    # holeFactor gives the peak normal stress at the hole edge; the edge is in
    # uniaxial stress, so the maximum shear stress is half of it
    return 0.5 * T * ro * holeFactor(ri, ro) / j


def holeFactor(ri, ro):
    c1 = 4
    c2 = -6.055 + 3.184 * radiusRatio(ri, ro) - 3.461 * radiusRatio(ri, ro) ** 2
    c3 = 32.764 - 30.121 * radiusRatio(ri, ro) + 39.887 * radiusRatio(ri, ro) ** 2
    c4 = (
        -38.330 + 51.542 * math.sqrt(radiusRatio(ri, ro)) - 27.483 * radiusRatio(ri, ro)
    )
    k = c1 + c2 * holeRatio(ro) + c3 * holeRatio(ro) ** 2 + c4 * holeRatio(ro) ** 3
    return k


def holeRatio(ro):
    return rHole / ro


def radiusRatio(ri, ro):
    return ri / (ro)


def angleOfTwist(radii):
    return T * L / (J(radii[0], radii[1]) * shearModulus)


def J(ri, ro):
    return math.pi / 2 * (ro**4 - ri**4)


# objective function
def crossSection(radii):
    return math.pi * (radii[1] ** 2 - radii[0] ** 2)


constraints = [
    {"type": "ineq", "fun": lambda radii: radii[1] - radii[0]},
    {"type": "ineq", "fun": lambda radii: radii[0] - 1e-6},
    {"type": "ineq", "fun": lambda radii: radii[1] - 1e-6},
    {"type": "ineq", "fun": lambda radii: radii[1] - radii[0] - 5e-3},
    {"type": "ineq", "fun": lambda radii: maxTwist - angleOfTwist(radii)},
    {"type": "ineq", "fun": lambda radii: maxStress - stress(radii)},
    # validity range of the holeFactor, not design requirements
    {"type": "ineq", "fun": lambda radii: 0.8 - radii[0] / radii[1]},
    {"type": "ineq", "fun": lambda radii: 0.4 - holeRatio(radii[1])},
]

res = minimize(
    crossSection,
    [70e-3, 100e-3],
    method="SLSQP",
    constraints=constraints,
)

ri, ro = res.x
print(f"Converged:      {res.success} ({res.message})")
print(f"Inner radius:   {ri * 1e3:.2f} mm")
print(f"Outer radius:   {ro * 1e3:.2f} mm")
print(f"Wall thickness: {(ro - ri) * 1e3:.2f} mm (min 5 mm)")
print(f"Cross-section:  {crossSection(res.x) * 1e6:.0f} mm^2")
print(f"Twist:          {math.degrees(angleOfTwist(res.x)):.3f} deg (max 1.5 deg)")
print(f"Shear stress:   {stress(res.x) / 1e6:.1f} MPa (max 140 MPa)")

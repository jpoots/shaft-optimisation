# Hollow Shaft Optimisation

Minimises the cross-sectional area (and so the mass) of a hollow steel transmission shaft under torsion, using SciPy's SLSQP optimiser. The shaft has a 30 mm transverse hole through its full diameter, so the stress check includes the stress concentration at the hole.

The problem is based on a university CAE coursework brief.

## The problem

| Item | Value |
|---|---|
| Torque, T | 4,955 N·m |
| Length, L | 1 m |
| Transverse hole | 30 mm diameter, through the full shaft |
| Material | S235 steel: shear modulus G = 95 GPa, shear yield strength 140 MPa |
| Original design | Outer radius 100 mm, inner radius 70 mm |

Design variables: inner radius r<sub>i</sub> and outer radius r<sub>o</sub>.

Objective: minimise the cross-sectional area, π(r<sub>o</sub>² − r<sub>i</sub>²). The length is fixed, so this also minimises volume and mass.

Constraints:

| Constraint | Expression | Source |
|---|---|---|
| Shear stress at the hole | τ<sub>max</sub> ≤ 140 MPa | Brief |
| Angle of twist | θ = TL / (GJ) ≤ 1.5° | Brief |
| Wall thickness | r<sub>o</sub> − r<sub>i</sub> ≥ 5 mm | Brief |
| Wall ratio | r<sub>i</sub> / r<sub>o</sub> ≤ 0.8 | Validity range of the stress concentration fit |
| Hole ratio | hole diameter / outer diameter ≤ 0.4 | Validity range of the stress concentration fit |

where J = π/2 (r<sub>o</sub>⁴ − r<sub>i</sub>⁴) is the polar second moment of area.

## Stress at the hole

The nominal torsional shear stress is T r<sub>o</sub> / J. The hole raises it locally, by a stress concentration factor K from a published curve fit for a tube with a transverse hole in torsion:

K = C<sub>1</sub> + C<sub>2</sub>(d/D) + C<sub>3</sub>(d/D)² + C<sub>4</sub>(d/D)³

where d/D is the hole diameter over the outer diameter, and C<sub>2</sub> to C<sub>4</sub> depend on the inner over outer diameter. K is the ratio of the peak normal stress at the edge of the hole to the nominal shear stress. The edge of the hole is in uniaxial stress, so the maximum shear stress there is half the peak normal stress:

τ<sub>max</sub> = 0.5 K T r<sub>o</sub> / J

The fit is only valid for r<sub>i</sub> / r<sub>o</sub> ≤ 0.8 and d/D ≤ 0.4, so these are included as constraints. They are limits of the data, not design requirements.

## Results

| | Original design | SLSQP from the original design |
|---|---|---|
| Inner radius | 70 mm | 0.12 mm |
| Outer radius | 100 mm | 37.50 mm |
| Cross-section | 16,022 mm² | 4,418 mm² |
| Mass | 125.8 kg | 34.7 kg |
| Twist | | 0.96° |
| Shear stress at the hole | | 135.7 MPa |
| Binding constraint | | Hole ratio (d/D = 0.4) |

Masses use a density of 7,850 kg/m³ and ignore the small volume of the hole.

**SLSQP is a local, gradient-based optimiser.** From the original design it stops at a non-optimal point while reporting success. A multi-start (running from several starting points and keeping the best) or a global method would avoid this.

## Running it

```
pip install scipy
python3 shaft_optimisation.py
```

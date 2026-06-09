"""
Error Budget — force contribution hierarchy at the April 2029 close approach.

Method: a BASELINE model is compared against models with ONE extra force
activated at a time. The 3D position difference (meters) vs BASELINE measures
each force's contribution.

BASELINE: Sun + 8 planets + Moon + Pluto (ASSIST "PLANETS" covers Moon+Pluto)
          + Earth harmonics (J2,J3,J4) + Sun harmonics (J2) + GR_SIMPLE.

  TEST 1  Yarkovsky : BASELINE + NON_GRAVITATIONAL (A1=5.0e-13, A2=-2.9e-14)
  TEST 2  Asteroids : BASELINE + 16 massive asteroids (sb441-n16.bsp)
  TEST 3  GR_EIH    : BASELINE with GR_SIMPLE replaced by GR_EIH (11 bodies)
  TEST 4  FULL      : all of the above (= "ASSIST ONA", Descargas/Apophis.ipynb)

Evaluation epoch: April close approach, t = 10694.5 days (2462239.5 - jd_ref).
"""
import pathlib
import numpy as np
import rebound
import assist

DATA = pathlib.Path(__file__).parent.parent / "data"
AU2M = 149_597_870_700.0          # meters per AU

ephem  = assist.Ephem(str(DATA / "de440.bsp"), str(DATA / "sb441-n16.bsp"))
JD_REF = ephem.jd_ref             # 2451545.0

t_initial = 2.4621385359989386e6 - JD_REF
t_ca      = 10694.5               # April 13 close approach epoch

# Yarkovsky parameters (A1, A2, A3)
A1, A2, A3 = 5.0e-13, -2.9e-14, 0.0

# Apophis initial conditions (Holman et al. 2023), heliocentric -> barycentric
apophis_helio = rebound.Particle(
    x=-5.5946538550488512e-1, y= 8.5647564757574512e-1, z= 3.0415066217102493e-1,
    vx=-1.3818324735921638e-2, vy=-6.0088275597939191e-3, vz=-2.5805044631309632e-3)


def run(forces, yark=False, eih=False):
    """Integrate Apophis to t_ca with the given ASSIST force set; return xyz (AU)."""
    sim = rebound.Simulation()
    sim.add(apophis_helio + ephem.get_particle("sun", t_initial))
    sim.t = t_initial
    sim.ri_ias15.min_dt = 1e-4
    ex = assist.Extras(sim, ephem)
    ex.forces = forces
    if eih:
        ex.gr_eih_sources = 11
    if yark:
        ex.particle_params = np.array([A1, A2, A3])
    ex.integrate_or_interpolate(t_ca)
    return np.array(sim.particles[0].xyz)


BASE = ["SUN", "PLANETS", "EARTH_HARMONICS", "SUN_HARMONICS", "GR_SIMPLE"]

pos_base = run(BASE)
pos_t1   = run(BASE + ["NON_GRAVITATIONAL"], yark=True)
pos_t2   = run(BASE + ["ASTEROIDS"])
pos_t3   = run(["SUN", "PLANETS", "EARTH_HARMONICS", "SUN_HARMONICS", "GR_EIH"], eih=True)
pos_t4   = run(["SUN", "PLANETS", "ASTEROIDS", "NON_GRAVITATIONAL",
                "EARTH_HARMONICS", "SUN_HARMONICS", "GR_EIH"], yark=True, eih=True)

def dm(p):
    return np.linalg.norm(p - pos_base) * AU2M   # meters

rows = [
    ("TEST 1  Yarkovsky (A1,A2)",        dm(pos_t1)),
    ("TEST 2  Asteroids (16 bodies)",    dm(pos_t2)),
    ("TEST 3  GR_EIH (vs GR_SIMPLE)",    dm(pos_t3)),
    ("TEST 4  FULL (ASSIST ONA)",        dm(pos_t4)),
]
linear_sum = rows[0][1] + rows[1][1] + rows[2][1]

print("\n" + "=" * 60)
print(f"ERROR BUDGET — April close approach (t = {t_ca} days)")
print("BASELINE: Sun+Planets+Moon+Pluto + Earth/Sun harmonics + GR_SIMPLE")
print("=" * 60)
print(f"{'Force activated':<34}{'Δr vs BASELINE [m]':>24}")
print("-" * 60)
for name, d in rows:
    print(f"{name:<34}{d:>24.2f}")
print("-" * 60)
print(f"{'Σ (TEST 1+2+3, linear)':<34}{linear_sum:>24.2f}")
print(f"{'TEST 4 (FULL, nonlinear)':<34}{rows[3][1]:>24.2f}")
print("=" * 60)
print("Note: TEST4 != linear sum -> forces couple nonlinearly near the encounter.")

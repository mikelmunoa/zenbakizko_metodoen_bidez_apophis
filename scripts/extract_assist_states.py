"""
ASSIST eredu osoaren EGOERAK (posizioa + abiadura) atera hainbat unetan:
  - hurbilketa-unea (perigeoa) aurkitu
  - perigeo - 1 egun  (hurbilketaren AURRETIK)
  - perigeo + 1 egun  (hurbilketaren OSTEAN)
  - 2030-01-01        (urte amaierako erreferentzia)

Irteera: egoerak AU, AU/egun-tan (barizentriko), GrAL-en hasteko.
"""
import pathlib
import numpy as np
import rebound
import assist

DATA = pathlib.Path(__file__).parent.parent / "data"
ephem = assist.Ephem(str(DATA / "de440.bsp"), str(DATA / "sb441-n16.bsp"))
JD_REF = ephem.jd_ref

t_initial = 2.4621385359989386e6 - JD_REF
t_final   = 2.4625030372426095e6 - JD_REF   # 2030-01-01

apophis = rebound.Particle(
    x=-5.5946538550488512e-1, y= 8.5647564757574512e-1, z= 3.0415066217102493e-1,
    vx=-1.3818324735921638e-2, vy=-6.0088275597939191e-3, vz=-2.5805044631309632e-3,
) + ephem.get_particle("sun", t_initial)

def fresh_sim():
    s = rebound.Simulation(); s.add(apophis); s.t = t_initial; s.ri_ias15.min_dt = 1e-4
    e = assist.Extras(s, ephem); e.gr_eih_sources = 11
    e.particle_params = np.array([5.0e-13, -2.9e-14, 0.0])
    return s, e

# ── 1. Perigeoa aurkitu (geozentriko distantzia minimoa) ─────────────────────
# OHARRA: t ABSOLUTUA da (jd_ref-etik egun). Hurbilketa 2029-04-13 ≈ t=10694.5.
sim, ex = fresh_sim()
best_t, best_d = None, 1e9
for t in np.arange(10693.0, 10697.0, 0.002):
    ex.integrate_or_interpolate(float(t))
    ap = np.array(sim.particles[0].xyz)
    ea = np.array(ephem.get_particle("earth", float(t)).xyz)
    d = np.linalg.norm(ap - ea)
    if d < best_d:
        best_d, best_t = d, float(t)
AU2KM = 149_597_870.7
print(f"# Perigeoa: t = {best_t:.4f} egun, geozentriko distantzia = {best_d*AU2KM:.0f} km")
print(f"#           (2029-04-13 ingurua, ~{best_d*AU2KM/6371:.1f} Lur-erradio)")

t_before = round(best_t - 1.0, 4)
t_after  = round(best_t + 1.0, 4)

# ── 2. Egoerak atera (gorakako ordenan) ──────────────────────────────────────
sim, ex = fresh_sim()
print("\n# Egoerak (AU, AU/egun) — barizentriko J2000:")
for label, t in [("AURRETIK (peri-1d)", t_before),
                 ("PERIGEOA",           round(best_t, 4)),
                 ("OSTEAN   (peri+1d)",  t_after),
                 ("2030-01-01",          t_final)]:
    ex.integrate_or_interpolate(float(t))
    p = sim.particles[0]
    print(f"{label:20s} t={t:.6f}  "
          f"r=[{p.x:.15e}, {p.y:.15e}, {p.z:.15e}]  "
          f"v=[{p.vx:.15e}, {p.vy:.15e}, {p.vz:.15e}]")

"""
apophis_sinplea_urte_osoa.txt sortu zuzen — simulazio freskoa.

Jatorria: assist/jupyter_examples/Apophis.ipynb
Indarrak: SUN + PLANETS + EARTH_HARMONICS + SUN_HARMONICS + GR_SIMPLE
Yarkovsky: Ez
t_final:   2.4625030372426095e6 - jd_ref  (2030-01-01)

Konponketa: Cell 13 → Cell 14 exekuzio-ordena arazoa saihesten da simulazio
freskoa erabiliz (integrate_or_interpolate deitu aurretik sim.t = t_initial).
"""
import pathlib
import numpy as np
import rebound
import assist

DATA = pathlib.Path("/home/mikel/Ingenieritza_Informatikoa/GrAL/GrAL_github/data")
OUT  = DATA / "apophis_sinplea_urte_osoa.txt"

ephem    = assist.Ephem(str(DATA / "de440.bsp"), str(DATA / "sb441-n16.bsp"))
JD_REF   = ephem.jd_ref

t_initial = 2.4621385359989386e6 - JD_REF
t_final   = 2.4625030372426095e6 - JD_REF   # 2030-01-01

apophis_helio = rebound.Particle(
    x=-5.5946538550488512e-1, y= 8.5647564757574512e-1, z= 3.0415066217102493e-1,
    vx=-1.3818324735921638e-2, vy=-6.0088275597939191e-3, vz=-2.5805044631309632e-3,
)
apophis_bary = apophis_helio + ephem.get_particle("sun", t_initial)

N_samples = 10_000
times = np.linspace(t_initial, t_final, N_samples, endpoint=True)

sim = rebound.Simulation()
sim.add(apophis_bary)
sim.t = t_initial
sim.ri_ias15.min_dt = 0.001

extras = assist.Extras(sim, ephem)
extras.forces = ["SUN", "PLANETS", "EARTH_HARMONICS", "SUN_HARMONICS", "GR_SIMPLE"]

print(f"Indarrak: {extras.forces}")
print(f"t_initial = {t_initial:.10f} egun")
print(f"t_final   = {t_final:.10f} egun  (2030-01-01)")
print(f"Iraupena  = {t_final - t_initial:.4f} egun")
print(f"Puntu kop.: {N_samples}")
print(f"Irteerako fitxategia: {OUT}")
print()

with open(OUT, "w") as f:
    for i, t in enumerate(times):
        extras.integrate_or_interpolate(t)
        x, y, z = sim.particles[0].xyz
        f.write(f"{t:.15e}[{x:.16e}  {y:.16e}  {z:.16e}]\n")
        if i % 1000 == 0:
            print(f"  {i}/{N_samples}  t={t:.4f}")

print(f"\nEginda — {N_samples} lerro idatzia: {OUT}")

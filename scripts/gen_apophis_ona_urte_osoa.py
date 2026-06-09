"""
apophis_ona_urte_osoa.txt sortu — ASSIST eredu OSOA (ONA), urte osorako (2030-01-01).

Indarrak: SUN+PLANETS+ASTEROIDS+EARTH_HARMONICS+SUN_HARMONICS+GR_EIH(11) + Yarkovsky.
(= Descargas/Apophis.ipynb-eko 'ASSIST ONA', baina t_final=2030.)

Beharrezkoa da: apophis_urte_osoa.txt gaizki izendatuta dago (sinplifikatua da),
beraz hau da urte osorako EGIAZKO eredu osoa, GrAL+Yarkovsky alderatzeko.

Yarkovsky: A1=5.0e-13, A2=-2.9e-14, A3=0  (GrAL-en physics tuplearekin bat).
Formatua: <t_egun>[x_AU  y_AU  z_AU]  (Julia split() bateragarria).
"""
import pathlib
import numpy as np
import rebound
import assist

DATA = pathlib.Path(__file__).parent.parent / "data"
OUT  = DATA / "apophis_ona_urte_osoa.txt"

ephem   = assist.Ephem(str(DATA / "de440.bsp"), str(DATA / "sb441-n16.bsp"))
JD_REF  = ephem.jd_ref

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
sim.ri_ias15.min_dt = 1e-4

extras = assist.Extras(sim, ephem)
# Ez da extras.forces ezartzen -> indar guztiak (eredu osoa)
extras.gr_eih_sources = 11
extras.particle_params = np.array([5.0e-13, -2.9e-14, 0.0])   # Yarkovsky (GrAL-ekin bat)

print(f"Indarrak (lehenetsiak): {extras.forces}")
print(f"GR_EIH iturriak: {extras.gr_eih_sources}")
print(f"Yarkovsky A1=5.0e-13  A2=-2.9e-14")
print(f"t_final = {t_final:.10f} egun (2030-01-01)")
print(f"Irteera: {OUT}")
print()

with open(OUT, "w") as f:
    for i, t in enumerate(times):
        extras.integrate_or_interpolate(t)
        x, y, z = sim.particles[0].xyz
        f.write(f"{t:.15e}[{x:.16e}  {y:.16e}  {z:.16e}]\n")
        if i % 1000 == 0:
            print(f"  {i}/{N_samples}  t={t:.4f}")

print(f"\nEginda — {N_samples} lerro: {OUT}")

"""
apophis_pos_april.txt sortu — Descargas/Apophis.ipynb konfigurazioa.

Indarrak (lehenetsiak, ez da extras.forces ezartzen):
    SUN + PLANETS + ASTEROIDS + NON_GRAVITATIONAL + EARTH_HARMONICS + SUN_HARMONICS + GR_EIH
GR_EIH iturriak: 11 (Eguzki + 10 planeta)
Yarkovsky:       A1=4.999999873689e-13, A2=-2.901085508711e-14, A3=0
min_dt:          0.001 egun

Irteerako formatua: <t_egun>[x_AU  y_AU  z_AU]  (espazio bidez bereizita — Julia split() bateragarria)
"""
import pathlib
import numpy as np
import rebound
import assist

DATA = pathlib.Path(__file__).parent.parent / "data"
OUT  = DATA / "apophis_pos_april.txt"

ephem   = assist.Ephem(str(DATA / "de440.bsp"), str(DATA / "sb441-n16.bsp"))
JD_REF  = ephem.jd_ref   # 2451545.0

t_initial = 2.4621385359989386e6 - JD_REF
t_final   = 2462239.50000        - JD_REF

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
sim.ri_ias15.min_dt = 0.001   # Descargas/Apophis.ipynb-ko berdina

extras = assist.Extras(sim, ephem)
# Ez da extras.forces ezartzen → indar guztiak lehenetsiak:
#   SUN + PLANETS + ASTEROIDS + NON_GRAVITATIONAL + EARTH_HARMONICS + SUN_HARMONICS + GR_EIH
extras.gr_eih_sources = 11                                           # Eguzki + 10 planeta
extras.particle_params = np.array([4.999999873689e-13,              # Yarkovsky A1
                                   -2.901085508711e-14,              # Yarkovsky A2
                                   0.0])                             # A3=0

print(f"t_initial = {t_initial:.15f} egun")
print(f"t_final   = {t_final:.15f} egun")
print(f"Iraupena  = {t_final - t_initial:.6f} egun")
print(f"Puntu kop.: {N_samples}")
print(f"Indarrak (lehenetsiak): {extras.forces}")
print(f"GR_EIH iturriak: {extras.gr_eih_sources}")
print(f"Yarkovsky A1={extras._particle_params[0]:.6e}  A2={extras._particle_params[1]:.6e}")
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

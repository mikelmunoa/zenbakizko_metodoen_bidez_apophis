"""
ASSIST GR_SIMPLE posizio-fitxategia sortu.
Irteerako formatua: apophis_pos_april.txt-ren berdina → <t_egun>[x_AU  y_AU  z_AU]

Hiru konfigurazio:
  1. GR_EIH  + ASTEROIDS (erreferentzia — apophis_pos_april.txt-ren berdina)
  2. GR_SIMPLE            (D-D-ren baliokidea ASSIST-en arabera)
  3. GR_EIH, iturri=1    (EIH Eguzki bakarrik = Schwarzschild)
"""
import sys, pathlib
import numpy as np
import rebound
import assist

ASSIST_DATA = pathlib.Path("/home/mikel/Documentos/Ingenieritza_Informatikoa/GrAL/assist/data")
APRIL_FILE  = pathlib.Path("/home/mikel/Ingenieritza_Informatikoa/GrAL/GrAL_github/data/apophis_pos_april.txt")
OUT_DIR     = pathlib.Path("/home/mikel/Ingenieritza_Informatikoa/GrAL/GrAL_github/data")

# ── Denbora-puntuak april fitxategitik irakurri ───────────────────────────────
times_days = []
with open(APRIL_FILE) as f:
    for line in f:
        k = line.index('[')
        times_days.append(float(line[:k]))
times_days = np.array(times_days)
t_ini = times_days[0]
t_end = times_days[-1]
print(f"Denbora-tartea: {t_ini:.4f} → {t_end:.4f} egun ({len(times_days)} puntu)")

# ── Hasierako baldintzak ──────────────────────────────────────────────────────
ephem = assist.Ephem(str(ASSIST_DATA/"de440.bsp"), str(ASSIST_DATA/"sb441-n16.bsp"))
JD_REF = ephem.jd_ref   # J2000 = 2451545.0

apophis_helio = rebound.Particle(
    x=-5.594653855048851e-1, y= 8.564756475757451e-1, z= 3.041506621710249e-1,
    vx=-1.381832473592163e-2, vy=-6.008827559793919e-3, vz=-2.580504463130963e-3
)
sun_ini = ephem.get_particle("sun", t_ini)
apophis_bary = apophis_helio + sun_ini

def run_assist(forces_list, gr_eih_src=11, label=""):
    """ASSIST exekutatu eta posizio-array bat itzuli (n×3, AU)."""
    sim = rebound.Simulation()
    sim.add(apophis_bary)
    sim.t = t_ini
    sim.ri_ias15.min_dt = 1e-4   # egun — zehaztasuna

    ext = assist.Extras(sim, ephem)
    ext.forces = forces_list
    if "GR_EIH" in forces_list:
        ext.gr_eih_sources = gr_eih_src

    pos = np.zeros((len(times_days), 3))
    for i, t in enumerate(times_days):
        ext.integrate_or_interpolate(t)
        p = sim.particles[0]
        pos[i] = [p.x, p.y, p.z]
        if i % 1000 == 0:
            print(f"  {label}  {i}/{len(times_days)}  t={t:.2f}")
    return pos

# ── 1. GR_EIH + ASTEROIDS (erreferentzia berdina) ────────────────────────────
print("\n[1/3] GR_EIH (11 iturri) + ASTEROIDS ...")
pos_eih = run_assist(["SUN","PLANETS","ASTEROIDS","GR_EIH"], gr_eih_src=11, label="EIH+AST")

# ── 2. GR_SIMPLE (D-D baliokidea) ────────────────────────────────────────────
print("\n[2/3] GR_SIMPLE (Eguzki bakarrik, Schwarzschild) ...")
pos_simple = run_assist(["SUN","PLANETS","GR_SIMPLE"], label="GR_SIMPLE")

# ── 3. GR_EIH iturri=1 (Schwarzschild EIH bidez) ─────────────────────────────
print("\n[3/3] GR_EIH (1 iturri = Schwarzschild EIH) ...")
pos_eih1 = run_assist(["SUN","PLANETS","GR_EIH"], gr_eih_src=1, label="EIH1")

# ── Fitxategiak idatzi ────────────────────────────────────────────────────────
def write_pos(fname, times, pos):
    with open(fname, "w") as f:
        for t, (x, y, z) in zip(times, pos):
            f.write(f"{t:.15e}[{x:.16e}  {y:.16e}  {z:.16e}]\n")
    print(f"Idatzia: {fname}  ({len(times)} lerro)")

write_pos(OUT_DIR/"apophis_assist_gr_eih.txt",    times_days, pos_eih)
write_pos(OUT_DIR/"apophis_assist_gr_simple.txt", times_days, pos_simple)
write_pos(OUT_DIR/"apophis_assist_gr_eih1.txt",   times_days, pos_eih1)

# ── Azken puntuko diferentziak ────────────────────────────────────────────────
AU2KM = 149597870.7
d_simple_eih = np.linalg.norm(pos_simple[-1] - pos_eih[-1]) * AU2KM * 1000
d_eih1_eih   = np.linalg.norm(pos_eih1[-1]  - pos_eih[-1]) * AU2KM * 1000
d_simple_eih1= np.linalg.norm(pos_simple[-1] - pos_eih1[-1]) * AU2KM * 1000

print(f"\n── Azken puntuko diferentziak (hurbilpen-unea) ──────────────────")
print(f"  GR_SIMPLE vs GR_EIH(11)   : {d_simple_eih:.2f} m")
print(f"  GR_EIH(1) vs GR_EIH(11)   : {d_eih1_eih:.2f} m")
print(f"  GR_SIMPLE vs GR_EIH(1)    : {d_simple_eih1:.4f} m")
print("\nEginda.")

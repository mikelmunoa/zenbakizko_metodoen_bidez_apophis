# ──────────────────────────────────────────────────────────────────────────────
# Sistema lehenetsia: gorputzak eta GM-ak ID-etatik (SPICE bidez), erabiltzaileak
# GM balio bakoitza eskuz idatzi gabe (ASSIST-en "PLANETS" parametroaren antzera).
# GM-ak gm_de440.tpc-tik hartzen dira (bodvrd) — DE440 goiburuko balio autoritatiboak.
# ──────────────────────────────────────────────────────────────────────────────

"Eguzkia + 8 planeta (barizentriko) + Ilargia + Pluton (NAIF ID-ak)."
const PLANET_IDS = [10, 399, 301, 1, 2, 4, 5, 6, 7, 8, 9]

"sb441-n16 asteroide-multzoa (16 nagusiak), ASSIST-en berdina (NAIF ID-ak)."
const ASTEROID_IDS = [2_000_001, 2_000_002, 2_000_003, 2_000_004, 2_000_007, 2_000_010,
                      2_000_015, 2_000_016, 2_000_031, 2_000_052, 2_000_065, 2_000_087,
                      2_000_088, 2_000_107, 2_000_511, 2_000_704]

"sb441-n16 asteroideen GM-ak (km³/s²), ASTEROID_IDS ordena berean (ASSIST/sb441 balioak)."
const ASTEROID_GMS = [62.6289, 13.6659, 1.9206, 17.2882, 1.1399, 5.6251, 2.0230, 1.5897,
                      1.0794, 2.6830, 0.9381, 2.1682, 1.1898, 1.4437, 3.8945, 2.8304]

"""
Eredu fisiko lehenetsia: GR_EIH (11 iturri) + Eguzki J2 + Lur J2/J3/J4.
Yarkovsky objektuaren araberakoa denez, ez dago hemen; gehitzeko NamedTuple-a bateratu:
`(; DEFAULT_PHYSICS..., yarkovsky=true, A1=..., A2=..., A3=...)`.
"""
const DEFAULT_PHYSICS = (gr=false, gr_eih=true, n_eih=length(PLANET_IDS),
                         j2_sun=true, j2_earth=true, j3_earth=true, j4_earth=true)

"""
    planet_system(tiv; datadir="data", asteroids=false) -> (ids, mus, bodies)

Sistema osoa eraiki ID-etatik, GM balioak eskuz idatzi gabe. GM-ak `gm_de440.tpc`-tik
hartzen dira (`bodvrd`, DE440 goiburuko balioak), eta efemeride-koefizienteak `tiv`
denbora-tartean sortzen dira. `asteroids=true` bada, sb441-n16 16 asteroideak gehitzen ditu.

SPICE kernelak aldez aurretik kargatuta egon behar dira (`naif0012.tls`, `de440.bsp`,
eta `asteroids=true` bada `sb441-n16.bsp`); `gm_de440.tpc` funtzioak berak kargatzen du.

# Adibidea
```julia
furnsh("data/naif0012.tls", "data/de440.bsp", "data/sb441-n16.bsp")
ids, mus, bodies = planet_system((et_0, et_end); asteroids=true)
p = (mus=mus, ids=ids, bodies=bodies,
     physics=(; DEFAULT_PHYSICS..., yarkovsky=true, A1=5.0e-13, A2=-2.9e-14, A3=0.0))
```
"""
function planet_system(tiv::Tuple{<:Real,<:Real}; datadir::AbstractString="data",
                       asteroids::Bool=false,
                       coeffs_json::AbstractString=joinpath(datadir, "coeffs_system.json"),
                       coeffs_csv::AbstractString=joinpath(datadir, "coeffs_system.csv"))
    furnsh(joinpath(datadir, "gm_de440.tpc"))          # GM-ak ID-etik (bodvrd)
    ids = asteroids ? vcat(PLANET_IDS, ASTEROID_IDS) : copy(PLANET_IDS)
    mus = Float64[bodvrd(string(i), "GM", 1)[1] for i in PLANET_IDS]
    asteroids && append!(mus, ASTEROID_GMS)
    fallback = Dict{Int,Tuple{Int,Int}}(399 => (13, 8))   # Lurra (399) EMB-tik eratorria
    if asteroids
        for id in ASTEROID_IDS
            fallback[id] = (13, 8)
        end
    end
    create_coeffs_file(coeffs_json, coeffs_csv, ids, fill(tiv, length(ids)); fallback_params=fallback)
    bodies = [BodyCoeffs(coeffs_json, coeffs_csv, id, tiv) for id in ids]
    return ids, mus, bodies
end

#=
GrAL-en posizioa bi unetan kalkulatu eta idatzi.
  - Gerturapena: 2029-04-13
  - 2030 hasiera

UNITATEAK: GrAL barnean km, km/s; emaitza AU-tara (/ AU2KM).
DENBORAK:  t_egun → et (× DAY2S). t_final erabili (EZ ASSIST-en sim.t).
=#

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GravitationSimulation, LittleEphemeris
using SciMLBase, IRKGaussLegendre
using SPICE, Printf

const AU2KM = 149_597_870.7   # km/AU
const DAY2S = 86_400.0        # s/egun

DATA = joinpath(@__DIR__, "..", "data")
furnsh(joinpath(DATA,"naif0012.tls"), joinpath(DATA,"de440.bsp"), joinpath(DATA,"sb441-n16.bsp"))

# ── Denborak ─────────────────────────────────────────────────────────────────
const JD_REF = 2_451_545.0
t_initial = 2.4621385359989386e6 - JD_REF
const et_0 = t_initial * DAY2S

t_ger = 10694.5                # gerturapena   = 2462239.5            - jd_ref
t_30  = 10958.0372426095       # 2030 hasiera  = 2462503.0372426095  - jd_ref

# ── Hasiera-baldintzak (Holman et al. 2023): AU, AU/egun → km, km/s ──────────
helio_pos = [-5.5946538550488512e-1,  8.5647564757574512e-1,  3.0415066217102493e-1]
helio_vel = [-1.3818324735921638e-2, -6.0088275597939191e-3, -2.5805044631309632e-3]
sun_st, _ = spkez(Int32(10), et_0, "J2000", "NONE", Int32(0))
u0 = [helio_pos[1]*AU2KM + sun_st[1],        helio_pos[2]*AU2KM + sun_st[2],        helio_pos[3]*AU2KM + sun_st[3],
      helio_vel[1]*AU2KM/DAY2S + sun_st[4],  helio_vel[2]*AU2KM/DAY2S + sun_st[5],  helio_vel[3]*AU2KM/DAY2S + sun_st[6]]

# ── Gorputzak, GM-ak eta indarrak ───────────────────────────────────────────
body_ids = [10, 399, 301, 1, 2, 4, 5, 6, 7, 8]
body_mus = [1.32712440041279e11, 398600.436, 4902.800, 22031.869, 324858.592,
            42828.376, 1.267127641e8, 3.7940584842e7, 5.7945564e6, 6.836527101e6]
ast_ids  = [2_000_001,2_000_002,2_000_003,2_000_004,2_000_007,2_000_010,
            2_000_015,2_000_016,2_000_031,2_000_052,2_000_065,2_000_087,
            2_000_088,2_000_107,2_000_511,2_000_704]
ast_mus  = [62.6289,13.6659,1.9206,17.2882,1.1399,5.6251,2.0230,1.5897,
            1.0794,2.6830,0.9381,2.1682,1.1898,1.4437,3.8945,2.8304]

tiv = (et_0 - 5*DAY2S, t_30*DAY2S + 5*DAY2S)
all_ids = vcat(body_ids, ast_ids)
fallback = Dict{Int,Tuple{Int,Int}}(id => (13,8) for id in ast_ids); fallback[399] = (13,8)
create_coeffs_file(joinpath(DATA,"coeffs_bi.json"), joinpath(DATA,"coeffs_bi.csv"),
                   all_ids, fill(tiv, length(all_ids)); fallback_params=fallback)
bodies = [BodyCoeffs(joinpath(DATA,"coeffs_bi.json"), joinpath(DATA,"coeffs_bi.csv"), id, tiv) for id in all_ids]
physics = (gr=true, j2_sun=true, j2_earth=true, j3_earth=true, j4_earth=true,
           yarkovsky=true, A1=5.0e-13, A2=-2.9e-14, A3=0.0)
p = (mus=vcat(body_mus,ast_mus), ids=all_ids, bodies=bodies, physics=physics)

# ── Posizioa kalkulatu (et_0 → t_egun), AU-tan ───────────────────────────────
solve(ODEProblem(f_master!, u0, (et_0, et_0+3600.0), p), IRKGL16();
      reltol=1e-12, abstol=1e-12, save_everystep=false)   # JIT warm-up

function gral_pos(t_egun)
    et = t_egun * DAY2S
    sol = solve(ODEProblem(f_master!, u0, (et_0, et), p), IRKGL16();
                reltol=1e-12, abstol=1e-12, saveat=[et])
    sol.u[end][1:3] ./ AU2KM
end

pos_ger = gral_pos(t_ger)
pos_30  = gral_pos(t_30)

@printf("Gerturapena   (t=%.10f egun): [%.15f, %.15f, %.15f]\n", t_ger, pos_ger...)
@printf("2030 hasiera  (t=%.10f egun): [%.15f, %.15f, %.15f]\n", t_30,  pos_30...)

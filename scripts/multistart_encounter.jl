using GravitationSimulation, LittleEphemeris, SciMLBase, IRKGaussLegendre
using SPICE, LinearAlgebra, Printf
const AU2KM=149_597_870.7; const DAY2S=86_400.0
furnsh("data/naif0012.tls","data/de440.bsp","data/sb441-n16.bsp")
et_0=(2.4621385359989386e6-2_451_545.0)*DAY2S
t_2030=10958.0372426095; et_2030=t_2030*DAY2S

# ASSIST eredu osoaren egoerak (AU, AU/egun) — extract_assist_states.py-tik
# (r, v, t_egun)
states = [
 ("AURRETIK (peri-1d)", 10694.408, [-9.255005157964148e-01,-3.596294399999496e-01,-1.568937511349967e-01], [8.858947897497478e-03,-1.294921964345618e-02,-4.566081474579382e-03]),
 ("PERIGEOA",           10695.408, [-9.163515191154751e-01,-3.724967833036590e-01,-1.613858267088595e-01], [1.033901204637969e-02,-1.312874650101358e-02,-4.532157931106631e-03]),
 ("OSTEAN (peri+1d)",   10696.408, [-9.061536142715118e-01,-3.863636036449206e-01,-1.665811441945212e-01], [1.021032177029702e-02,-1.387918877709171e-02,-5.253369156543048e-03]),
]
# ASSIST 2030 erreferentzia (AU)
assist_2030 = [1.757586529869240e-02, 1.219667338845853e+00, 4.783965553708925e-01]

# Holman IC (jatorrizko hasiera, t_initial) — erreferentziarako
sun_st,_=spkez(Int32(10),et_0,"J2000","NONE",Int32(0))
hp=[-5.5946538550488512e-1,8.5647564757574512e-1,3.0415066217102493e-1]
hv=[-1.3818324735921638e-2,-6.0088275597939191e-3,-2.5805044631309632e-3]
u0_holman=[hp[1]*AU2KM+sun_st[1],hp[2]*AU2KM+sun_st[2],hp[3]*AU2KM+sun_st[3],hv[1]*AU2KM/DAY2S+sun_st[4],hv[2]*AU2KM/DAY2S+sun_st[5],hv[3]*AU2KM/DAY2S+sun_st[6]]

bids=[10,399,301,1,2,4,5,6,7,8]; bmus=[1.32712440041279e11,398600.436,4902.800,22031.869,324858.592,42828.376,1.267127641e8,3.7940584842e7,5.7945564e6,6.836527101e6]
aids=[2_000_001,2_000_002,2_000_003,2_000_004,2_000_007,2_000_010,2_000_015,2_000_016,2_000_031,2_000_052,2_000_065,2_000_087,2_000_088,2_000_107,2_000_511,2_000_704]
amus=[62.6289,13.6659,1.9206,17.2882,1.1399,5.6251,2.0230,1.5897,1.0794,2.6830,0.9381,2.1682,1.1898,1.4437,3.8945,2.8304]
tiv=(et_0-5*DAY2S,et_2030+5*DAY2S); allids=vcat(bids,aids); fb=Dict{Int,Tuple{Int,Int}}(id=>(13,8) for id in aids); fb[399]=(13,8)
create_coeffs_file("data/coeffs_ms.json","data/coeffs_ms.csv",allids,fill(tiv,length(allids));fallback_params=fb)
bod=[BodyCoeffs("data/coeffs_ms.json","data/coeffs_ms.csv",id,tiv) for id in allids]
ph=(gr=true,j2_sun=true,j2_earth=true,j3_earth=true,j4_earth=true,yarkovsky=true,A1=5.0e-13,A2=-2.9e-14,A3=0.0)
p=(mus=vcat(bmus,amus),ids=allids,bodies=bod,physics=ph)
solve(ODEProblem(f_master!,u0_holman,(et_0,et_0+3600.0),p),IRKGL16();reltol=1e-12,abstol=1e-12,save_everystep=false)
run2030(ui,eti)=(s=solve(ODEProblem(f_master!,ui,(eti,et_2030),p),IRKGL16();reltol=1e-12,abstol=1e-12,save_everystep=false); s.u[end][1:3]./AU2KM)

println("\nGrAL+Yarkovsky → 2030, ASSIST-en egoeratik abiatuta:")
println(rpad("Hasiera-unea",22), rpad("Δ 2030 vs ASSIST",16), "perigeoa zeharkatu?")
println("-"^58)
# Jatorrizko hasiera (Holman, t_initial): perigeoa zeharkatzen du
g_full=run2030(u0_holman,et_0)
d_full=norm(g_full.-assist_2030)*AU2KM
@printf("%-22s %-16s %s\n","Holman IC (hasiera)", d_full<1000 ? @sprintf("%.1f m",d_full*1000) : @sprintf("%.1f km",d_full), "BAI (100 egun lehenago)")
for (lab,t,r,v) in states
    ui=[r[1]*AU2KM,r[2]*AU2KM,r[3]*AU2KM, v[1]*AU2KM/DAY2S,v[2]*AU2KM/DAY2S,v[3]*AU2KM/DAY2S]
    g=run2030(ui,t*DAY2S)
    d=norm(g.-assist_2030)*AU2KM
    cross = t < 10695.408 ? "BAI" : "EZ"
    @printf("%-22s %-16s %s\n", lab, d<1000 ? @sprintf("%.1f m",d*1000) : @sprintf("%.1f km",d), cross)
end
println("-"^58)
println("(Δ = GrAL+Yark posizioa 2030ean vs ASSIST eredu osoa. Hasiera-egoera ASSIST-ena denez,\n hasierako diferentzia 0 da; 2030eko Δ eredu-diferentzia (GR_SIMPLE vs GR_EIH) da, perigeoak anplifikatua.)")

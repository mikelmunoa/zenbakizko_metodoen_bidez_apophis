using StaticArrays # Gomendatua: svector-ak erabiltzeko

# ── Konpilazio-denborako konstanteak ──────────────────────────────────────────
# Module-mailako const-ak konpiladoreari bermatzen dio hauek behin kalkulatzen
# direla eta inline erabil daitezkeela, baita trig funtzioak direla ere.
const _C2       = 299792.458^2              # km²/s²  (c²)
const _J2_SUN   = 2.1961391516529825e-7    # DE440 balioa (ASSIST-ekin bat)
const _R_SUN    = 696000.0                 # km — DE440 ASUN (ASSIST-ekin bat)
const _J2_EARTH = 0.00108262539            # DE440 balioa (ASSIST-ekin bat)
const _J3_EARTH = -2.53241e-6              # DE440 balioa (ASSIST-ekin bat)
const _J4_EARTH = -1.619898e-6             # DE440 balioa (ASSIST-ekin bat)
const _R_EARTH  = 6378.1366               # km — DE440 RE (ASSIST-ekin bat)
# Eguzkiaren polo-ardatza: RA=286.13°, Dec=63.87° (Park et al. 2021, DE441)
const _CAS = cos(286.13 * π / 180.0)   # cos(RA_Sun)
const _SAS = sin(286.13 * π / 180.0)   # sin(RA_Sun)
const _CDS = cos(63.87  * π / 180.0)   # cos(Dec_Sun)
const _SDS = sin(63.87  * π / 180.0)   # sin(Dec_Sun)
# Yarkovsky / indar ez-grabitatorioa (ASSIST NON_GRAVITATIONAL): unitate-bihurketa
const _AU2KM     = 149_597_870.7            # km/AU (DE440)
const _DAY2S     = 86_400.0                 # s/egun
const _YARK_CONV = _AU2KM^3 / _DAY2S^2      # km³/s² — A[AU/egun²]·g·û → km/s²

"""
    _find_required_indices(ids)

NAIF ID hauen indizeak bilatu:
- 10: Eguzkia
- 399: Lurra (barizentrikoa, SSB-rekiko) — `f_master!`-entzat nahitaezkoa
- 301: Ilargia (barizentrikoa, SSB-rekiko)

Oharra: body 3 (EMB) ere onartzen da funtzio zaharrentzat (`f_all!`, `f_all_rel!`, ...).
`f_master!` bertan body 399 besterik ez du onartzen.
"""
@inline function _find_required_indices(ids)
    idx_sun   = 0
    idx_earth = 0  # 399 (Lurra, barizentrikoa) edo 3 (EMB, funtzio zaharrentzat)
    idx_moon  = 0
    @inbounds for j in eachindex(ids)
        idj = ids[j]
        idj == 10                && (idx_sun   = j)
        (idj == 399 || idj == 3) && (idx_earth = j)
        idj == 301               && (idx_moon  = j)
    end
    (idx_sun == 0 || idx_earth == 0 || idx_moon == 0) &&
        throw(ArgumentError("p.ids-en NAIF IDak falta dira: beharrezkoak 10 (Sun), 399 (Earth), 301 (Moon)"))
    return idx_sun, idx_earth, idx_moon
end

"""
f_all!(du, u, p, t)

Satelite baten mugimendu ekuazioak kalkulatzen ditu N-gorputzen eremuan.
- `u`: [x, y, z, vx, vy, vz]
- `p`: `mus` (grabitate parametroak), `bodies` (posizioak/funtzioak) eta `ids` (NAIF ID-ak)
"""
function f_all!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko (destructuring)
    # @views erabiliz ez dugu memoria berririk erreserbatzen
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    ax, ay, az = 0.0, 0.0, 0.0

    # ============================================================
    # GARRANTZITSUA: Earth(t) = EMB (NAIF 3), ez Lurra bera!
    # Moon(t) = Ilargiaren posizioa EMB-rekiko (NAIF 301).
    # Lurraren posizio erreala: EMB - (mu_M/mu_E) * Moon_rel_EMB
    # Ilargiaren posizio bary:  EMB + Moon_rel_EMB
    # p.bodies ordena: [Sun, Earth(EMB), Moon(rel-EMB), ...]
    # ============================================================
    idx_sun, idx_emb, idx_moon = _find_required_indices(p.ids)
    emb_pos = p.bodies[idx_emb](t, 1)
    moon_rel = p.bodies[idx_moon](t, 1)
    mu_ratio = p.mus[idx_moon] / p.mus[idx_emb]

    # 2. Grabitatearen kalkulua (Newton hutsa)
    # zip + Tuple erabiltzen da: konpiladoreak luzera ezagutzen du,
    # loop-a unroll egin dezake, eta 0 heap alokazio sortzen dira.
    @inbounds for (i, (mu_i, body_i)) in enumerate(zip(p.mus, p.bodies))
        if i == idx_emb
            bx = emb_pos[1] - mu_ratio * moon_rel[1]
            by = emb_pos[2] - mu_ratio * moon_rel[2]
            bz = emb_pos[3] - mu_ratio * moon_rel[3]
        elseif i == idx_moon
            bx = emb_pos[1] + moon_rel[1]
            by = emb_pos[2] + moon_rel[2]
            bz = emb_pos[3] + moon_rel[3]
        else
            bx, by, bz = body_i(t, 1)
        end

        dx = bx - x
        dy = by - y
        dz = bz - z

        # muladd: FMA instrukzioa erabiltzen du (azkarragoa eta zehatzagoa)
        r_sq = muladd(dx, dx, muladd(dy, dy, dz * dz))

        inv_r = 1.0 / sqrt(r_sq)
        common = mu_i * inv_r * inv_r * inv_r

        ax = muladd(dx, common, ax)
        ay = muladd(dy, common, ay)
        az = muladd(dz, common, az)
    end

    # 3. Deribatuak esleitu
    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] = ax
    du[5] = ay
    du[6] = az

    return nothing
end



"""
    f_all_rel!(du, u, p, t)

Satelite baten mugimendu ekuazioak N-gorputzen eremuan,
Eguzkiaren efektu erlatibistarekin (Damour & Deruelle, 1985).

    a_REL = (GM☉ / (c² |r|³)) * [(4GM☉/|r| - |v|²)r + 4(r·v)v]

non r eta v koordenatu heliozentrikoak diren.
- `u`: [x, y, z, vx, vy, vz] (barizentrikoak)
- `p`: `mus` (grabitate parametroak), `bodies` (posizioak/funtzioak) eta `ids` (NAIF ID-ak)
"""
function f_all_rel!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko (destructuring)
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    ax, ay, az = 0.0, 0.0, 0.0

    # Earth(EMB) eta Moon(rel-EMB) -> Earth bary + Moon bary
    idx_sun, idx_emb, idx_moon = _find_required_indices(p.ids)
    emb_state = p.bodies[idx_emb](t)
    moon_rel = p.bodies[idx_moon](t)
    mu_ratio = p.mus[idx_moon] / p.mus[idx_emb]

    # 2. Grabitatearen kalkulua (Newtondarra + erlatibista Eguzkiarentzat)
    @inbounds for (i, (mu_i, body_i)) in enumerate(zip(p.mus, p.bodies))
        if i == idx_emb
            state_i = emb_state .- mu_ratio .* moon_rel
        elseif i == idx_moon
            state_i = emb_state .+ moon_rel
        else
            state_i = body_i(t)
        end
        bx, by, bz = state_i[1], state_i[2], state_i[3]

        dx = bx - x
        dy = by - y
        dz = bz - z

        r_sq = muladd(dx, dx, muladd(dy, dy, dz * dz))

        inv_r = 1.0 / sqrt(r_sq)
        common = mu_i * inv_r * inv_r * inv_r

        ax = muladd(dx, common, ax)
        ay = muladd(dy, common, ay)
        az = muladd(dz, common, az)

        # =======================================================
        # EFEKTU ERLATIBISTA (Damour & Deruelle, 1985)
        # Bakarrik Eguzkiarentzat (NAIF 10)
        # =======================================================
        if i == idx_sun
            # Eguzkiaren abiadura (SVector-aren 4., 5. eta 6. osagaiak)
            bvx, bvy, bvz = state_i[4], state_i[5], state_i[6]

            # Apophisen egoera HELIOZENTRIKOA (Apophis - Eguzkia)
            rx_h = x - bx
            ry_h = y - by
            rz_h = z - bz
            vx_h = vx - bvx
            vy_h = vy - bvy
            vz_h = vz - bvz

            # Magnitudeak eta biderkadura eskalarra
            r_h_norm = sqrt(muladd(rx_h, rx_h, muladd(ry_h, ry_h, rz_h * rz_h)))
            v_h_sq   = muladd(vx_h, vx_h, muladd(vy_h, vy_h, vz_h * vz_h))
            r_dot_v  = muladd(rx_h, vx_h, muladd(ry_h, vy_h, rz_h * vz_h))

            # a_REL = (GM☉ / (c² |r|³)) * [(4GM☉/|r| - |v|²)r + 4(r·v)v]
            coef  = mu_i / (_C2 * r_h_norm^3)
            term1 = muladd(4.0, mu_i / r_h_norm, -v_h_sq)

            ax += coef * muladd(term1, rx_h, 4.0 * r_dot_v * vx_h)
            ay += coef * muladd(term1, ry_h, 4.0 * r_dot_v * vy_h)
            az += coef * muladd(term1, rz_h, 4.0 * r_dot_v * vz_h)
        end
    end

    # 3. Deribatuak esleitu
    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] = ax
    du[5] = ay
    du[6] = az

    return nothing
end




"""
    f_all_rel_J!(du, u, p, t)

Satelite baten mugimendu ekuazioak N-gorputzen eremuan:
- Newton grabitazioa gorputz guztientzat
- Efektu erlatibista Eguzkiarentzat (Damour & Deruelle, 1985)
- Harmoniko grabitatorioak: Eguzkiaren J2 eta Lurraren J2, J3, J4 (ASSIST, Park et al. 2021)

- `u`: [x, y, z, vx, vy, vz] (barizentrikoak, km eta km/s)
- `p`: `mus` (grabitate parametroak), `bodies` (posizioak/funtzioak) eta `ids` (NAIF ID-ak)
"""
function f_all_rel_J!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # Konstanteak moduluko const-etatik (_C2, _J2_SUN, _R_SUN, _J2_EARTH, ...)

    # ============================================================
    # KOORDENATU-BIRAKETA: Jn harmonikoak gorputzaren biraketa-
    # ardatzarekiko (polo) definitzen dira, ez ekliptikaren z-rekiko!
    #
    # Lurraren oblikuotasuna: ε = 23.4393° (J2000)
    # Ekliptika → Ekuatorial: R_x(ε)
    #   dx_eq = dx
    #   dy_eq =  cosε·dy - sinε·dz
    #   dz_eq =  sinε·dy + cosε·dz
    #
    # Ekuatorial → Ekliptika: R_x(-ε) = R_x(ε)ᵀ
    #   ax_ecl = ax_eq
    #   ay_ecl =  cosε·ay_eq + sinε·az_eq
    #   az_ecl = -sinε·ay_eq + cosε·az_eq
    #
    # Eguzkiaren J2 oso txikia da (2.2e-7) eta poloa ~7.25°-ra dago
    # ekliptikatik, beraz ekliptikaren z-a erabiltzen dugu (hurbilketa ona).
    # ============================================================
    cosε = 0.9174820620691818   # cos(23.4393°)
    sinε = 0.3977771559319137   # sin(23.4393°)

    ax, ay, az = 0.0, 0.0, 0.0

    # Earth(EMB) eta Moon(rel-EMB) -> Earth bary + Moon bary
    idx_sun, idx_emb, idx_moon = _find_required_indices(p.ids)
    emb_state = p.bodies[idx_emb](t)
    moon_rel = p.bodies[idx_moon](t)
    mu_ratio = p.mus[idx_moon] / p.mus[idx_emb]

    # 2. Grabitatearen kalkulua
    @inbounds for (i, (mu_i, body_i)) in enumerate(zip(p.mus, p.bodies))
        if i == idx_emb
            state_i = emb_state .- mu_ratio .* moon_rel
        elseif i == idx_moon
            state_i = emb_state .+ moon_rel
        else
            state_i = body_i(t)
        end
        bx, by, bz = state_i[1], state_i[2], state_i[3]

        dx = bx - x   # gorputzetik partikulara: -(partikula - gorputza)
        dy = by - y
        dz = bz - z

        r_sq = muladd(dx, dx, muladd(dy, dy, dz * dz))
        inv_r = 1.0 / sqrt(r_sq)
        inv_r2 = inv_r * inv_r
        inv_r3 = inv_r2 * inv_r

        # --- Newton ---
        common = mu_i * inv_r3
        ax = muladd(dx, common, ax)
        ay = muladd(dy, common, ay)
        az = muladd(dz, common, az)

        # =======================================================
        # EGUZKIA (NAIF 10): Erlatibista + J2
        # =======================================================
        if i == idx_sun
            # --- Efektu erlatibista (Damour & Deruelle) ---
            bvx, bvy, bvz = state_i[4], state_i[5], state_i[6]
            rx_h = x - bx;  ry_h = y - by;  rz_h = z - bz
            vx_h = vx - bvx; vy_h = vy - bvy; vz_h = vz - bvz

            r_h_norm = sqrt(muladd(rx_h, rx_h, muladd(ry_h, ry_h, rz_h * rz_h)))
            v_h_sq   = muladd(vx_h, vx_h, muladd(vy_h, vy_h, vz_h * vz_h))
            r_dot_v  = muladd(rx_h, vx_h, muladd(ry_h, vy_h, rz_h * vz_h))

            coef  = mu_i / (_C2 * r_h_norm^3)
            term1 = muladd(4.0, mu_i / r_h_norm, -v_h_sq)

            ax += coef * muladd(term1, rx_h, 4.0 * r_dot_v * vx_h)
            ay += coef * muladd(term1, ry_h, 4.0 * r_dot_v * vy_h)
            az += coef * muladd(term1, rz_h, 4.0 * r_dot_v * vz_h)

            # # --- Eguzkiaren J2 ---
            # # J2_Sun = 2.2e-7 oso txikia da, eta ekliptikaren z erabiltzen
            # # dugu hurbilketa gisa (Eguzkiaren poloa ~7.25° ekliptikatik)
            # z2_r2 = dz^2 * inv_r2
            # J2_coef = 1.5 * mu_i * J2_Sun * R_Sun^2 * inv_r3 * inv_r2
            # fxy = J2_coef * (1.0 - 5.0 * z2_r2)
            # fz  = J2_coef * (3.0 - 5.0 * z2_r2)
            # ax += dx * fxy
            # ay += dy * fxy
            # az += dz * fz

        # =======================================================
        # LURRA (NAIF 3, EMB zuzenduta): J2, J3, J4 — POLE ERREFERENTZIAN
        # =======================================================
        elseif i == idx_emb
            # Lurraren poloa: RAe=0°, Dece=90° (J2000 ekuatorian)
            # Biraketa sinplifikatua: cos(0°)=1, sin(0°)=0, cos(90°)≈0, sin(90°)=1
            #   dxp_e = dy,  dyp_e = -dx,  dzp_e = dz
            # Itzulera: ax -= resy_e,  ay += resx_e,  az += resz_e
            dxp_e = -dy
            dyp_e =  dx
            dzp_e = -dz

            # Kalkuluak Lurraren ekuatorialean
            r = sqrt(r_sq)
            costheta2_e = dzp_e^2 * inv_r2

            # --- Lurraren J2 ---
            J2_prefac_e = 3.0 * _J2_EARTH * _R_EARTH^2 / (r_sq * r_sq * r * 2.0)
            J2_fac_e = 5.0 * costheta2_e - 1.0

            resx_e = mu_i * J2_prefac_e * J2_fac_e * dxp_e
            resy_e = mu_i * J2_prefac_e * J2_fac_e * dyp_e
            resz_e = mu_i * J2_prefac_e * (J2_fac_e - 2.0) * dzp_e

            # --- Lurraren J3 ---
            J3_prefac_e = 5.0 * _J3_EARTH * _R_EARTH^3 / (r_sq * r_sq * r * 2.0)
            J3_fac_e = 3.0 - 7.0 * costheta2_e

            resx_e += -mu_i * J3_prefac_e * (1.0 / r_sq) * J3_fac_e * dxp_e * dzp_e
            resy_e += -mu_i * J3_prefac_e * (1.0 / r_sq) * J3_fac_e * dyp_e * dzp_e
            resz_e += -mu_i * J3_prefac_e * (6.0 * costheta2_e - 7.0 * costheta2_e^2 - 0.6)

            # --- Lurraren J4 ---
            J4_prefac_e = 5.0 * _J4_EARTH * _R_EARTH^4 / (r_sq * r_sq * r_sq * r * 8.0)
            J4_fac_e = 63.0 * costheta2_e^2 - 42.0 * costheta2_e + 3.0

            resx_e += mu_i * J4_prefac_e * J4_fac_e * dxp_e
            resy_e += mu_i * J4_prefac_e * J4_fac_e * dyp_e
            resz_e += mu_i * J4_prefac_e * (J4_fac_e + 12.0 - 28.0 * costheta2_e) * dzp_e

            # Itzulera sinplifikatua (RAe=0°, Dece=90°):
            ax -= resy_e
            ay += resx_e
            az += resz_e
        end
    end

    # 3. Deribatuak esleitu
    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] = ax
    du[5] = ay
    du[6] = az

    return nothing
end


"""
    dd_perturbation(u, sun_state, mu_sun, [c2])

Damour-Deruelle erlatibista perturbazio azelerazioa kalkulatzen du (Eguzki-efektua soilik).
`f_master!`-en barneko kalkulua isolatzen du, proba eta analisi koadernoentzat.

- `u`         : partikulararen egoera barizentrikoa [x,y,z,vx,vy,vz] (km, km/s)
- `sun_state` : Eguzkiaren egoera barizentrikoa [x,y,z,vx,vy,vz] (km, km/s)
- `mu_sun`    : Eguzkiaren grabitazio parametroa (km³/s²)
- `c2`        : argiaren abiadura karratua (km/s)²  [default: 299792.458²]

Itzultzen du: SVector{3,Float64} — perturbazio azelerazioa (km/s²)
"""
function dd_perturbation(u, sun_state, mu_sun, c2 = 299792.458^2)
    rx_h = u[1] - sun_state[1]
    ry_h = u[2] - sun_state[2]
    rz_h = u[3] - sun_state[3]
    vx_h = u[4] - sun_state[4]
    vy_h = u[5] - sun_state[5]
    vz_h = u[6] - sun_state[6]

    r_h     = sqrt(muladd(rx_h, rx_h, muladd(ry_h, ry_h, rz_h * rz_h)))
    v_sq    = muladd(vx_h, vx_h, muladd(vy_h, vy_h, vz_h * vz_h))
    r_dot_v = muladd(rx_h, vx_h, muladd(ry_h, vy_h, rz_h * vz_h))

    coef  = mu_sun / (c2 * r_h^3)
    term1 = muladd(4.0, mu_sun / r_h, -v_sq)

    return SVector(
        coef * muladd(term1, rx_h, 4.0 * r_dot_v * vx_h),
        coef * muladd(term1, ry_h, 4.0 * r_dot_v * vy_h),
        coef * muladd(term1, rz_h, 4.0 * r_dot_v * vz_h),
    )
end


"""
    f_master!(du, u, p, t)

Funtzio nagusia: satelite baten mugimendu ekuazioak N-gorputzen eremuan,
fisika-efektu aukeragarriekin — ASSIST-en moduan.

`p` parametroetan `physics` eremu aukerazkoa (NamedTuple):

    physics = (
        gr        = true,   # Damour-Deruelle erlatibitatearen efektua (Eguzkia)
        j2_sun    = false,  # Eguzkiaren J2 harmonikoa (Park et al. 2021)
        j2_earth  = false,  # Lurraren J2 harmonikoa
        j3_earth  = false,  # Lurraren J3 harmonikoa
        j4_earth  = false,  # Lurraren J4 harmonikoa
        yarkovsky = false,  # Yarkovsky/indar ez-grabitatorioa (ASSIST NON_GRAVITATIONAL)
        A1 = 0.0, A2 = 0.0, A3 = 0.0,   # Yarkovsky parametroak (AU/egun²); yarkovsky=true bada
    )

`yarkovsky` aukeran bada (`get`-ekin irakurria), `A1/A2/A3` ere egon behar dira
`physics`-en. Apophis: A1=5.0e-13, A2=-2.9e-14, A3=0.0.

`physics` ematen ez bada, balioak aurrenak dira: GR aktibo, harmonikorik gabe
(hots, `f_all_rel!`-en baliokidea).

Adibidea:
    p = (mus     = [mu_S, mu_E, mu_M],
         ids     = [10, 399, 301],          # body 399 (Lurra, barizentrikoa) nahitaez
         physics = (gr=true, j2_sun=true, j2_earth=true, j3_earth=true, j4_earth=false))
    RK4(u0, t0, t_end, h, p, f_master!, m_out)

OHARRA: `p.ids`-ek 10 (Eguzkia), 399 edo 3 (Lurra/EMB), eta 301 (Ilargia) eduki behar ditu.
`p.bodies` (LittleEphemeris BodyCoeffs) beharrezkoa da: posizioak `p.bodies[i](t)` bidez jasotzen dira.
  - ids[idx_earth] == 399 → Earth barizentrikoa zuzenean; Moon barizentrikoa zuzenean
  - ids[idx_earth] == 3   → EMB barizentrikoa; Ilargia EMB-rekiko (LittleEphemeris konbentzioa)
Bi kasuetan J2/J3/J4 harmonikoak aplikatzen dira idx_earth gorputzarentzat.
"""
Base.@constprop :aggressive function f_master!(du, u, p, t)
    x, y, z    = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # ── Fisika konfigurazioa (aurrekoak: GR bai, harmonikorik ez) ────────────
    cfg = hasproperty(p, :physics) ? p.physics :
          (gr=true, j2_sun=false, j2_earth=false, j3_earth=false, j4_earth=false)

    # ── Konstanteak moduluko const-etatik (_C2, _J2_SUN, _R_SUN, _CAS, _SAS, ...) ──
    # Lurraren poloa: RAe=0°, Dece=90° → dxp_e=dy, dyp_e=-dx, dzp_e=dz
    # Itzulera: ax -= resy_e,  ay += resx_e,  az += resz_e

    ax, ay, az = 0.0, 0.0, 0.0

    idx_sun, idx_earth, idx_moon = _find_required_indices(p.ids)

    # Earth/Moon posizioak p.bodies[i](t) bidez (0-alokazio SVector, ez spkez):
    #   ids[idx_earth] == 399 → Earth eta Moon biak barizentrikoak zuzenean (gomendatua)
    #   ids[idx_earth] == 3   → EMB barizentrikoa; Moon BARIZENTRIKOA ere bai
    #     (LittleEphemeris-ek spkez(301,...,0) erabiltzen du → Moon barizentrikoa,
    #      ez EMB-rekiko erlatiboa. moon_rel = moon_bary - emb bidez zuzentzen da.)
    if p.ids[idx_earth] == 3
        emb_state   = p.bodies[idx_earth](t)
        moon_bary   = p.bodies[idx_moon](t)
        moon_rel    = moon_bary .- emb_state          # Moon EMB-rekiko benetakoa
        mu_ratio    = p.mus[idx_moon] / (p.mus[idx_earth] + p.mus[idx_moon])
        earth_state = emb_state .- mu_ratio .* moon_rel
        moon_state  = moon_bary
    else  # 399
        earth_state = p.bodies[idx_earth](t)
        moon_state  = p.bodies[idx_moon](t)
    end

    @inbounds for i in eachindex(p.mus)
        mu_i = p.mus[i]
        if i == idx_earth
            state_i = earth_state
        elseif i == idx_moon
            state_i = moon_state
        else
            state_i = p.bodies[i](t)
        end
        bx, by, bz = state_i[1], state_i[2], state_i[3]

        dx = bx - x;  dy = by - y;  dz = bz - z
        r_sq   = muladd(dx, dx, muladd(dy, dy, dz * dz))
        r      = sqrt(r_sq)
        inv_r  = 1.0 / r
        inv_r2 = inv_r * inv_r
        inv_r3 = inv_r2 * inv_r

        # ── Newton ──────────────────────────────────────────────────────────
        common = mu_i * inv_r3
        ax = muladd(dx, common, ax)
        ay = muladd(dy, common, ay)
        az = muladd(dz, common, az)

        # ── Eguzkiaren efektuak (NAIF 10) ────────────────────────────────────
        if i == idx_sun
            bvx, bvy, bvz = state_i[4], state_i[5], state_i[6]
            rx_h = x - bx;  ry_h = y - by;  rz_h = z - bz
            vx_h = vx - bvx; vy_h = vy - bvy; vz_h = vz - bvz
            r_h  = sqrt(muladd(rx_h, rx_h, muladd(ry_h, ry_h, rz_h * rz_h)))

            # Damour-Deruelle
            if cfg.gr
                v_sq  = muladd(vx_h, vx_h, muladd(vy_h, vy_h, vz_h * vz_h))
                rdv   = muladd(rx_h, vx_h, muladd(ry_h, vy_h, rz_h * vz_h))
                coef  = mu_i / (_C2 * r_h^3)
                t1    = muladd(4.0, mu_i / r_h, -v_sq)
                ax   += coef * muladd(t1, rx_h, 4.0 * rdv * vx_h)
                ay   += coef * muladd(t1, ry_h, 4.0 * rdv * vy_h)
                az   += coef * muladd(t1, rz_h, 4.0 * rdv * vz_h)
            end

            # Eguzkiaren J2
            if cfg.j2_sun
                # ASSIST konbentzioa: d = asteroidea − Eguzkia (rx_h), EZ Eguzkia − ast (dx).
                # J2 azelerazioa d-rekiko BAKOITIA da → zeinu egokia rx_h erabiliz.
                dxp = muladd(-rx_h, _SAS, ry_h * _CAS)
                dyp = muladd(-rx_h, _CAS * _SDS, muladd(-ry_h, _SAS * _SDS, rz_h * _CDS))
                dzp = muladd( rx_h, _CAS * _CDS, muladd( ry_h, _SAS * _CDS, rz_h * _SDS))
                cth2   = dzp * dzp * inv_r2
                J2pref = 3.0 * _J2_SUN * _R_SUN^2 * inv_r3 * inv_r2 / 2.0
                J2f    = muladd(5.0, cth2, -1.0)
                resx   = mu_i * J2pref * J2f * dxp
                resy   = mu_i * J2pref * J2f * dyp
                resz   = mu_i * J2pref * muladd(J2f, 1.0, -2.0) * dzp
                ax += muladd(-resx, _SAS, muladd(-resy, _CAS * _SDS, resz * _CAS * _CDS))
                ay += muladd( resx, _CAS, muladd(-resy, _SAS * _SDS, resz * _SAS * _CDS))
                az += muladd(resy, _CDS, resz * _SDS)
            end

            # ── Yarkovsky / indar ez-grabitatorioa (ASSIST NON_GRAVITATIONAL) ──
            #   a = (1/r²)·(A1·r̂ + A2·t̂ + A3·n̂);  h = d×v (normala), t = h×d (transbertsala)
            #   A AU/egun² → km/s²:  × _YARK_CONV/r_h²  (r̂,t̂,n̂ eskala-independenteak)
            if get(cfg, :yarkovsky, false)
                A1y = cfg.A1;  A2y = cfg.A2;  A3y = cfg.A3
                inv_rh = 1.0 / r_h
                ur1 = rx_h * inv_rh;  ur2 = ry_h * inv_rh;  ur3 = rz_h * inv_rh   # r̂
                hx = ry_h * vz_h - rz_h * vy_h                                     # h = d×v
                hy = rz_h * vx_h - rx_h * vz_h
                hz = rx_h * vy_h - ry_h * vx_h
                inv_h = 1.0 / sqrt(muladd(hx, hx, muladd(hy, hy, hz * hz)))
                un1 = hx * inv_h;  un2 = hy * inv_h;  un3 = hz * inv_h             # n̂
                tx = hy * rz_h - hz * ry_h                                         # t = h×d
                ty = hz * rx_h - hx * rz_h
                tz = hx * ry_h - hy * rx_h
                inv_t = 1.0 / sqrt(muladd(tx, tx, muladd(ty, ty, tz * tz)))
                ut1 = tx * inv_t;  ut2 = ty * inv_t;  ut3 = tz * inv_t             # t̂
                gc = _YARK_CONV * inv_rh * inv_rh                                  # (1/r²)·bihurketa
                ax += gc * muladd(A1y, ur1, muladd(A2y, ut1, A3y * un1))
                ay += gc * muladd(A1y, ur2, muladd(A2y, ut2, A3y * un2))
                az += gc * muladd(A1y, ur3, muladd(A2y, ut3, A3y * un3))
            end

        # ── Lurraren efektuak (NAIF 399) ─────────────────────────────────────
        elseif i == idx_earth
            if cfg.j2_earth || cfg.j3_earth || cfg.j4_earth
                # RAe=0°, Dece=90°; d = ast − Lurra (ASSIST konbentzioa, EZ Lurra − ast).
                # J2/J4 d-rekiko bakoitiak dira → zeinu egokia d = ast−Lurra erabiliz.
                dxp_e  = -dy
                dyp_e  =  dx
                dzp_e  = -dz
                cth2_e = dzp_e * dzp_e * inv_r2

                resx_e = 0.0;  resy_e = 0.0;  resz_e = 0.0

                if cfg.j2_earth
                    J2pref_e = 3.0 * _J2_EARTH * _R_EARTH^2 * inv_r3 * inv_r2 / 2.0
                    J2f_e    = muladd(5.0, cth2_e, -1.0)
                    resx_e  += mu_i * J2pref_e * J2f_e * dxp_e
                    resy_e  += mu_i * J2pref_e * J2f_e * dyp_e
                    resz_e  += mu_i * J2pref_e * muladd(J2f_e, 1.0, -2.0) * dzp_e
                end

                if cfg.j3_earth
                    J3pref_e = 5.0 * _J3_EARTH * _R_EARTH^3 * inv_r3 * inv_r2 / 2.0
                    J3f_e    = muladd(-7.0, cth2_e, 3.0)
                    resx_e  += -mu_i * J3pref_e * inv_r2 * J3f_e * dxp_e * dzp_e
                    resy_e  += -mu_i * J3pref_e * inv_r2 * J3f_e * dyp_e * dzp_e
                    resz_e  += -mu_i * J3pref_e * muladd(muladd(-7.0, cth2_e, 6.0), cth2_e, -0.6)
                end

                if cfg.j4_earth
                    J4pref_e = 5.0 * _J4_EARTH * _R_EARTH^4 * inv_r3 * inv_r2 * inv_r2 / 8.0
                    J4f_e    = muladd(muladd(63.0, cth2_e, -42.0), cth2_e, 3.0)
                    resx_e  += mu_i * J4pref_e * J4f_e * dxp_e
                    resy_e  += mu_i * J4pref_e * J4f_e * dyp_e
                    resz_e  += mu_i * J4pref_e * muladd(J4f_e + 12.0, 1.0, -28.0 * cth2_e) * dzp_e
                end

                # RAe=0°, Dece=90° itzulera: ax -= resy_e, ay += resx_e, az += resz_e
                ax -= resy_e
                ay += resx_e
                az += resz_e
            end
        end
    end

    du[1] = vx;  du[2] = vy;  du[3] = vz
    du[4] = ax;  du[5] = ay;  du[6] = az
    return nothing
end


function f_BigFloat!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko (destructuring)
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # BigFloat aritmetika: dena BigFloat
    ax, ay, az = BigFloat(0), BigFloat(0), BigFloat(0)

    # Earth(EMB) eta Moon(rel-EMB) -> Earth bary + Moon bary
    idx_sun, idx_emb, idx_moon = _find_required_indices(p.ids)
    emb_pos = p.bodies[idx_emb](t)
    moon_rel = p.bodies[idx_moon](t)
    mu_ratio = p.mus[idx_moon] / p.mus[idx_emb]

    # 2. Grabitatearen kalkulua
    @inbounds for (i, (mu_i, body_i)) in enumerate(zip(p.mus, p.bodies))
        if i == idx_emb
            bx = emb_pos[1] - mu_ratio * moon_rel[1]
            by = emb_pos[2] - mu_ratio * moon_rel[2]
            bz = emb_pos[3] - mu_ratio * moon_rel[3]
        elseif i == idx_moon
            bx = emb_pos[1] + moon_rel[1]
            by = emb_pos[2] + moon_rel[2]
            bz = emb_pos[3] + moon_rel[3]
        else
            bx, by, bz = body_i(t)
        end

        dx = bx - x
        dy = by - y
        dz = bz - z

        r_sq = muladd(dx, dx, muladd(dy, dy, dz * dz))

        inv_r = BigFloat(1) / sqrt(r_sq)
        common = mu_i * inv_r * inv_r * inv_r

        ax = muladd(dx, common, ax)
        ay = muladd(dy, common, ay)
        az = muladd(dz, common, az)
    end

    # 3. Deribatuak esleitu
    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] = ax
    du[5] = ay
    du[6] = az

    return nothing
end


"""
    f_all_rel_J2_Sun!(du, u, p, t)

Satelite baten mugimendu ekuazioak N-gorputzen eremuan:
- Newton grabitazioa gorputz guztientzat
- Efektu erlatibista Eguzkiarentzat (Damour & Deruelle, 1985)
- Eguzkiaren J2 harmonikoa rigorous frame rotation-ekin (Park et al. 2021)

Eguzkiaren poloa (RA, Dec) zuzen bihurturiko koordenatuetan kalkulatzen da J2 harmonikoa.

- `u`: [x, y, z, vx, vy, vz] (barizentrrikoak, km eta km/s)
- `p`: `mus` (grabitate parametroak), `bodies` (posizioak/funtzioak) eta `ids` (NAIF ID-ak)
"""
function f_all_rel_J2_Sun!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # Konstanteak moduluko const-etatik (_C2, _J2_SUN, _R_SUN, _CAS, _SAS, _CDS, _SDS)
    
    ax, ay, az = 0.0, 0.0, 0.0

    # Earth(EMB) eta Moon(rel-EMB) -> Earth bary + Moon bary
    idx_sun, idx_emb, idx_moon = _find_required_indices(p.ids)
    emb_state = p.bodies[idx_emb](t)
    moon_rel = p.bodies[idx_moon](t)
    mu_ratio = p.mus[idx_moon] / p.mus[idx_emb]

    # 2. Grabitatearen kalkulua
    @inbounds for (i, (mu_i, body_i)) in enumerate(zip(p.mus, p.bodies))
        if i == idx_emb
            state_i = emb_state .- mu_ratio .* moon_rel
        elseif i == idx_moon
            state_i = emb_state .+ moon_rel
        else
            state_i = body_i(t)
        end
        bx, by, bz = state_i[1], state_i[2], state_i[3]

        dx = bx - x
        dy = by - y
        dz = bz - z

        r_sq = muladd(dx, dx, muladd(dy, dy, dz * dz))
        r = sqrt(r_sq)
        inv_r = 1.0 / r
        inv_r2 = inv_r * inv_r
        inv_r3 = inv_r2 * inv_r

        # --- Newton ---
        common = mu_i * inv_r3
        ax = muladd(dx, common, ax)
        ay = muladd(dy, common, ay)
        az = muladd(dz, common, az)

        # =======================================================
        # EGUZKIA (NAIF 10): Erlatibista + J2 (Solar frame rotation)
        # =======================================================
        if i == idx_sun
            # --- Efektu erlatibista (Damour & Deruelle) ---
            bvx, bvy, bvz = state_i[4], state_i[5], state_i[6]
            rx_h = x - bx;  ry_h = y - by;  rz_h = z - bz
            vx_h = vx - bvx; vy_h = vy - bvy; vz_h = vz - bvz

            r_h_norm = sqrt(muladd(rx_h, rx_h, muladd(ry_h, ry_h, rz_h * rz_h)))
            v_h_sq   = muladd(vx_h, vx_h, muladd(vy_h, vy_h, vz_h * vz_h))
            r_dot_v  = muladd(rx_h, vx_h, muladd(ry_h, vy_h, rz_h * vz_h))

            coef  = mu_i / (_C2 * r_h_norm^3)
            term1 = muladd(4.0, mu_i / r_h_norm, -v_h_sq)

            ax += coef * muladd(term1, rx_h, 4.0 * r_dot_v * vx_h)
            ay += coef * muladd(term1, ry_h, 4.0 * r_dot_v * vy_h)
            az += coef * muladd(term1, rz_h, 4.0 * r_dot_v * vz_h)

            # --- Eguzkiaren J2 (Solar frame rotation) ---
            # ASSIST konbentzioa: d = ast − Eguzkia (rx_h), EZ Eguzkia − ast (dx). J2 d-rekiko bakoitia.
            dxp = muladd(-rx_h, _SAS, ry_h * _CAS)
            dyp = muladd(-rx_h, _CAS * _SDS, muladd(-ry_h, _SAS * _SDS, rz_h * _CDS))
            dzp = muladd( rx_h, _CAS * _CDS, muladd( ry_h, _SAS * _CDS, rz_h * _SDS))

            costheta2 = dzp * dzp * inv_r2
            J2_prefac = 3.0 * _J2_SUN * _R_SUN^2 / (r_sq * r_sq * r * 2.0)
            J2_fac    = muladd(5.0, costheta2, -1.0)

            resx = mu_i * J2_prefac * J2_fac * dxp
            resy = mu_i * J2_prefac * J2_fac * dyp
            resz = mu_i * J2_prefac * muladd(J2_fac, 1.0, -2.0) * dzp

            resxp = muladd(-resx, _SAS, muladd(-resy, _CAS * _SDS, resz * _CAS * _CDS))
            resyp = muladd( resx, _CAS, muladd(-resy, _SAS * _SDS, resz * _SAS * _CDS))
            reszp = muladd(resy, _CDS, resz * _SDS)

            ax += resxp
            ay += resyp
            az += reszp
        end
    end

    # 3. Deribatuak esleitu
    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] = ax
    du[5] = ay
    du[6] = az

    return nothing
end


"""
    f_all_rel_J2_Sun_Earth!(du, u, p, t)

Satelite baten mugimendu ekuazioak N-gorputzen eremuan:
- Newton grabitazioa gorputz guztientzat
- Efektu erlatibista Eguzkiarentzat (Damour & Deruelle, 1985)
- Eguzkiaren J2 eta Lurraren J2, J3, J4 harmonikoak (pole-based, Park et al. 2021)

Eguzkiaren eta Lurraren poloak zuzen bihurturiko koordenatuetan kalkulatzen dira J harmonikoak.

- `u`: [x, y, z, vx, vy, vz] (barizentrrikoak, km eta km/s)
- `p`: `mus` (grabitate parametroak), `bodies` (posizioak/funtzioak) eta `ids` (NAIF ID-ak)
"""
function f_all_rel_J2_Sun_Earth!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # Konstanteak moduluko const-etatik
    # Lurraren poloa: RAe=0°, Dece=90° → dxp_e=dy, dyp_e=-dx, dzp_e=dz (sinplifikatua)
    
    ax, ay, az = 0.0, 0.0, 0.0

    # Earth(EMB) eta Moon(rel-EMB) -> Earth bary + Moon bary
    idx_sun, idx_emb, idx_moon = _find_required_indices(p.ids)
    emb_state = p.bodies[idx_emb](t)
    moon_rel = p.bodies[idx_moon](t)
    mu_ratio = p.mus[idx_moon] / p.mus[idx_emb]

    # 2. Grabitatearen kalkulua
    @inbounds for (i, (mu_i, body_i)) in enumerate(zip(p.mus, p.bodies))
        if i == idx_emb
            state_i = emb_state .- mu_ratio .* moon_rel
        elseif i == idx_moon
            state_i = emb_state .+ moon_rel
        else
            state_i = body_i(t)
        end
        bx, by, bz = state_i[1], state_i[2], state_i[3]

        dx = bx - x
        dy = by - y
        dz = bz - z

        r_sq = muladd(dx, dx, muladd(dy, dy, dz * dz))
        r = sqrt(r_sq)
        inv_r = 1.0 / r
        inv_r2 = inv_r * inv_r
        inv_r3 = inv_r2 * inv_r
        # --- Newton ---
        common = mu_i * inv_r3
        ax = muladd(dx, common, ax)
        ay = muladd(dy, common, ay)
        az = muladd(dz, common, az)

        # =======================================================
        # EGUZKIA (NAIF 10): Erlatibista + J2
        # =======================================================
        if i == idx_sun
            # --- Efektu erlatibista (Damour & Deruelle) ---
            bvx, bvy, bvz = state_i[4], state_i[5], state_i[6]
            rx_h = x - bx;  ry_h = y - by;  rz_h = z - bz
            vx_h = vx - bvx; vy_h = vy - bvy; vz_h = vz - bvz

            r_h_norm = sqrt(muladd(rx_h, rx_h, muladd(ry_h, ry_h, rz_h * rz_h)))
            v_h_sq   = muladd(vx_h, vx_h, muladd(vy_h, vy_h, vz_h * vz_h))
            r_dot_v  = muladd(rx_h, vx_h, muladd(ry_h, vy_h, rz_h * vz_h))

            coef  = mu_i / (_C2 * r_h_norm^3)
            term1 = muladd(4.0, mu_i / r_h_norm, -v_h_sq)

            ax += coef * muladd(term1, rx_h, 4.0 * r_dot_v * vx_h)
            ay += coef * muladd(term1, ry_h, 4.0 * r_dot_v * vy_h)
            az += coef * muladd(term1, rz_h, 4.0 * r_dot_v * vz_h)

            # --- Eguzkiaren J2 (Solar frame rotation) ---
            # ASSIST konbentzioa: d = ast − Eguzkia (rx_h), EZ Eguzkia − ast (dx). J2 d-rekiko bakoitia.
            dxp = muladd(-rx_h, _SAS, ry_h * _CAS)
            dyp = muladd(-rx_h, _CAS * _SDS, muladd(-ry_h, _SAS * _SDS, rz_h * _CDS))
            dzp = muladd( rx_h, _CAS * _CDS, muladd( ry_h, _SAS * _CDS, rz_h * _SDS))

            costheta2 = dzp * dzp * inv_r2
            J2_prefac = 3.0 * _J2_SUN * _R_SUN^2 / (r_sq * r_sq * r * 2.0)
            J2_fac    = muladd(5.0, costheta2, -1.0)

            resx = mu_i * J2_prefac * J2_fac * dxp
            resy = mu_i * J2_prefac * J2_fac * dyp
            resz = mu_i * J2_prefac * muladd(J2_fac, 1.0, -2.0) * dzp

            resxp = muladd(-resx, _SAS, muladd(-resy, _CAS * _SDS, resz * _CAS * _CDS))
            resyp = muladd( resx, _CAS, muladd(-resy, _SAS * _SDS, resz * _SAS * _CDS))
            reszp = muladd(resy, _CDS, resz * _SDS)

            ax += resxp
            ay += resyp
            az += reszp

        # =======================================================
        # LURRA (NAIF 3, EMB zuzenduta): J2, J3, J4
        # =======================================================
        elseif i == idx_emb
            # Lurraren poloa: RAe=0°, Dece=90° (J2000 ekuatorian)
            # Biraketa sinplifikatua: dxp_e=dy, dyp_e=-dx, dzp_e=dz
            dxp_e = -dy
            dyp_e =  dx
            dzp_e = -dz

            costheta2_e = dzp_e * dzp_e * inv_r2

            # --- Earth J2 ---
            J2_prefac_e = 3.0 * _J2_EARTH * _R_EARTH^2 / (r_sq * r_sq * r * 2.0)
            J2_fac_e    = muladd(5.0, costheta2_e, -1.0)

            resx_e = mu_i * J2_prefac_e * J2_fac_e * dxp_e
            resy_e = mu_i * J2_prefac_e * J2_fac_e * dyp_e
            resz_e = mu_i * J2_prefac_e * muladd(J2_fac_e, 1.0, -2.0) * dzp_e

            # --- Earth J3 ---
            J3_prefac_e = 5.0 * _J3_EARTH * _R_EARTH^3 / (r_sq * r_sq * r * 2.0)
            J3_fac_e    = muladd(-7.0, costheta2_e, 3.0)

            resx_e += -mu_i * J3_prefac_e * inv_r2 * J3_fac_e * dxp_e * dzp_e
            resy_e += -mu_i * J3_prefac_e * inv_r2 * J3_fac_e * dyp_e * dzp_e
            resz_e += -mu_i * J3_prefac_e * muladd(muladd(-7.0, costheta2_e, 6.0), costheta2_e, -0.6)

            # --- Earth J4 ---
            J4_prefac_e = 5.0 * _J4_EARTH * _R_EARTH^4 / (r_sq * r_sq * r_sq * r * 8.0)
            J4_fac_e    = muladd(muladd(63.0, costheta2_e, -42.0), costheta2_e, 3.0)

            resx_e += mu_i * J4_prefac_e * J4_fac_e * dxp_e
            resy_e += mu_i * J4_prefac_e * J4_fac_e * dyp_e
            resz_e += mu_i * J4_prefac_e * muladd(J4_fac_e + 12.0, 1.0, -28.0 * costheta2_e) * dzp_e

            # Itzulera sinplifikatua (RAe=0°, Dece=90°):
            ax -= resy_e
            ay += resx_e
            az += resz_e
        end
    end

    # 3. Deribatuak esleitu
    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] = ax
    du[5] = ay
    du[6] = az

    return nothing
end

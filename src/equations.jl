using StaticArrays # Gomendatua: svector-ak erabiltzeko

"""
f_all!(du, u, p, t)

Satelite baten mugimendu ekuazioak kalkulatzen ditu N-gorputzen eremuan.
- `u`: [x, y, z, vx, vy, vz]
- `p`: `mus` (grabitate parametroak) eta `bodies` (posizioak/funtzioak)
"""
function f_all!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko (destructuring)
    # @views erabiliz ez dugu memoria berririk erreserbatzen
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    ax, ay, az = 0.0, 0.0, 0.0

    # 2. Grabitatearen kalkulua (Newton hutsa)
    # zip + Tuple erabiltzen da: konpiladoreak luzera ezagutzen du,
    # loop-a unroll egin dezake, eta 0 heap alokazio sortzen dira.
    @inbounds for (mu_i, body_i) in zip(p.mus, p.bodies)
        bx, by, bz = body_i(t)

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
- `p`: `mus` (grabitate parametroak) eta `bodies` (posizioak/funtzioak)
"""
function f_all_rel!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko (destructuring)
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # Argiaren abiadura km/s-tan eta bere karratua
    c2 = 299792.458^2

    ax, ay, az = 0.0, 0.0, 0.0

    # 2. Grabitatearen kalkulua (Newtondarra + erlatibista Eguzkiarentzat)
    @inbounds for (i, (mu_i, body_i)) in enumerate(zip(p.mus, p.bodies))
        state_i = body_i(t)  # SVector{6}: [x, y, z, vx, vy, vz]
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
        # Bakarrik Eguzkiarentzat (i == 1)
        # =======================================================
        if i == 1
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
            r_h_norm = sqrt(rx_h^2 + ry_h^2 + rz_h^2)
            v_h_sq = vx_h^2 + vy_h^2 + vz_h^2
            r_dot_v = rx_h * vx_h + ry_h * vy_h + rz_h * vz_h

            # a_REL = (GM☉ / (c² |r|³)) * [(4GM☉/|r| - |v|²)r + 4(r·v)v]
            coef = mu_i / (c2 * r_h_norm^3)
            term1 = (4.0 * mu_i / r_h_norm) - v_h_sq

            ax += coef * (term1 * rx_h + 4.0 * r_dot_v * vx_h)
            ay += coef * (term1 * ry_h + 4.0 * r_dot_v * vy_h)
            az += coef * (term1 * rz_h + 4.0 * r_dot_v * vz_h)
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
- `p`: `mus` (grabitate parametroak) eta `bodies` (posizioak/funtzioak)
       Ordena: Eguzkia (1), Lurra (2), Ilargia (3), ...
"""
function f_all_rel_J!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # Konstanteak
    c2 = 299792.458^2           # Argiaren abiaduraren karratua (km/s)²

    # Eguzkiaren J2 (Park et al. 2021)
    J2_Sun = 2.2e-7
    R_Sun  = 695700.0           # km

    # Lurraren Jn koefizienteak (DE441 headerra, Park et al. 2021)
    J2_Earth = 1.08e-3
    J3_Earth = -2.53e-6
    J4_Earth = -1.62e-6
    R_Earth  = 6378.137         # km

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

    # 2. Grabitatearen kalkulua
    @inbounds for (i, (mu_i, body_i)) in enumerate(zip(p.mus, p.bodies))
        state_i = body_i(t)
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
        # EGUZKIA (i == 1): Erlatibista + J2
        # =======================================================
        if i == 1
            # --- Efektu erlatibista (Damour & Deruelle) ---
            bvx, bvy, bvz = state_i[4], state_i[5], state_i[6]
            rx_h = x - bx;  ry_h = y - by;  rz_h = z - bz
            vx_h = vx - bvx; vy_h = vy - bvy; vz_h = vz - bvz

            r_h_norm = sqrt(rx_h^2 + ry_h^2 + rz_h^2)
            v_h_sq = vx_h^2 + vy_h^2 + vz_h^2
            r_dot_v = rx_h * vx_h + ry_h * vy_h + rz_h * vz_h

            coef = mu_i / (c2 * r_h_norm^3)
            term1 = (4.0 * mu_i / r_h_norm) - v_h_sq

            ax += coef * (term1 * rx_h + 4.0 * r_dot_v * vx_h)
            ay += coef * (term1 * ry_h + 4.0 * r_dot_v * vy_h)
            az += coef * (term1 * rz_h + 4.0 * r_dot_v * vz_h)

            # --- Eguzkiaren J2 ---
            # J2_Sun = 2.2e-7 oso txikia da, eta ekliptikaren z erabiltzen
            # dugu hurbilketa gisa (Eguzkiaren poloa ~7.25° ekliptikatik)
            z2_r2 = dz^2 * inv_r2
            J2_coef = 1.5 * mu_i * J2_Sun * R_Sun^2 * inv_r3 * inv_r2
            fxy = J2_coef * (1.0 - 5.0 * z2_r2)
            fz  = J2_coef * (3.0 - 5.0 * z2_r2)
            ax += dx * fxy
            ay += dy * fxy
            az += dz * fz

        # =======================================================
        # LURRA (i == 2): J2, J3, J4 — EKUATORE ERREFERENTZIAN
        # =======================================================
        elseif i == 2
            # Ekliptikatik Lurraren ekuatorialera biratu: R_x(ε)
            # x ardatza partekatua da (udaberri puntua), beraz dx ez da aldatzen
            dx_eq = dx
            dy_eq = cosε * dy - sinε * dz
            dz_eq = sinε * dy + cosε * dz

            # r berdina da bi erreferentzia-markoetan (biraketa isometrikoa)
            dz_eq2_r2 = dz_eq^2 * inv_r2
            dz_eq4_r4 = dz_eq2_r2 * dz_eq2_r2

            # Harmonikoen azelerazioa Lurraren ekuatorialean kalkulatu
            ax_eq = 0.0
            ay_eq = 0.0
            az_eq = 0.0

            # --- Lurraren J2 ---
            J2_coef = 1.5 * mu_i * J2_Earth * R_Earth^2 * inv_r3 * inv_r2
            fxy_2 = J2_coef * (1.0 - 5.0 * dz_eq2_r2)
            fz_2  = J2_coef * (3.0 - 5.0 * dz_eq2_r2)
            ax_eq += dx_eq * fxy_2
            ay_eq += dy_eq * fxy_2
            az_eq += dz_eq * fz_2

            # --- Lurraren J3 ---
            inv_r7 = inv_r3 * inv_r2 * inv_r2
            J3_coef_xy = -2.5 * mu_i * J3_Earth * R_Earth^3 * inv_r7
            J3_fxy = J3_coef_xy * dz_eq * (3.0 - 7.0 * dz_eq2_r2)
            ax_eq += dx_eq * J3_fxy
            ay_eq += dy_eq * J3_fxy

            J3_coef_z = 0.5 * mu_i * J3_Earth * R_Earth^3 * inv_r3 * inv_r2
            az_eq += J3_coef_z * (3.0 - 30.0 * dz_eq2_r2 + 35.0 * dz_eq4_r4)

            # --- Lurraren J4 ---
            J4_coef_xy = -1.875 * mu_i * J4_Earth * R_Earth^4 * inv_r7  # 15/8 = 1.875
            J4_fxy = J4_coef_xy * (21.0 * dz_eq4_r4 - 14.0 * dz_eq2_r2 + 1.0)
            ax_eq += dx_eq * J4_fxy
            ay_eq += dy_eq * J4_fxy

            J4_coef_z = -0.625 * mu_i * J4_Earth * R_Earth^4 * inv_r7  # 5/8 = 0.625
            az_eq += dz_eq * J4_coef_z * (15.0 - 70.0 * dz_eq2_r2 + 63.0 * dz_eq4_r4)

            # Ekuatorialetik ekliptikara biratu: R_x(-ε) = R_x(ε)ᵀ
            ax +=  ax_eq
            ay +=  cosε * ay_eq + sinε * az_eq
            az += -sinε * ay_eq + cosε * az_eq
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


function f_BigFloat!(du, u, p, t)
    # 1. Egoera errazago irakurtzeko (destructuring)
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # BigFloat aritmetika: dena BigFloat
    ax, ay, az = BigFloat(0), BigFloat(0), BigFloat(0)

    # 2. Grabitatearen kalkulua
    @inbounds for (mu_i, body_i) in zip(p.mus, p.bodies)
        bx, by, bz = body_i(t)

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

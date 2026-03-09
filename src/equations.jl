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

    # 2. Grabitatearen kalkulua
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

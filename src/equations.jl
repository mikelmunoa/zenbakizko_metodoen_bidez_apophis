using LinearAlgebra
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

    # 2. Grabitatearen kalkulua (Loop optimizatua)
    @inbounds for i in 1:length(p.mus)
        # Gorputzaren posizioa lortu
        # Pentsatu: p.bodies[i] finko bat izan daiteke edo funtzio bat: p.bodies(t)[i]
        bx, by, bz = p.bodies[i](t)

        dx = bx - x
        dy = by - y
        dz = bz - z

        r_sq = dx^2 + dy^2 + dz^2

        # Singulartasunak saihesteko (softening), aukerakoa:
        # r_sq += 1e-12

        # Doitasun handiko eta abiadura optimizatuko grabitatea
        inv_r = 1.0 / sqrt(r_sq)
        common = p.mus[i] * (inv_r^3)

        ax += dx * common
        ay += dy * common
        az += az + (dz * common) # Kontuz: az metatu behar da
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

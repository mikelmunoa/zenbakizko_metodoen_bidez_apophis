@testset "RK4 — sistema sinplea" begin
    # Osziladore harmonikoa: d²x/dt² = -ω²x
    # Soluzio analitikoa: x(t) = cos(ωt), v(t) = -ω·sin(ωt)
    ω = 1.0
    u0 = [1.0, 0.0]  # Vector{Float64}: RK4-k mutable array behar du

    function harmonic_osc!(du, u, p, t)
        du[1] = u[2]
        du[2] = -p.omega^2 * u[1]
    end

    p = (omega = ω,)
    T = 2π   # periodo oso bat
    h = 1e-3

    tt, uu = RK4(u0, 0.0, T, h, p, harmonic_osc!)

    # Soluzio analitikoaren konpara: x(t)=cos(t), v(t)=-sin(t)
    # tt[end] ez da zehazki 2π (T/h erredondeatua), baina errore numerikoa txikia da
    @test uu[end][1] ≈ cos(tt[end])  atol=1e-9
    @test uu[end][2] ≈ -sin(tt[end]) atol=1e-9
end

@testset "RK4 — energia kontserbazioa" begin
    ω = 1.0
    u0 = [1.0, 0.0]

    function harmonic_osc2!(du, u, p, t)
        du[1] = u[2]
        du[2] = -p.omega^2 * u[1]
    end

    p = (omega = ω,)
    tt, uu = RK4(u0, 0.0, 10.0, 1e-3, p, harmonic_osc2!)

    # Energia = 0.5*(v² + ω²x²) = 0.5 hasieran
    E0 = 0.5 * (uu[1][2]^2 + ω^2 * uu[1][1]^2)
    Ef = 0.5 * (uu[end][2]^2 + ω^2 * uu[end][1]^2)

    @test abs(Ef - E0) / E0 < 1e-8
end

@testset "RK4 — irteeraren tamaina" begin
    u0 = [1.0, 0.0]

    function trivial!(du, u, p, t)
        du[1] = 0.0
        du[2] = 0.0
    end

    T, h, m = 10.0, 0.1, 2
    p_empty = NamedTuple()
    tt, uu = RK4(u0, 0.0, T, h, p_empty, trivial!, m)

    # m=2 batekin: N_steps = round(T / (h*m)) + 1 elementu
    expected = round(Int, T / (h * m)) + 1
    @test length(uu) == expected
    @test length(tt) == expected
end

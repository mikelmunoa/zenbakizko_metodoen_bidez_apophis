@testset "dd_perturbation" begin
    AU = 1.495978707e8  # km

    # Apophis 1 AU-ra +X ardatzean, abiadura Y ardatzean (zirkular hurbilketa)
    sun_state = SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    apo_state = SVector(AU, 0.0, 0.0, 0.0, 30.0, 0.0)
    mu_sun    = 1.32712440018e11  # km³/s²

    dd = dd_perturbation(apo_state, sun_state, mu_sun)

    @test length(dd) == 3

    # D-D << Newton (gu ~1/1000 Newton baino txikiagoa)
    newton_mag = mu_sun / AU^2
    dd_mag = norm(dd)
    @test dd_mag > 0
    @test dd_mag < newton_mag * 1e-3

    # Magnitude erreferentzia: D-D ~1e-13 km/s² 1 AU-tan
    @test 1e-15 < dd_mag < 1e-11

    # Y eta Z osagaiak zero (r·v=0 eta v=vy konfigurazioan)
    @test abs(dd[3]) < 1e-30

    # Zuzenketa positiboa X-ean (grabitazioa arintzen du — perihelio biraketa)
    @test dd[1] > 0
end

@testset "dd_perturbation — v=0 kasua" begin
    AU     = 1.495978707e8
    mu_sun = 1.32712440018e11

    sun_state = SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    apo_state = SVector(AU, 0.0, 0.0, 0.0, 0.0, 0.0)  # abiadura = 0

    dd = dd_perturbation(apo_state, sun_state, mu_sun)

    # r·v = 0, v² = 0 → Y eta Z = 0
    @test dd[2] ≈ 0.0 atol=1e-30
    @test dd[3] ≈ 0.0 atol=1e-30
    # 4*mu/r > 0 → dd[1] > 0
    @test dd[1] > 0
end

@testset "f_master! — oinarrizko deialdi" begin
    @test_nowarn dd_perturbation(
        SVector(1e8, 0.0, 0.0, 0.0, 20.0, 0.0),
        SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        1.32712440018e11
    )
end

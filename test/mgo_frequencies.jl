@testset "MgO frequency regression at selected q-points" begin
    basepath = joinpath(DATA_DIR, "MgO")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ifc2_path = joinpath(basepath, "infile.forceconstant")

    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    uc = CrystalStructure(ucposcar_path)
    @test ifc2.has_polar_data

    ew = LatticeDynamicsToolkit.precompute_ewald_parameters(ifc2, uc)

    qdir_gamma = SVector{3,Float64}(1.0, 0.0, 0.0)
    q_gamma = SVector{3,Float64}(0.0, 0.0, 0.0)
    q_x = SVector{3,Float64}(0.5, 0.0, 0.0)
    q_l = SVector{3,Float64}(0.5, 0.5, 0.5)
    q_w = SVector{3,Float64}(0.5, 0.25, 0.0)

    function freqs_thz_at(q_frac::SVector{3,Float64})
        D = dynmat_q(ifc2, uc, q_frac; include_polar=true, qdir_gamma=qdir_gamma, ewald=ew)
        vals, _ = get_modes(Hermitian(D), Val{LatticeDynamicsToolkit.is_gamma(q_frac)}())
        ω = sign.(vals) .* sqrt.(abs.(vals))
        return ω .* LatticeDynamicsToolkit.frequency_Hartree_to_THz
    end

    # Reference values from current Julia MgO pipeline.
    # These act as regression guards while refactoring long-range code.
    expected_gamma = [
        0.0,
        0.0,
        0.0,
        1.3073669962668042e+01,
        1.3073669962668042e+01,
        2.1622190121976370e+01,
    ]
    expected_x = [
        8.7203751517091530e+00,
        8.7203751517091580e+00,
        1.0942832616142928e+01,
        1.0942832616142928e+01,
        1.6504173376748990e+01,
        1.6599101628140750e+01,
    ]
    expected_l = [
        8.7203751517091530e+00,
        8.7203751517091580e+00,
        1.0942832616142919e+01,
        1.0942832616142926e+01,
        1.6504173376748984e+01,
        1.6599101628140748e+01,
    ]
    expected_w = [
        8.4456354959476000e+00,
        1.0074514271461311e+01,
        1.2495247174084655e+01,
        1.3349031812068871e+01,
        1.3426262225798522e+01,
        1.6459296696135173e+01,
    ]

    γ = freqs_thz_at(q_gamma)
    x = freqs_thz_at(q_x)
    l = freqs_thz_at(q_l)
    w = freqs_thz_at(q_w)

    # Keep tolerance modest to catch regressions while allowing tiny numeric drift.
    @test γ ≈ expected_gamma atol=3e-3 rtol=0
    @test x ≈ expected_x atol=3e-3 rtol=0
    @test l ≈ expected_l atol=3e-3 rtol=0
    @test w ≈ expected_w atol=3e-3 rtol=0
end


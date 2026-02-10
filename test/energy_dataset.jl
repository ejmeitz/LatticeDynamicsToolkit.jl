@testset "Energy Dataset" begin
    basepath = joinpath(DATA_DIR, "SW")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ssposcar_path = joinpath(basepath, "infile.ssposcar")

    Random.seed!(1234)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    n_configs = 50_000
    temperature = 1300.0
    β = 1 / (LatticeDynamicsToolkit.kB_eV * temperature)
    settings = ClassicalConfigSettings(n_configs, temperature)
    ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, temperature, 3)


    tep_energies = make_energy_dataset(
        settings,
        uc,
        sc;
        ifc2 = ifc2,
        ifc3 = ifc3,
        ifc4 = ifc4,
        verbose = false
    )

    @test length(tep_energies) == n_configs
    @test all(all.(isfinite, tep_energies))

    V2 = getindex.(tep_energies, 1)
    V3 = getindex.(tep_energies, 2)
    V4 = getindex.(tep_energies, 3)

    # By Equipartition Theorem
    @test mean(V2) ≈ (1.5 * (length(sc) - 1) * LatticeDynamicsToolkit.kB_eV * temperature) atol=1e-1
    # Odd function integraded on symmetric interval
    @test mean(V3) ≈ 0.0 atol=1e-3
    # Results from cumulant expansion
    @test mean(V4) / length(sc) ≈  0.0021535 atol=1e-3
    @test (-0.5*β*mean(V3.^2) / length(sc)) ≈ -0.0016334 atol=1e-3

end


@testset "TEP Energy Consistency" begin
    # IFC2 and AmorphousIFC2 should give identical energies for the same displacements
    basepath = joinpath(DATA_DIR, "SW")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ssposcar_path = joinpath(basepath, "infile.ssposcar")

    ifc2 = read_ifc2(joinpath(basepath, "1300K_3UC", "infile.forceconstant"), ucposcar_path)
    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    ifc2_sc = remap(sc, uc, ifc2)[1]

    dynmat = dynmat_gamma(ifc2_sc, sc)
    freqs_sq, phi = get_modes(dynmat, Val{true}())
    freqs = sqrt.(freqs_sq)

    settings = ConfigSettings(4, 1300.0, Classical)
    u = canonical_configs(settings, freqs, phi, sc.m; verbose = false)

    e2_true = map(eachcol(u)) do u_flat
        u_svec = reinterpret(SVector{3, Float64}, u_flat)
        energies(u_svec, ifc2_sc)[1]
    end

    dense_ifc2 = DenseIFC2(ifc2_sc, sc)
    amorphous_ifc2 = AmorphousIFC2(dense_ifc2)

    e2_test = map(eachcol(u)) do u_flat
        u_svec = reinterpret(SVector{3, Float64}, u_flat)
        energies(u_svec, amorphous_ifc2)[1]
    end

    @test e2_true ≈ e2_test
end
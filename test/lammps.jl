# LAMMPS-dependent tests: run only when LAMMPSExt is loaded

if has_lammps()
    @testset "Energy Calculator (LAMMPS)" begin
        basepath = joinpath(DATA_DIR, "SW")
        ucposcar_path = joinpath(basepath, "infile.ucposcar")
        ssposcar_path = joinpath(basepath, "infile.ssposcar")

        ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, 1300.0, 3)
        uc = CrystalStructure(ucposcar_path)
        sc = CrystalStructure(ssposcar_path)

        n_configs = 100
        settings = ClassicalConfigSettings(n_configs, 1300.0)
        sw_pot = joinpath(basepath, "Si.sw")
        pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]
        make_calc = sc -> LAMMPSCalculator(sc, pot_cmds)

        tep_energies, V = make_energy_dataset(
            settings, uc, sc, make_calc;
            ifc2 = ifc2, ifc3 = ifc3, ifc4 = ifc4,
            verbose = false
        )

        # TEP energies are tested elsewhere
        # This is not testing correctness...fix that
        @test length(V) == n_configs
        @test all(isfinite, V)
        
    end

    # @testset "TI SW (LAMMPS)" begin
    #     basepath = joinpath(DATA_DIR, "SW")
    #     ucposcar_path = joinpath(basepath, "infile.ucposcar")
    #     ssposcar_path = joinpath(basepath, "infile.ssposcar")

    #     ifc2, _, _ = load_ifcs(basepath, ucposcar_path, 1300.0, 3)
    #     uc = CrystalStructure(ucposcar_path)
    #     sc = CrystalStructure(ssposcar_path)

    #     sw_pot = joinpath(basepath, "Si.sw")
    #     pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]
    #     settings = TISettings(1300.0, pot_cmds, 1000, 500; n_lambda = 7)

    #     F = sum(ThermodynmicIntegration(ifc2, sc, uc, settings))
    #     @test F ≈ -4.6745638 atol=1e-4
    # end

    # @testset "TI LJ (LAMMPS)" begin
    #     basepath = joinpath(DATA_DIR, "LJ")
    #     ucposcar_path = joinpath(basepath, "infile.ucposcar")
    #     ssposcar_path = joinpath(basepath, "infile.ssposcar")

    #     ifc2, _, _ = load_ifcs(basepath, ucposcar_path, 80.0, 4)
    #     uc = CrystalStructure(ucposcar_path)
    #     sc = CrystalStructure(ssposcar_path)

    #     pot_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]
    #     settings = TISettings(80.0, pot_cmds, 2000, 1000; n_lambda = 9)

    #     F = sum(ThermodynmicIntegration(ifc2, sc, uc, settings))
    #     @test F ≈ -0.0820711 atol=1e-5
    # end
else
    @testset "LAMMPS tests skipped (LAMMPS not loaded)" begin
        @test_skip "LAMMPS extension not available"
    end
end

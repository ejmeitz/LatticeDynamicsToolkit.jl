@testset "IO" begin
    basepath = joinpath(DATA_DIR, "SW")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")

    # POSCAR parsing
    _, x_frac_uc_file, L_uc_file = LatticeDynamicsToolkit.read_poscar_data(ucposcar_path)

    a = 5.43
    x_frac_uc = [SVector(0.0, 0.0, 0.0), SVector(0.25, 0.25, 0.25)]
    L_uc = a .* [
        0.500000000000  0.500000000000 0.000000000000;
        0.000000000000  0.500000000000 0.500000000000;
        0.500000000000  0.000000000000 0.500000000000
    ]
    L_uc = collect(transpose(L_uc))

    @test x_frac_uc ≈ x_frac_uc_file
    @test L_uc_file ≈ L_uc

    # Force constant reading
    ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, 1300.0, 3)
    @test ifc2.na == 2
    @test ifc3.na == 2
    @test ifc4.na == 2
end

@testset "Remap" begin
    basepath = joinpath(DATA_DIR, "SW")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ssposcar_path = joinpath(basepath, "1300K_5UC_remapped", "infile.ssposcar")

    ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, 1300.0, 3)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    new_ifc2, new_ifc3, new_ifc4 = remap(sc, uc, ifc2, ifc3, ifc4)

    @test new_ifc2.na == length(sc)
    @test new_ifc3.na == length(sc)
    @test new_ifc4.na == length(sc)

    # Spot-check energies after remap
    u = [@SVector rand(3) for i in 1:new_ifc2.na]
    e2, e3, e4 = energies(u, new_ifc2; fc3 = new_ifc3, fc4 = new_ifc4)
    @test e2 > 0
    @test isfinite(e2)

    # Compare remapped IFC2 energies to TDEP-remapped file (when available)
    remapped_path = joinpath(basepath, "1300K_5UC_remapped", "outfile.forceconstant_remapped")
    if isfile(remapped_path)
        tdep_ifc2_remapped = read_ifc2(remapped_path, ssposcar_path)
        e2_tdep, _, _ = energies(u, tdep_ifc2_remapped)
        @test e2 ≈ e2_tdep rtol = 1e-10
    end
end
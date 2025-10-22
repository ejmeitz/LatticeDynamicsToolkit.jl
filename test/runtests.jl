using Test
using TDEP_IFCs


data_dir = abspath(joinpath(@__DIR__, "..", "data"))

function load_sw_ifcs(ucposcar_path)
    ifc2_path = joinpath(data_dir, "100K_3UC", "infile.forceconstant")
    ifc3_path = joinpath(data_dir, "100K_3UC", "infile.forceconstant_thirdorder")
    ifc4_path = joinpath(data_dir, "100K_3UC", "infile.forceconstant_fourthorder")

    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    ifc3 = read_ifc3(ifc3_path, ucposcar_path)
    ifc4 = read_ifc4(ifc4_path, ucposcar_path)
    
    return ifc2, ifc3, ifc4
end


@testset "IO" begin 

    ucposcar_path = joinpath(data_dir, "infile.ucposcar")
    x_frac_uc_file, L_uc_file = read_poscar_data(ucposcar_path)

    a = 5.43
    x_frac_uc = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]
    L_uc = a .* [
        0.500000000000  0.500000000000 0.000000000000;
        0.000000000000  0.500000000000 0.500000000000;
        0.500000000000  0.000000000000 0.500000000000
        ]
    L_uc = collect(transpose(L_uc))

    @test x_frac_uc ≈ x_frac_uc_file
    @test L_uc_file ≈ L_uc

    ifc2, ifc3, ifc4 = load_sw_ifcs(ucposcar_path)

end

@testset "Remap" begin
    
    ucposcar_path = joinpath(data_dir, "infile.ucposcar")
    ssposcar_path = joinpath(data_dir, "100K_5UC_remapped", "infile.ssposcar")

    ifc2, ifc3, ifc4 = load_sw_ifcs(ucposcar_path)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    new_ifc2, new_ifc3, new_ifc4 = remap(sc, uc, ifc2, ifc3, ifc4)

    u = [@SVector rand(3) for i in 1:new_ifc2.na]
    println(TDEPToolkit.energies(u, new_ifc2; fc3 = new_ifc3))

    ifc2_remapped_path = joinpath(data_dir, "100K_5UC_remapped", "outfile.forceconstant_remapped")
    ifc3_remapped_path = joinpath(data_dir, "100K_5UC_remapped", "outfile.forceconstant_thirdorder_remapped")

    tdep_ifc2_remapped = read_ifc2(ifc2_remapped_path, ssposcar_path)
    tdep_ifc3_remapped = read_ifc3(ifc3_remapped_path, ssposcar_path)

    println(TDEPToolkit.energies_faithful(u, tdep_ifc2_remapped; fc3 = ifc3_remapped_path))

end


@testset "Energy Dataset" begin
    
    data_dir = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data"
    ucposcar_path = joinpath(data_dir, "infile.ucposcar")
    ssposcar_path = joinpath(data_dir, "100K_3UC", "infile.ssposcar")

    ifc2, ifc3, ifc4 = load_sw_ifcs(ucposcar_path)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    n_configs = 10_000
    temperature = 1300.0
    settings = ClassicalConfigSettings(n_configs, temperature)

    dynmat = TDEPToolkit.dynmat_gamma(ifc2_remapped, sc)
    freqs_sq, phi = TDEPToolkit.get_modes(dynmat)
    freqs_Thz = sqrt.(freqs_sq) .* TDEPToolkit.frequency_Hartree_to_THz

    tep_energies = TDEPToolkit.make_energy_dataset(
        settings,
        uc,
        sc;
        ifc2 = ifc2,
        ifc3 = ifc3,
        ifc4 = ifc4
    )

end
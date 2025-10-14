using Test
using TDEP_IFCs


data_dir = abspath(joinpath(@__DIR__, "..", "data"))

function load_sw_ifcs(ucposcar_path)
    ifc2_path = joinpath(data_dir, "infile.forceconstant")
    ifc3_path = joinpath(data_dir, "infile.forceconstant_thirdorder")
    ifc4_path = joinpath(data_dir, "infile.forceconstant_fourthorder")

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
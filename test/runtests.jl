using Test
using TDEP_IFCs
using LinearALgebra

data_dir = abspath(joinpath(@__DIR__, "..", "data"))

function load_ifcs(basepath, ucposcar_path, T, N::Integer)
    ifc2_path = joinpath(basepath, "$(Int(T))K_$(N)UC", "infile.forceconstant")
    ifc3_path = joinpath(basepath, "$(Int(T))K_$(N)UC", "infile.forceconstant_thirdorder")
    ifc4_path = joinpath(basepath, "$(Int(T))K_$(N)UC", "infile.forceconstant_fourthorder")

    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    ifc3 = read_ifc3(ifc3_path, ucposcar_path)
    ifc4 = read_ifc4(ifc4_path, ucposcar_path)
    
    return ifc2, ifc3, ifc4
end

@testset "IO" begin 

    basepath = joinpath(data_dir, "SW")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
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

    ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, 1300.0, 3)

end

@testset "Remap" begin
    
    basepath = joinpath(data_dir, "SW")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ssposcar_path = joinpath(basepath, "100K_5UC_remapped", "infile.ssposcar")

    ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, 1300.0, 3)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    new_ifc2, new_ifc3, new_ifc4 = remap(sc, uc, ifc2, ifc3, ifc4)

    u = [@SVector rand(3) for i in 1:new_ifc2.na]
    println(TDEPToolkit.energies(u, new_ifc2; fc3 = new_ifc3))

    ifc2_remapped_path = joinpath(basepath, "100K_5UC_remapped", "outfile.forceconstant_remapped")
    ifc3_remapped_path = joinpath(basepath, "100K_5UC_remapped", "outfile.forceconstant_thirdorder_remapped")

    tdep_ifc2_remapped = read_ifc2(ifc2_remapped_path, ssposcar_path)
    tdep_ifc3_remapped = read_ifc3(ifc3_remapped_path, ssposcar_path)

    println(TDEPToolkit.energies_faithful(u, tdep_ifc2_remapped; fc3 = ifc3_remapped_path))

end


@testset "Energy Dataset" begin
    
    basepath = joinpath(data_dir, "SW")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ssposcar_path = joinpath(basepath, "100K_3UC", "infile.ssposcar")

    ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, 1300.0, 3)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    n_configs = 10_000
    temperature = 1300.0
    settings = ClassicalConfigSettings(n_configs, temperature)

    tep_energies = make_energy_dataset(
        settings,
        uc,
        sc;
        ifc2 = ifc2,
        ifc3 = ifc3,
        ifc4 = ifc4
    )

end

#TODO TEST WITH THE 8 ATOM UNIT CELL TO VERIFY I GET SAME RESULTS
@testset "Energy Calculator" begin
    
    basepath = joinpath(data_dir, "SW")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ssposcar_path = joinpath(basepath, "infile.ssposcar")

    T = 1300.0
    ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, 1300.0, 3)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    n_configs = 100_000
    settings = ClassicalConfigSettings(n_configs, T)

    sw_pot = joinpath(basepath, "Si.sw")
    pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]

    make_calc = (sc) -> LAMMPSCalculator(sc, pot_cmds)

    tep_energies, V = make_energy_dataset(
        settings,
        uc,
        sc,
        make_calc;
        ifc2 = ifc2,
        ifc3 = ifc3,
        ifc4 = ifc4,
    )

end

@testset "TI SW" begin
    
    # basepath = joinpath(data_dir, "SW)")
    basepath = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data\SW"
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ssposcar_path = joinpath(basepath, "infile.ssposcar")

    T = 1300.0
    ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, 1300.0, 3)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    sw_pot = joinpath(basepath, "Si.sw")
    pot_cmds = ["pair_style sw", "pair_coeff * * \"$(sw_pot)\" Si"]

    nsteps = 150_000
    nsteps_equil = 50_000
    n_lambda = 9
    settings = TISettings(T, pot_cmds, nsteps, nsteps_equil; n_lambda = n_lambda)

    F = ThermodynmicIntegration(
        ifc2, sc, uc, settings
    )

end

@testset "TI LJ" begin
    
    # basepath = joinpath(data_dir, "LJ")
    basepath = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data\LJ"
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ssposcar_path = joinpath(basepath, "infile.ssposcar")

    T = 80.0
    ifc2, ifc3, ifc4 = load_ifcs(basepath, ucposcar_path, T, 4)

    uc = CrystalStructure(ucposcar_path)
    sc = CrystalStructure(ssposcar_path)

    pot_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.010423 3.4", "pair_modify shift yes"]

    nsteps = 150_000
    nsteps_equil = 50_000
    n_lambda = 9
    settings = TISettings(T, pot_cmds, nsteps, nsteps_equil; n_lambda = n_lambda)

    F = ThermodynmicIntegration(
        ifc2, sc, uc, settings
    )

end

@testset "Fit Amorphous" begin


    # ============ FILL THESE IN ============
    poscar_path = "C:/Users/ejmei/repos/TDEP_IFCs.jl/data/aSi/infile.ssposcar"
    # poscar_path = "C:/Users/ejmei/repos/TDEP_IFCs.jl/data/SW/infile.ssposcar"
    positions_path = "C:/Users/ejmei/repos/TDEP_IFCs.jl/data/aSi/infile.positions"
    # positions_path = "C:/Users/ejmei/repos/TDEP_IFCs.jl/data/SW/1300K_3UC/infile.positions"
    forces_path = "C:/Users/ejmei/repos/TDEP_IFCs.jl/data/aSi/infile.forces"
    # forces_path = "C:/Users/ejmei/repos/TDEP_IFCs.jl/data/SW/1300K_3UC/infile.forces"
    n_timesteps = 300 #300  # number of MD snapshots
    r_cut_angstrom = 6.0  # cutoff in Angstrom
    # =======================================

    # Convert cutoff to bohr (internal units)
    r_cut = r_cut_angstrom * LatticeDynamicsToolkit.A_to_bohr

    # Load structure (skip symmetry analysis for amorphous)
    println("Loading structure...")
    crystal = CrystalStructure(poscar_path; compute_symmetry=false)
    na = length(crystal)
    println("  Atoms: $na")

    # Read MD data
    println("Reading positions and forces...")
    positions = read_positions(positions_path, crystal, n_timesteps)
    forces = read_forces(forces_path, na, n_timesteps)
    println("  Positions shape: $(size(positions))")
    println("  Forces shape: $(size(forces))")

    # Fit IFCs
    println("\nFitting IFCs...")
    ifc = fit_ifc2(crystal, positions, forces; r_cut=r_cut, λ=1e-6)
    println(ifc)

    # Convert to dense matrix
    println("\nConverting to dense matrix...")
    Φ = Matrix(ifc)
    println("  Matrix size: $(size(Φ))")

    # ============ TESTS ============

    # Test 1: Check symmetry (Φ should equal Φ^T)
    symmetry_error = maximum(abs.(Φ - Φ'))
    println("1. Symmetry check: max|Φ - Φ^T| = $(symmetry_error)")
    println("   ", symmetry_error < 1e-10 ? "✓ PASSED" : "✗ FAILED")

    # Test 2: Check ASR (row sums should be zero)
    # For each atom i, sum over all j: Σ_j Φ_ij = 0
    asr_errors = Float64[]
    for i in 1:na
        for α in 1:3
            row = 3*(i-1) + α
            row_sum = sum(Φ[row, :])
            push!(asr_errors, abs(row_sum))
        end
    end
    max_asr_error = maximum(asr_errors)
    println("2. ASR check: max|Σ_j Φ_ij| = $(max_asr_error)")
    println("   ", max_asr_error < 1e-10 ? "✓ PASSED" : "✗ FAILED")

    # Test 3: Check block-wise symmetry Φ_ij = Φ_ji^T (Newton's 3rd law)
    block_sym_errors = Float64[]
    for i in 1:na
        ri = 3*(i-1)
        for j in i+1:na
            rj = 3*(j-1)
            Φ_ij = Φ[ri+1:ri+3, rj+1:rj+3]
            Φ_ji = Φ[rj+1:rj+3, ri+1:ri+3]
            push!(block_sym_errors, maximum(abs.(Φ_ij - Φ_ji')))
        end
    end
    if !isempty(block_sym_errors)
        max_block_sym = maximum(block_sym_errors)
        println("3. Newton's 3rd law: max|Φ_ij - Φ_ji^T| = $(max_block_sym)")
        println("   ", max_block_sym < 1e-10 ? "✓ PASSED" : "✗ FAILED")
    end

    # Test 4: Eigenvalue check (should have 3 zero modes for translations)
    println("\n4. Eigenvalue analysis...")
    eigenvalues, _ = eigen(Symmetric(Φ))
    sorted_eigs = sort(eigenvalues)
    println("   3 smallest eigenvalues: $(sorted_eigs[1:3])")
    println("   3 largest eigenvalues: $(sorted_eigs[end-2:end])")
    n_negative = count(x -> x < -1e-8, eigenvalues)
    println("   Negative eigenvalues: $n_negative")
end

@testset "TEP Energy Consistency" begin

    # Load crystalline silicion IFCs and calculate energy of equilibrium structure
    # Convert those crystalline IFCs to AmorphousIFC2 and calcaulte same energy

    ifc_path = "C:/Users/ejmei/repos/TDEP_IFCs.jl/data/SW/1300K_3UC/infile.forceconstant"
    ssposcar_path = "C:/Users/ejmei/repos/TDEP_IFCs.jl/data/SW/infile.ssposcar"
    ucposcar_path = "C:/Users/ejmei/repos/TDEP_IFCs.jl/data/SW/infile.ucposcar"

    ifc2 = read_ifc2(ifc_path, ucposcar_path)
    ifc2_sc = remap(sc, uc, ifc2)[1]
    sc = CrystalStructure(ssposcar_path)
    uc = CrystalStructure(ucposcar_path)

    
    dynmat = dynmat_gamma(ifc2_sc, sc)
    freqs_sq, phi = get_modes(dynmat, Val{true}())
    freqs = sqrt.(freqs_sq)

    settings = ConfigSettings(4, 1300.0, Classical)
    u = canonical_configs(settings, freqs, phi, sc.m)

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


    # ifc2_remapped = remap(sc, uc, ifc2)[1]
    # dynmat = TDEPToolkit.dynmat_gamma(ifc2_remapped, sc)
    # freqs_sq, phi = TDEPToolkit.get_modes(dynmat)
    # freqs_Thz = sqrt.(freqs_sq) .* TDEPToolkit.frequency_Hartree_to_THz
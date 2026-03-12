
@testset "Long-range Ewald: MgO comparison with TDEP" begin
    
    # Load MgO data
    basepath = joinpath(DATA_DIR, "MgO")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ifc2_path = joinpath(basepath, "infile.forceconstant")
    
    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    uc = CrystalStructure(ucposcar_path)
    
    @test ifc2.has_polar_data
    
    # Set up Ewald parameters using the lambda from the file
    ew = LatticeDynamicsToolkit.EwaldParameters()
    LatticeDynamicsToolkit.set_ewald_parameters!(ew, uc, ifc2.polar.eps; strategy=3, lambda_forced=ifc2.polar.lambda)
    
    # Test q-points: Gamma, X, L, and a few others
    q_points = [
        ("Gamma", SVector{3,Float64}(0.0, 0.0, 0.0)),
        ("X", SVector{3,Float64}(0.5, 0.0, 0.0)),
        ("L", SVector{3,Float64}(0.5, 0.5, 0.5)),
        ("W", SVector{3,Float64}(0.5, 0.25, 0.0)),
        ("K", SVector{3,Float64}(0.375, 0.375, 0.0)),
    ]
    
    println("\n" * "="^80)
    println("MgO Long-range Dynamical Matrix Comparison")
    println("="^80)
    println("Unit cell: $(length(uc)) atoms")
    println("Dielectric tensor:")
    println(ifc2.polar.eps)
    println("Lambda (Ewald parameter): $(ifc2.polar.lambda)")
    println("Born charges:")
    for a in 1:length(uc)
        println("  Atom $a:")
        println(ifc2.polar.born_Z[:, :, a])
    end
    println()
    
    for (name, q_frac) in q_points
        D_lr = LatticeDynamicsToolkit.longrange_dynamical_matrix(ew, uc, q_frac, ifc2.polar.born_Z, ifc2.polar.eps)
        
        # Print key values for comparison
        println("q-point: $name = ($(q_frac[1]), $(q_frac[2]), $(q_frac[3]))")
        println("  D[1,1,1,1] (real, imag): ($(real(D_lr[1,1,1,1])), $(imag(D_lr[1,1,1,1])))")
        println("  D[1,1,2,2] (real, imag): ($(real(D_lr[1,1,2,2])), $(imag(D_lr[1,1,2,2])))")
        println("  D[1,2,1,2] (real, imag): ($(real(D_lr[1,2,1,2])), $(imag(D_lr[1,2,1,2])))")
        
        # Check Hermiticity
        max_herm_error = 0.0
        for a1 in 1:length(uc), a2 in 1:length(uc), i in 1:3, j in 1:3
            err = abs(D_lr[i, j, a1, a2] - conj(D_lr[j, i, a2, a1]))
            max_herm_error = max(max_herm_error, err)
        end
        println("  Max Hermiticity error: $max_herm_error")
        @test max_herm_error < 1e-10
        
        # Print full matrix for Gamma (most important)
        if name == "Gamma"
            println("  Full D matrix (real part):")
            for a1 in 1:length(uc), a2 in 1:length(uc)
                println("    Block [a1=$a1, a2=$a2]:")
                println(real(D_lr[:, :, a1, a2]))
            end
        end
        println()
    end
    
    # Test with reconly=true (TDEP's default)
    println("Testing with reconly=true (TDEP default):")
    q_test = SVector{3,Float64}(0.1, 0.0, 0.0)
    D_full = LatticeDynamicsToolkit.longrange_dynamical_matrix(ew, uc, q_test, ifc2.polar.born_Z, ifc2.polar.eps; reconly=false)
    D_rec = LatticeDynamicsToolkit.longrange_dynamical_matrix(ew, uc, q_test, ifc2.polar.born_Z, ifc2.polar.eps; reconly=true)
    println("  D[1,1,1,1] full: $(D_full[1,1,1,1])")
    println("  D[1,1,1,1] rec:  $(D_rec[1,1,1,1])")
    println("  Difference: $(D_full[1,1,1,1] - D_rec[1,1,1,1])")
    println()
end

@testset "Long-range Ewald precompute parity in harmonic workflows" begin
    basepath = joinpath(DATA_DIR, "MgO")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ifc2_path = joinpath(basepath, "infile.forceconstant")

    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    uc = CrystalStructure(ucposcar_path)
    @test ifc2.has_polar_data

    ew = LatticeDynamicsToolkit.EwaldParameters()
    LatticeDynamicsToolkit.set_ewald_parameters!(ew, uc, ifc2.polar.eps; strategy=3, lambda_forced=ifc2.polar.lambda)

    q = SVector{3,Float64}(0.375, 0.375, 0.0)

    D_on_demand = dynmat_q(ifc2, uc, q; ewald=nothing)
    D_precomputed = dynmat_q(ifc2, uc, q; ewald=ew)
    @test D_on_demand ≈ D_precomputed atol=1e-10 rtol=1e-10

    Dd_on_demand, dDdq_on_demand = dynmat_and_derivative_q(ifc2, uc, q; ewald=nothing)
    Dd_precomputed, dDdq_precomputed = dynmat_and_derivative_q(ifc2, uc, q; ewald=ew)
    @test Dd_on_demand ≈ Dd_precomputed atol=1e-10 rtol=1e-10
    @test dDdq_on_demand ≈ dDdq_precomputed atol=1e-9 rtol=1e-9

    mesh = SVector{3,Int}(2, 2, 2)
    nb = 3 * length(uc)
    ibz = SimpleMesh(uc, mesh)
    dd = LatticeDynamicsToolkit.DispersionDataSimple(uc, ifc2, ibz; n_threads=1)

    q1 = ibz.k_ibz[1]
    D1_ref, dD1_ref = dynmat_and_derivative_q(ifc2, uc, q1; ewald=nothing)
    freqs_sq_ref, phi_ref = get_modes(Hermitian(D1_ref), Val{LatticeDynamicsToolkit.is_gamma(q1)}())
    freqs_ref = LatticeDynamicsToolkit.negsqrt.(freqs_sq_ref)
    ix = sortperm(freqs_ref)
    freqs_ref_sorted = freqs_ref[ix]
    phi_ref_sorted = phi_ref[:, ix]
    vels_ref = LatticeDynamicsToolkit.group_velocities(freqs_ref_sorted, phi_ref_sorted, dD1_ref)
    @test dd.freqs[1] ≈ SVector{nb, Float64}(freqs_ref_sorted) atol=1e-10 rtol=1e-10
    @test dd.vels[1] ≈ SMatrix{3, nb, Float64}(vels_ref) atol=1e-9 rtol=1e-9

    f = ω -> ω
    thermo_val = LatticeDynamicsToolkit.sum_over_freqs(f, ibz, uc, ifc2; n_threads=1)
    thermo_ref = 0.0
    for i in eachindex(ibz.weights)
        q_i = ibz.k_ibz[i]
        D_i = dynmat_q(ifc2, uc, q_i; ewald=nothing)
        freqs_sq_i, _ = get_modes(Hermitian(D_i), Val{LatticeDynamicsToolkit.is_gamma(q_i)}())
        freqs_i = LatticeDynamicsToolkit.negsqrt.(freqs_sq_i)
        local_sum = 0.0
        for ω in freqs_i
            if ω > LatticeDynamicsToolkit.lo_freqtol
                local_sum += f(ω)
            end
        end
        thermo_ref += ibz.weights[i] * local_sum
    end
    @test thermo_val ≈ thermo_ref atol=1e-10 rtol=1e-10
end

@testset "Ewald strategy coverage and supercell polar energy consistency" begin
    basepath = joinpath(DATA_DIR, "MgO")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ifc2_path = joinpath(basepath, "infile.forceconstant")

    ifc2_uc = read_ifc2(ifc2_path, ucposcar_path)
    uc = CrystalStructure(ucposcar_path)
    sc = uc
    ifc2_sc = ifc2_uc

    @test ifc2_uc.has_polar_data
    @test ifc2_sc.has_polar_data

    ew1 = LatticeDynamicsToolkit.EwaldParameters()
    ew2 = LatticeDynamicsToolkit.EwaldParameters()
    ew3 = LatticeDynamicsToolkit.EwaldParameters()
    # Explicitly exercise strategy=1 and strategy=2 paths (strategy=3 is current production path).
    LatticeDynamicsToolkit.set_ewald_parameters!(ew1, uc, ifc2_uc.polar.eps; strategy=1)
    LatticeDynamicsToolkit.set_ewald_parameters!(ew2, uc, ifc2_uc.polar.eps; strategy=2)
    LatticeDynamicsToolkit.set_ewald_parameters!(ew3, uc, ifc2_uc.polar.eps; strategy=3, lambda_forced=ifc2_uc.polar.lambda)

    q = SVector{3,Float64}(0.5, 0.25, 0.0)
    D1 = LatticeDynamicsToolkit.longrange_dynamical_matrix(ew1, uc, q, ifc2_uc.polar.born_Z, ifc2_uc.polar.eps; reconly=true)
    D2 = LatticeDynamicsToolkit.longrange_dynamical_matrix(ew2, uc, q, ifc2_uc.polar.born_Z, ifc2_uc.polar.eps; reconly=true)
    D3 = LatticeDynamicsToolkit.longrange_dynamical_matrix(ew3, uc, q, ifc2_uc.polar.born_Z, ifc2_uc.polar.eps; reconly=true)

    # Ewald splitting parameter changes decomposition cost, but physical D(q) should be consistent.
    @test D1 ≈ D3 atol=1e-7 rtol=1e-6
    @test D2 ≈ D3 atol=1e-7 rtol=1e-6

    Φlr1 = LatticeDynamicsToolkit.supercell_longrange_forceconstant(ew1, ifc2_uc.polar.born_Z, ifc2_uc.polar.eps, sc, uc)
    Φlr2 = LatticeDynamicsToolkit.supercell_longrange_forceconstant(ew2, ifc2_uc.polar.born_Z, ifc2_uc.polar.eps, sc, uc)
    Φlr3 = LatticeDynamicsToolkit.supercell_longrange_forceconstant(ew3, ifc2_uc.polar.born_Z, ifc2_uc.polar.eps, sc, uc)

    Random.seed!(42)
    u = [SVector{3,Float64}(1e-3 * randn(), 1e-3 * randn(), 1e-3 * randn()) for _ in 1:length(sc)]

    function ep_from_forceconstants(uvec::Vector{SVector{3,Float64}}, fcp::Array{Float64,4})
        acc = 0.0
        na = length(uvec)
        @inbounds for a1 in 1:na, a2 in 1:na
            acc += dot(uvec[a1], fcp[:, :, a1, a2] * uvec[a2])
        end
        return 0.5 * acc
    end

    ep1 = ep_from_forceconstants(u, Φlr1)
    ep2 = ep_from_forceconstants(u, Φlr2)
    ep3 = ep_from_forceconstants(u, Φlr3)

    dense = DensePolarIFCs(Φlr3, length(sc))
    _, _, _, ep_dense = energies(u, ifc2_sc; fcp=dense, n_threads=1)

    # Validate supercell LRFC -> polar energy path used in epot.
    @test ep_dense ≈ ep3 atol=1e-10 rtol=1e-8
    @test ep1 ≈ ep3 atol=1e-8 rtol=1e-4
    @test ep2 ≈ ep3 atol=1e-8 rtol=1e-4
end


@testset "Long-range Ewald: MgO comparison with TDEP" begin
    
    # Load MgO data
    basepath = joinpath(DATA_DIR, "MgO")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ifc2_path = joinpath(basepath, "infile.forceconstant")
    
    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    uc = CrystalStructure(ucposcar_path)
    
    @test ifc2.has_polar_data "MgO IFC2 should have polar data"
    
    # Set up Ewald parameters using the lambda from the file
    ew = EwaldParameters()
    set_ewald_parameters!(ew, uc, ifc2.polar.eps; strategy=3, lambda_forced=ifc2.polar.lambda)
    
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
        D_lr = longrange_dynamical_matrix(ew, uc, q_frac, ifc2.polar.born_Z, ifc2.polar.eps)
        
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
    D_full = longrange_dynamical_matrix(ew, uc, q_test, ifc2.polar.born_Z, ifc2.polar.eps; reconly=false)
    D_rec = longrange_dynamical_matrix(ew, uc, q_test, ifc2.polar.born_Z, ifc2.polar.eps; reconly=true)
    println("  D[1,1,1,1] full: $(D_full[1,1,1,1])")
    println("  D[1,1,1,1] rec:  $(D_rec[1,1,1,1])")
    println("  Difference: $(D_full[1,1,1,1] - D_rec[1,1,1,1])")
    println()
end

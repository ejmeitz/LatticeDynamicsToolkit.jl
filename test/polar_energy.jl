@testset "Polar Energy Consistency (MgO)" begin
    basepath = joinpath(DATA_DIR, "MgO")
    ucposcar_path = joinpath(basepath, "infile.ucposcar")
    ifc2_path = joinpath(basepath, "infile.forceconstant")

    ifc2 = read_ifc2(ifc2_path, ucposcar_path)
    uc = CrystalStructure(ucposcar_path)
    sc = uc

    @test ifc2.has_polar_data

    ew = LatticeDynamicsToolkit.EwaldParameters()
    LatticeDynamicsToolkit.set_ewald_parameters!(ew, uc, ifc2.polar.eps; strategy=3, lambda_forced=ifc2.polar.lambda)

    Φlr = LatticeDynamicsToolkit.supercell_longrange_forceconstant(ew, ifc2, sc, uc)
    dense = DensePolarIFCs(ifc2, uc, sc)

    Random.seed!(17)
    u = [SVector{3,Float64}(1e-3 * randn(), 1e-3 * randn(), 1e-3 * randn()) for _ in 1:length(sc)]

    ep_direct = 0.0
    @inbounds for a1 in 1:length(sc), a2 in 1:length(sc)
        ep_direct += dot(u[a1], Φlr[:, :, a1, a2] * u[a2])
    end
    ep_direct *= 0.5

    _, _, _, ep_dense = energies(u, ifc2; fcp=dense, n_threads=1)
    @test ep_dense ≈ ep_direct atol=1e-10 rtol=1e-8
end

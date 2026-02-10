# @testset "Fit Amorphous" begin
#     basepath = joinpath(DATA_DIR, "aSi")
#     poscar_path = joinpath(basepath, "infile.ssposcar")
#     positions_path = joinpath(basepath, "infile.positions")
#     forces_path = joinpath(basepath, "infile.forces")
#     n_timesteps = 300
#     r_cut_angstrom = 6.0

#     r_cut = r_cut_angstrom * LatticeDynamicsToolkit.A_to_bohr

#     crystal = CrystalStructure(poscar_path; compute_symmetry=false)
#     na = length(crystal)

#     positions = read_positions(positions_path, crystal, n_timesteps)
#     forces = read_forces(forces_path, na, n_timesteps)

#     ifc, _ = fit_ifc2(crystal, positions, forces, r_cut; λ=1e-6)

#     @test ifc.na == na

#     Φ = Matrix(ifc)

#     # Symmetry: Φ = Φ'
#     @test maximum(abs.(Φ - Φ')) < 1e-8

#     # ASR: row sums zero (relaxed tol - fit uses penalty constraints)
#     for i in 1:na
#         for α in 1:3
#             row = 3*(i-1) + α
#             @test abs(sum(Φ[row, :])) < 5e-5
#         end
#     end

#     # Block symmetry Φ_ij = Φ_ji^T
#     for i in 1:na
#         ri = 3*(i-1)
#         for j in (i+1):na
#             rj = 3*(j-1)
#             Φ_ij = Φ[ri+1:ri+3, rj+1:rj+3]
#             Φ_ji = Φ[rj+1:rj+3, ri+1:ri+3]
#             @test maximum(abs.(Φ_ij - Φ_ji')) < 1e-8
#         end
#     end

#     # 3 zero modes (translations)
#     eigenvalues, _ = eigen(Symmetric(Φ))
#     n_negative = count(x -> x < -1e-8, eigenvalues)
#     @test n_negative == 0
# end

# #!/usr/bin/env julia

# using Pkg
# Pkg.activate(@__DIR__)

# using LinearAlgebra
# using Printf
# using StaticArrays
# using LatticeDynamicsToolkit
# import LatticeDynamicsToolkit: frequency_Hartree_to_THz

# const Q_GAMMA = SVector{3,Float64}(0.0, 0.0, 0.0)
# const Q_X = SVector{3,Float64}(0.5, 0.0, 0.0)
# const QDIR_GAMMA = SVector{3,Float64}(1.0, 0.0, 0.0)

# is_gamma_q(q::SVector{3,Float64}) = norm(q) < 1e-14

# function hermitian_residual(M::AbstractMatrix{ComplexF64})
#     return norm(M .- M')
# end

# function tensor4_herm_residual(T::Array{Float64,4})
#     na = size(T, 3)
#     acc = 0.0
#     for a1 in 1:na, a2 in 1:na, i in 1:3, j in 1:3
#         d = T[i, j, a1, a2] - T[j, i, a2, a1]
#         acc += d^2
#     end
#     return sqrt(acc)
# end

# function add_massweighted_lr!(D::Matrix{ComplexF64}, D_lr::Array{ComplexF64,4}, uc::CrystalStructure)
#     na = length(uc)
#     @inbounds for a2 in 1:na, a1 in 1:na
#         w = uc.invsqrtm[a1] * uc.invsqrtm[a2]
#         r1, r2 = 3*(a1-1)+1, 3*(a2-1)+1
#         for j in 1:3, i in 1:3
#             D[r1+i-1, r2+j-1] += D_lr[i, j, a1, a2] * w
#         end
#     end
#     return D
# end

# function born_onsite_correction(ew, ifc2::IFC2, uc::CrystalStructure)
#     na = length(uc)
#     q0 = SVector{3,Float64}(0.0, 0.0, 0.0)
#     Dg = LatticeDynamicsToolkit.longrange_dynamical_matrix(
#         ew, uc, q0, ifc2.polar.born_Z, ifc2.polar.eps; reconly=true
#     )
#     onsite = zeros(Float64, 3, 3, na)
#     @inbounds for a1 in 1:na
#         m = zeros(Float64, 3, 3)
#         for a2 in 1:na
#             m .+= real.(Dg[:, :, a1, a2])
#         end
#         onsite[:, :, a1] .= -m
#     end
#     return onsite
# end

# freqs_from_eigvals(vals::AbstractVector{<:Real}) = sign.(vals) .* sqrt.(abs.(vals))

# function write_metric(io::IO, category::AbstractString, qpoint::AbstractString, metric::AbstractString, value::Real)
#     println(io, string(category, ",", qpoint, ",", metric, ",", @sprintf("%.16e", value)))
# end

# function write_freqs(io::IO, qlabel::String, prefix::String, freqs_hartree::AbstractVector{<:Real})
#     for (i, ω) in enumerate(freqs_hartree)
#         write_metric(io, prefix, qlabel, "mode_$(i)", ω * frequency_Hartree_to_THz)
#     end
# end

# function write_matrix_norms_complex(io::IO, category::String, qlabel::String, M::AbstractMatrix{ComplexF64})
#     write_metric(io, category, qlabel, "trace_real", real(tr(M)))
#     write_metric(io, category, qlabel, "trace_imag", imag(tr(M)))
#     write_metric(io, category, qlabel, "frob", norm(M))
#     write_metric(io, category, qlabel, "herm_residual", hermitian_residual(M))
# end

# function write_ewald_probe(io::IO, ew, uc::CrystalStructure, ifc2::IFC2)
#     tol = 1e-20
#     λ = ifc2.polar.lambda
#     L_rec = inv(uc.L)'
#     pts = LatticeDynamicsToolkit._points_on_sphere(40)
#     kins = LatticeDynamicsToolkit.inscribed_sphere_in_box(L_rec)
#     k0 = kins * 1e-9
#     k1 = kins * 1e9
#     inv4 = 1.0 / (4.0 * λ * λ)

#     fk(r) = begin
#         s = 0.0
#         for i in 1:size(pts, 2)
#             v = SVector{3,Float64}(pts[1,i], pts[2,i], pts[3,i]) * r * LatticeDynamicsToolkit.TWOPI
#             kn = dot(v, ifc2.polar.eps * v)
#             s += exp(-kn * inv4) * kn
#         end
#         s - tol
#     end
#     y0 = fk(k0)
#     y1 = fk(k1)

#     write_metric(io, "ewald", "all", "lambda_from_file_ok", λ > 0 ? 1.0 : 0.0)
#     write_metric(io, "ewald", "all", "lambda", λ)
#     write_metric(io, "ewald", "all", "tol", tol)
#     write_metric(io, "ewald", "all", "k_inscribed", kins)
#     write_metric(io, "ewald", "all", "k0", k0)
#     write_metric(io, "ewald", "all", "k1", k1)
#     write_metric(io, "ewald", "all", "k_y0", y0)
#     write_metric(io, "ewald", "all", "k_y1", y1)
#     write_metric(io, "ewald", "all", "k_y0_times_y1", y0 * y1)
#     write_metric(io, "ewald", "all", "k_bracket_ok", y0 * y1 < 0 ? 1.0 : 0.0)

#     write_metric(io, "ewald", "all", "n_Rvector", size(ew.Rvec, 2))
#     write_metric(io, "ewald", "all", "n_Gvector", size(ew.Gvec, 2))
# end

# function run_polar_diagnostics_mgo()
#     base = joinpath(@__DIR__, "data", "MgO")
#     ucposcar = joinpath(base, "infile.ucposcar")
#     ssposcar = joinpath(base, "infile.ssposcar")
#     ifc2_path = joinpath(base, "infile.forceconstant")

#     uc = CrystalStructure(ucposcar)
#     sc = CrystalStructure(ssposcar)
#     ifc2_uc = read_ifc2(ifc2_path, ucposcar)

#     open(joinpath(@__DIR__, "outfile.polar_diagnostics_julia.csv"), "w") do io
#         println(io, "category,qpoint,metric,value")
#         write_metric(io, "meta", "all", "na_uc", length(uc))
#         write_metric(io, "meta", "all", "na_sc", length(sc))
#         write_metric(io, "meta", "all", "polar_enabled", ifc2_uc.has_polar_data ? 1.0 : 0.0)
#         write_metric(io, "meta", "all", "correctiontype", ifc2_uc.polar.correction_type)
#         write_metric(io, "eps", "all", "trace", tr(Matrix(ifc2_uc.polar.eps)))
#         write_metric(io, "eps", "all", "frob", norm(Matrix(ifc2_uc.polar.eps)))
#         write_metric(io, "born", "all", "born_frob", norm(ifc2_uc.polar.born_Z))

#         # Remap diagnostics (this is often where parity issues hide).
#         remap_ok = false
#         try
#             ifc2_sc = remap(sc, uc, ifc2_uc)[1]
#             fcp = DensePolarIFCs(ifc2_sc, uc, sc)
#             write_metric(io, "fcp_sc", "all", "frob", norm(fcp.fcp))
#             write_metric(io, "fcp_sc", "all", "herm_residual", tensor4_herm_residual(fcp.fcp))
#             remap_ok = true
#         catch
#             write_metric(io, "fcp_sc", "all", "remap_ok", 0.0)
#         end
#         remap_ok && write_metric(io, "fcp_sc", "all", "remap_ok", 1.0)

#         ew = LatticeDynamicsToolkit.precompute_ewald_parameters(ifc2_uc, uc)
#         born_onsite = born_onsite_correction(ew, ifc2_uc, uc)
#         write_ewald_probe(io, ew, uc, ifc2_uc)

#         for (qlabel, q) in (("GAMMA", Q_GAMMA), ("X", Q_X))
#             D_sr = Matrix(dynmat_q(ifc2_uc, uc, q; include_polar=false, ewald=ew))
#             D_skipna = copy(D_sr)
#             # For direct CSV parity with TDEP diagnostics, report q_cart without 2π.
#             q_cart = uc.L_inv' * q

#             # Build SR+Pt explicitly, independent of current dynmat_q gamma-path logic.
#             D_lr = LatticeDynamicsToolkit.longrange_dynamical_matrix(
#                 ew, uc, q, ifc2_uc.polar.born_Z, ifc2_uc.polar.eps;
#                 reconly=true,
#                 born_onsite=born_onsite,
#             )
#             add_massweighted_lr!(D_skipna, D_lr, uc)

#             D_with = copy(D_skipna)
#             if is_gamma_q(q)
#                 D_nac = zeros(ComplexF64, size(D_sr))
#                 LatticeDynamicsToolkit.add_nonanalytic_gamma!(D_nac, ifc2_uc, uc, LatticeDynamicsToolkit.q_cart_from_frac(uc, QDIR_GAMMA))
#                 D_with .+= D_nac
#             end

#             D_impl = Matrix(dynmat_q(ifc2_uc, uc, q; include_polar=true, qdir_gamma=QDIR_GAMMA, ewald=ew))

#             write_metric(io, "qpoint", qlabel, "q_frac_x", q[1])
#             write_metric(io, "qpoint", qlabel, "q_frac_y", q[2])
#             write_metric(io, "qpoint", qlabel, "q_frac_z", q[3])
#             write_metric(io, "qpoint", qlabel, "q_cart_x", q_cart[1])
#             write_metric(io, "qpoint", qlabel, "q_cart_y", q_cart[2])
#             write_metric(io, "qpoint", qlabel, "q_cart_z", q_cart[3])

#             write_matrix_norms_complex(io, "dm_with", qlabel, D_with)
#             write_matrix_norms_complex(io, "dm_skipna", qlabel, D_skipna)
#             write_matrix_norms_complex(io, "dm_sr", qlabel, D_sr)
#             write_metric(io, "dm_delta", qlabel, "frob_with_minus_skipna", norm(D_with .- D_skipna))
#             write_metric(io, "dm_delta", qlabel, "frob_skipna_minus_sr", norm(D_skipna .- D_sr))
#             write_metric(io, "dm_delta", qlabel, "frob_impl_minus_with", norm(D_impl .- D_with))

#             γ = Val(is_gamma_q(q))
#             ω_with = freqs_from_eigvals(get_modes(Hermitian(D_with), γ)[1])
#             ω_skipna = freqs_from_eigvals(get_modes(Hermitian(D_skipna), γ)[1])
#             ω_sr = freqs_from_eigvals(get_modes(Hermitian(D_sr), γ)[1])
#             write_freqs(io, qlabel, "freq_with", ω_with)
#             write_freqs(io, qlabel, "freq_skipna", ω_skipna)
#             write_freqs(io, qlabel, "freq_sr", ω_sr)
#         end
#     end

#     println("Wrote outfile.polar_diagnostics_julia.csv")
# end

# if abspath(PROGRAM_FILE) == @__FILE__
#     run_polar_diagnostics_mgo()
# end


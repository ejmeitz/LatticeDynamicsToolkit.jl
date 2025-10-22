export energies, make_energy_dataset

function energies(
        u::Vector{SVector{3,Float64}},
        fc2::IFC2;
        fc3::Union{Nothing,IFC3}=nothing,
        fc4::Union{Nothing,IFC4}=nothing,
        n_threads::Integer = Threads.nthreads()
    )

    na = length(u)

    maybe_check_na(i::Union{<:IFCs, Nothing}, na_true) = isnothing(i) ? true : na_true == i.na 

    if !all(maybe_check_na.((fc2, fc3, fc4), na))
        throw(ArgumentError("Displacements have $(na) atoms, but passed IFCs are built from different size cell. IFCs should be remapped to the supercell."))
    end

    e2, e3, e4 = @tasks for a1 in 1:na

        @set begin
            ntasks  = n_threads
            reducer = .+
        end

        @local v0 = MVector{3,Float64}(0,0,0)

        e2_local = 0.0
        e3_local = 0.0
        e4_local = 0.0
        
        v0 .= 0.0
        for pair in get_interactions(fc2, a1)
            a2  = pair.idxs[2]
            mul!(v0, pair.ifcs, u[a2], 1.0, 1.0)
        end
        e2_local += dot(u[a1], v0)

        # # Cubic term
        if fc3 !== nothing
            v0 .= 0.0
            for trip in get_interactions(fc3, a1)
                a2 = trip.idxs[2]; u2 = u[a2]
                a3 = trip.idxs[3]; u3 = u[a3]

                # einsum notation: V[i] = M[i, j, k] * u2[j] * u3[k]
                for i1 in 1:3, i2 in 1:3, i3 in 1:3
                    v0[i1] += trip.ifcs[i1,i2,i3] * u2[i2] * u3[i3]
                end

            end
            e3_local += dot(v0, u[a1])
        end

        # # quartic term
        if fc4 !== nothing
            v0 .= 0.0
            for quat in get_interactions(fc4, a1) # each quat is SVector{N, FC4Data{T}}
                a2 = quat.idxs[2]; u2 = u[a2]
                a3 = quat.idxs[3]; u3 = u[a3]
                a4 = quat.idxs[4]; u4 = u[a4]

                # einsum notation: V[i] = M[i, j, k, l] * u2[j] * u3[k] * u4[l]
                for i1 in 1:3, i2 in 1:3, i3 in 1:3, i4 in 1:3
                    v0[i1] += quat.ifcs[i1,i2,i3,i4] * u2[i2] * u3[i3] * u4[i4]
                end

            end
            e4_local += dot(v0, u[a1])
        end

        (e2_local, e3_local, e4_local)
    end

    e2 *= 0.5
    e3 /= 6.0
    e4 /= 24.0

    return e2, e3, e4
end

# Generate canonical_configrations and caluclate their energies from the Taylor series.
# Assumes IFCs are from unitcell if sc is passed
function make_energy_dataset(
        cc_settings::ConfigSettings,
        uc::CrystalStructure,
        sc::CrystalStructure;
        ifc2::IFC2, # required, but pass as kwarg
        ifc3::Union{Nothing, IFC3} = nothing,
        ifc4::Union{Nothing, IFC4} = nothing
    )

    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))
    
    @info "Remapping IFCs to Supercell"
    valid_ifcs_remapped = remap(sc, uc, valid_ifcs...)
    valid_ifcs_remapped_kwargs = build_kwargs(valid_ifcs_remapped...)
    
    return _make_energy_dataset(cc_settings, sc; valid_ifcs_remapped_kwargs...)
end

# Assumes IFCs are supercell already
function _make_energy_dataset(
    cc_settings::ConfigSettings,
    sc::CrystalStructure;
    ifc2::IFC2,
    ifc3::Union{Nothing, IFC3} = nothing,
    ifc4::Union{Nothing, IFC4} = nothing,
    n_threads::Integer = Threads.nthreads()
)
    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))

    remap_checks(sc, valid_ifcs...)

    dynmat = dynmat_gamma(ifc2, sc)
    freqs_sq, phi = get_modes(dynmat)
    freqs = sqrt.(freqs_sq)  # Will error for negative frequencies which I am ok with

    tep_energies = zeros(SVector{3, Float64}, cc_settings.n_configs)

    f = (config) -> energies(config, ifc2; fc3=ifc3, fc4=ifc4, n_threads=1)

    @info "Building Energy Dataset"
    canonical_configs(
        tep_energies,
        f,
        cc_settings,
        freqs,
        phi,
        sc.m;
        n_threads = n_threads
    )

    return tep_energies

end


# """
#     energies_and_forces_taylor(u, fc2; fc3=nothing, fc4=nothing)

# Compute Taylor-series energies (quadratic/cubic/quartic).

# Inputs
# - `u::AbstractMatrix{T}`: 3×na displacements (columns are atoms)
# - `fc2::IFCs{2,T}` (required)
# - `fc3::Union{Nothing,IFCs{3,T}}` (optional)
# - `fc4::Union{Nothing,IFCs{4,T}}` (optional)

# Returns
# - `(e2, e3, e4)`

# Conventions match your Fortran:
# - pair:  -½ * Σ_a u_a ⋅ ( Σ_pairs -m2 * u_b )
# - cubic: -⅓ * u_a ⋅ ( ½ * Σ_triplets Σ m3 * (u_b ⊗ u_c) )  and adds stress term
# - quart: -¼ * u_a ⋅ ( (1/6) * Σ_quartets Σ m4 * (u_b ⊗ u_c ⊗ u_d) )
# """
# function energies_faithful(
#         u::AbstractVector{SVector{3,Float64}},
#         fc2::IFC2;
#         fc3::Union{Nothing,IFC3}=nothing,
#         fc4::Union{Nothing,IFC4}=nothing
#     ) 

#     na = length(u)

#     maybe_check_na(i::Union{<:IFCs, Nothing}, na_true) = isnothing(i) ? true : na_true == i.na 

#     if !all(maybe_check_na.((fc2, fc3, fc4), na))
#         throw(ArgumentError("Displacements have $(na) atoms, but passed IFCs are built from different size cell."))
#     end

#     e2 = 0.0
#     e3 = 0.0
#     e4 = 0.0

#     v0 = MVector{3,Float64}(0,0,0)

#     # Quadratic term
#     @inbounds for a1 in 1:na
#         v0 .= 0.0
#         for pair in get_interactions(fc2, a1) # each pair is SVector{N, FC2Data{T}}
#             a2  = pair.idxs[2]
#             v0 -= pair.ifcs * u[a2]
#         end
#         e2 -= dot(u[a1], v0) * 0.5
#     end

#     # Cubic term
#     if fc3 !== nothing
#         @inbounds for a1 in 1:na
#             v0 .= 0.0
#             u1 = u[a1]
#             for trip in get_interactions(fc3, a1) # each pair is SVector{N, FC3Data{T}}
#                 a2 = trip.idxs[2]; u2 = u[a2]
#                 a3 = trip.idxs[3]; u3 = u[a3]

#                 @inbounds for i1 in 1:3, i2 in 1:3, i3 in 1:3
#                     v0[i1] -= trip.ifcs[i1,i2,i3] * u2[i2] * u3[i3]
#                 end
#             end
#             v0 .*= 0.5
#             e3 -= dot(v0, u1) / 3.0
#         end
#     end

#     # quartic term
#     if fc4 !== nothing
#         @inbounds for a1 in 1:na
#             v0 .= 0.0
#             u1 = u[a1]
#             for quat in get_interactions(fc4, a1) # each quat is SVector{N, FC4Data{T}}
#                 a2 = quat.idxs[2]; u2 = u[a2]
#                 a3 = quat.idxs[3]; u3 = u[a3]
#                 a4 = quat.idxs[4]; u4 = u[a4]

#                 @inbounds for i1 in 1:3, i2 in 1:3, i3 in 1:3, i4 in 1:3
#                     v0[i1] -= quat.ifcs[i1,i2,i3,i4] * u2[i2] * u3[i3] * u4[i4]
#                 end
#             end
#             v0 ./= 6.0
#             e4 -= dot(v0, u1) / 4.0
#         end
#     end

#     return e2, e3, e4
# end

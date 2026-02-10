export energies, make_energy_dataset

@inline function contract33!(out::MVector{3, T}, ifcs::SArray{Tuple{3,3,3},T,3,27},
                           u2::SVector{3,T}, u3::SVector{3,T}) where {T}
    @inbounds for i1 in 1:3
        # accumulate per row j = i2
        r1 = muladd(ifcs[i1,1,1],u3[1], muladd(ifcs[i1,1,2],u3[2], ifcs[i1,1,3]*u3[3]))
        r2 = muladd(ifcs[i1,2,1],u3[1], muladd(ifcs[i1,2,2],u3[2], ifcs[i1,2,3]*u3[3]))
        r3 = muladd(ifcs[i1,3,1],u3[1], muladd(ifcs[i1,3,2],u3[2], ifcs[i1,3,3]*u3[3]))
        out[i1] = muladd(u2[1], r1, muladd(u2[2], r2, muladd(u2[3], r3, out[i1])))
    end
    return out
end

#    r = A[:,:,:,i4]·u4  →  3×3;   t = r·u3  → 3;   v[i1]+= dot(u2,t)
@inline function contract44!(out::MVector{3,T}, quat::SArray{Tuple{3,3,3,3},T,4,81},
                                u2::SVector{3,T}, u3::SVector{3,T}, u4::SVector{3,T}) where {T}
    @inbounds for i1 in 1:3
        # acc over i2 after contracting i4 then i3
        acc = MVector{3,T}(0,0,0)
        for i2 in 1:3
            t = zero(T)
            for i3 in 1:3
                # r = dot(quat[i1,i2,i3,:], u4) with FMAs
                r = muladd(quat[i1,i2,i3,3], u4[3],
                    muladd(quat[i1,i2,i3,2], u4[2], quat[i1,i2,i3,1]*u4[1]))
                t = muladd(r, u3[i3], t)
            end
            acc[i2] = t
        end
        # v[i1] += dot(u2, acc)
        out[i1] = muladd(u2[3], acc[3], muladd(u2[2], acc[2], muladd(u2[1], acc[1], out[i1])))
    end
    return out
end


function energies(
        u::AbstractVector{SVector{3,Float64}},
        fc2::IFC2;
        fc3::Union{Nothing,IFC3}=nothing,
        fc4::Union{Nothing,IFC4}=nothing,
        fcp::Union{Nothing,DensePolarIFCs}=nothing,
        n_threads::Integer = Threads.nthreads()
    )

    na = length(u)

    maybe_check_na(i::Union{<:IFCs, Nothing}, na_true) = isnothing(i) ? true : na_true == i.na 

    if !all(maybe_check_na.((fc2, fc3, fc4), na))
        throw(ArgumentError("Displacements have $(na) atoms, but passed IFCs are built from different size cell. IFCs should be remapped to the supercell."))
    end

    has_fcp = fcp !== nothing
    fcp_mat = has_fcp ? fcp.fcp : nothing

    e2, e3, e4, ep = @tasks for a1 in 1:na

        @set begin
            ntasks  = n_threads
            reducer = .+
            scheduler = :static
        end

        @local v0 = MVector{3,Float64}(0,0,0)

        e2_local = 0.0
        e3_local = 0.0
        e4_local = 0.0
        ep_local = 0.0
        
        v0 .= 0.0
        for pair in get_interactions(fc2, a1)
            a2  = pair.idxs[2]
            mul!(v0, pair.ifcs, u[a2], 1.0, 1.0)
        end
        e2_local += dot(u[a1], v0)

        # Polar long-range term
        if has_fcp
            v0 .= 0.0
            for a2 in 1:na
                @inbounds for j in 1:3
                    v0[1] += fcp_mat[1, j, a1, a2] * u[a2][j]
                    v0[2] += fcp_mat[2, j, a1, a2] * u[a2][j]
                    v0[3] += fcp_mat[3, j, a1, a2] * u[a2][j]
                end
            end
            ep_local += dot(u[a1], v0)
        end

        # Cubic term
        if fc3 !== nothing
            v0 .= 0.0
            for trip in get_interactions(fc3, a1)
                a2 = trip.idxs[2]; u2 = u[a2]
                a3 = trip.idxs[3]; u3 = u[a3]

                contract33!(v0, trip.ifcs, u2, u3)

            end
            e3_local += dot(v0, u[a1])
        end

        # Quartic term
        if fc4 !== nothing
            v0 .= 0.0
            for quat in get_interactions(fc4, a1) # each quat is SVector{N, FC4Data{T}}
                a2 = quat.idxs[2]; u2 = u[a2]
                a3 = quat.idxs[3]; u3 = u[a3]
                a4 = quat.idxs[4]; u4 = u[a4]

                contract44!(v0, quat.ifcs, u2, u3, u4)
            end
            e4_local += dot(v0, u[a1])
        end

        (e2_local, e3_local, e4_local, ep_local)
    end

    e2 *= 0.5
    e3 /= 6.0
    e4 /= 24.0
    ep *= 0.5

    return e2, e3, e4, ep
end

function energies(
    u::AbstractVector{SVector{3,Float64}},
    fc2::AmorphousIFC2;
    fc3::Nothing = nothing,
    fc4::Nothing = nothing,
    fcp::Nothing = nothing,
    n_threads::Integer = Threads.nthreads()
)
    na = length(u)
    
    if fc2.na != na
        throw(ArgumentError("Displacements have $(na) atoms, but IFC2 has $(fc2.na) atoms."))
    end
    
    # Parallel reduction over atoms - each atom i owns pairs (i, j) where j > i
    # Compute E = 0.5 * u' * Φ * u directly without allocations
    # For each stored pair (i < j):
    #   - Off-diagonal: 2 * u_i' * Φ_ij * u_j (factor 2 for lower triangle symmetry)
    #   - Diagonal (ASR: Φ_ii = -Σ_j Φ_ij): -u_i' * Φ_ij * u_i - u_j' * Φ_ij' * u_j
    # All SMatrix/SVector ops are stack-allocated (no heap allocations)
    
    e2 = @tasks for i in 1:na
        @set begin
            ntasks = n_threads
            reducer = +
            scheduler = :static
        end
        
        e2_local = 0.0
        u_i = u[i]
        
        @inbounds for (k, j) in enumerate(fc2.neighbors[i])
            Φ_ij = fc2.blocks[i][k]
            u_j = u[j]
            
            # Off-diagonal contribution (×2 for symmetry with lower triangle)
            e2_local += 2.0 * dot(u_i, Φ_ij * u_j)
            
            # Diagonal contributions from ASR (Φ_ii = -Σ_{j≠i} Φ_ij)
            e2_local -= dot(u_i, Φ_ij * u_i)      # contribution to Φ_ii
            e2_local -= dot(u_j, Φ_ij' * u_j)     # contribution to Φ_jj
        end
        
        e2_local
    end
    
    e2 *= 0.5
    
    return e2, 0.0, 0.0, 0.0
end


# Generate canonical_configrations and caluclate their energies from the Taylor series.
# Assumes IFCs are from unitcell if sc is passed
function make_energy_dataset(
        cc_settings::ConfigSettings,
        uc::CrystalStructure,
        sc::CrystalStructure;
        ifc2::Union{IFC2, AmorphousIFC2}, # required, but pass as kwarg
        ifc3::Union{Nothing, IFC3} = nothing,
        ifc4::Union{Nothing, IFC4} = nothing,
        n_threads::Integer = Threads.nthreads(),
        verbose::Bool = true
    )

    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))

    if isa(ifc2, AmorphousIFC2) && (ifc3 !== nothing || ifc4 !== nothing)
        error(ArgumentError("Does not make sense to use AmorphousIFC2 with higher order IFCs to build energy dataset"))
    end
    
    verbose && @info "Remapping IFCs to Supercell"
    valid_ifcs_remapped = remap(sc, uc, valid_ifcs...)
    valid_ifcs_remapped_kwargs = build_kwargs(valid_ifcs_remapped...)

    # Build dense long-range polar data once (TDEP fcp equivalent)
    fcp = nothing
    if isa(ifc2, IFC2) && valid_ifcs_remapped_kwargs.ifc2.has_polar_data
        fcp = DensePolarIFCs(valid_ifcs_remapped_kwargs.ifc2, uc, sc)
    end

    return _make_energy_dataset(
            cc_settings,
            sc;
            valid_ifcs_remapped_kwargs...,
            fcp = fcp,
            n_threads = n_threads,
            verbose = verbose
        )
end

# Assumes IFCs are supercell already
function _make_energy_dataset(
    cc_settings::ConfigSettings,
    sc::CrystalStructure;
    ifc2::Union{IFC2, AmorphousIFC2},
    ifc3::Union{Nothing, IFC3} = nothing,
    ifc4::Union{Nothing, IFC4} = nothing,
    fcp::Union{Nothing,DensePolarIFCs} = nothing,
    n_threads::Integer = Threads.nthreads(),
    verbose::Bool = true
)
    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))

    remap_checks(sc, valid_ifcs...)

    dynmat = dynmat_gamma(ifc2, sc)
    freqs_sq, phi = get_modes(dynmat, Val{true}())
    freqs = sqrt.(freqs_sq)  # Will error for negative frequencies which I am ok with

    tep_energies = zeros(SVector{4, Float64}, cc_settings.n_configs)

    # Harmonic part for TEP is e2 + ep (short-range + polar); ep is 0 when fcp is nothing
    f = (config) -> begin
        e2, e3, e4, ep = energies(config, ifc2; fc3=ifc3, fc4=ifc4, fcp=fcp, n_threads=1)
        return SVector{4,Float64}(e2, e3, e4, ep)
    end

    verbose && @info "Building Energy Dataset"
    canonical_configs!(
        tep_energies,
        f,
        cc_settings,
        freqs,
        phi,
        sc.m;
        n_threads = n_threads,
        verbose = verbose
    )

    return Hartree_to_eV .* tep_energies

end


export energies

"""
    energies_and_forces_taylor(u, fc2; fc3=nothing, fc4=nothing)

Compute Taylor-series energies (quadratic/cubic/quartic).

Inputs
- `u::AbstractMatrix{T}`: 3×na displacements (columns are atoms)
- `fc2::IFCs{2,T}` (required)
- `fc3::Union{Nothing,IFCs{3,T}}` (optional)
- `fc4::Union{Nothing,IFCs{4,T}}` (optional)

Returns
- `(e2, e3, e4)`

Conventions match your Fortran:
- pair:  -½ * Σ_a u_a ⋅ ( Σ_pairs -m2 * u_b )
- cubic: -⅓ * u_a ⋅ ( ½ * Σ_triplets Σ m3 * (u_b ⊗ u_c) )  and adds stress term
- quart: -¼ * u_a ⋅ ( (1/6) * Σ_quartets Σ m4 * (u_b ⊗ u_c ⊗ u_d) )
"""
function energies_faithful(
        u::AbstractVector{SVector{3,T}},
        fc2::IFCs{2,T};
        fc3::Union{Nothing,IFCs{3,T}}=nothing,
        fc4::Union{Nothing,IFCs{4,T}}=nothing
    ) where {T <: AbstractFloat}

    na = length(u)

    maybe_check_na(i::Union{<:IFCs, Nothing}, na_true) = isnothing(i) ? true : na_true == i.na 

    if !all(maybe_check_na.((fc2, fc3, fc4), na))
        throw(ArgumentError("Displacements have $(na) atoms, but passed IFCs are built from different size cell."))
    end

    e2 = zero(T)
    e3 = zero(T)
    e4 = zero(T)

    v0 = MVector{3,T}(0,0,0)

    # Quadratic term
    @inbounds for a1 in 1:na
        v0 .= T(0.0)
        for pair in get_interactions(fc2, a1) # each pair is SVector{N, FC2Data{T}}
            a2  = pair.idxs[2]
            v0 -= pair.ifcs * u[a2]
        end
        e2 -= dot(u[a1], v0) * T(0.5)
    end

    # Cubic term
    if fc3 !== nothing
        @inbounds for a1 in 1:na
            v0 .= T(0.0)
            u1 = u[a1]
            for trip in get_interactions(fc3, a1) # each pair is SVector{N, FC3Data{T}}
                a2 = trip.idxs[2]; u2 = u[a2]
                a3 = trip.idxs[3]; u3 = u[a3]

                @inbounds for i1 in 1:3, i2 in 1:3, i3 in 1:3
                    v0[i1] -= trip.ifcs[i1,i2,i3] * u2[i2] * u3[i3]
                end
            end
            v0 .*= T(0.5)
            e3 -= dot(v0, u1) / T(3)
        end
    end

    # quartic term
    if fc4 !== nothing
        @inbounds for a1 in 1:na
            v0 .= T(0.0)
            u1 = u[a1]
            for quat in get_interactions(fc4, a1) # each quat is SVector{N, FC4Data{T}}
                a2 = quat.idxs[2]; u2 = u[a2]
                a3 = quat.idxs[3]; u3 = u[a3]
                a4 = quat.idxs[4]; u4 = u[a4]

                @inbounds for i1 in 1:3, i2 in 1:3, i3 in 1:3, i4 in 1:3
                    v0[i1] -= quat.ifcs[i1,i2,i3,i4] * u2[i2] * u3[i3] * u4[i4]
                end
            end
            v0 ./= T(6)
            e4 -= dot(v0, u1) / T(4)
        end
    end

    return e2, e3, e4
end

#! TRY LOOP VECTORIZATION??
#! Try storing IFCs as Tensorial types and using their einsum macros??
function energies_optimized(
        u::AbstractVector{SVector{3,T}},
        fc2::IFCs{2,T};
        fc3::Union{Nothing,IFCs{3,T}}=nothing,
        fc4::Union{Nothing,IFCs{4,T}}=nothing,
        n_threads::Integer = Threads.nthreads()
    ) where {T <: AbstractFloat}

    na = length(u)

    maybe_check_na(i::Union{<:IFCs, Nothing}, na_true) = ifelse(isnothing(i), true, na_true == i.na) 

    if !all(maybe_check_na.((fc2, fc3, fc4), na))
        throw(ArgumentError("Displacements have $(na) atoms, but passed IFCs are built from different size cell. IFCs should be remapped to the supercell."))
    end

    e2_local = zeros(T, n_threads)
    e3_local = zeros(T, n_threads)
    e4_local = zeros(T, n_threads)

    # Quadratic term
    e2, e3, e4 = @inbounds @tasks for a1 in 1:na

        @set begin
            ntasks  = n_threads
            reducer = .+
        end
        @local v0 = MVector{3,T}(0,0,0)

        e2_local = zero(T)
        e3_local = zero(T)
        e4_local = zero(T)

        v0 .= T(0.0)
        for pair in get_interactions(fc2, a1) # each pair is SVector{N, FC2Data{T}}
            a2  = pair.idxs[2]
            v0 += pair.ifcs * u[a2]
        end
        e2_local += dot(u[a1], v0)

        # Cubic term
        if fc3 !== nothing
            v0 .= T(0.0)
            u1 = u[a1]
            for trip in get_interactions(fc3, a1) # each pair is SVector{N, FC3Data{T}}
                a2 = trip.idxs[2]; u2 = u[a2]
                a3 = trip.idxs[3]; u3 = u[a3]

                #! use some einsum thing: V[i] = M[i, j, k] * u2[j] * u3[k]
                @inbounds for i1 in 1:3, i2 in 1:3, i3 in 1:3
                    v0[i1] += trip.ifcs[i1,i2,i3] * u2[i2] * u3[i3]
                end
            end
            e3_local += dot(v0, u1)
        end

        # quartic term
        if fc4 !== nothing
            v0 .= T(0.0)
            u1 = u[a1]
            for quat in get_interactions(fc4, a1) # each quat is SVector{N, FC4Data{T}}
                a2 = quat.idxs[2]; u2 = u[a2]
                a3 = quat.idxs[3]; u3 = u[a3]
                a4 = quat.idxs[4]; u4 = u[a4]

                #! use some einsum thing: V[i] = M[i, j, k, l] * u2[j] * u3[k] * u4[l]
                @inbounds for i1 in 1:3, i2 in 1:3, i3 in 1:3, i4 in 1:3
                    v0[i1] += quat.ifcs[i1,i2,i3,i4] * u2[i2] * u3[i3] * u4[i4]
                end
            end
            e4_local += dot(v0, u1)
        end

        (e2_local, e3_local, e4_local)
    end

    e2 *= T(0.5)
    e3 /= T(6)
    e4 /= T(24)

    return e2, e3, e4
end

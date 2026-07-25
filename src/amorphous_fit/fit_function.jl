using LinearAlgebra
using SparseArrays
using StaticArrays
using BSplineKit

abstract type AbstractIFC2FunctionModel end


"""
    IFC2PairEntry

Metadata for one independent off-diagonal IFC pair block.

`col_offset` is the first column in the original full IFC vector `x`
corresponding to this 3x3 block.

    col = col_offset + (β - 1) * 3 + α - 1

so `col_offset` should be 1-based.
"""
struct IFC2PairEntry
    col_offset::Int
    sp::Tuple{Symbol,Symbol}
    rij::SVector{3,Float64}
    T::Float64
end

@inline function ifc2_col(col_offset::Int, α::Int, β::Int)
    return col_offset + (β - 1) * 3 + α - 1
end

function build_ifc2_pair_entries(
    crys::CrystalStructure,
    neighbors::Vector{Vector{NeighborAtom}},
    pair_to_col,
    T::Float64,
)
    entries = IFC2PairEntry[]

    for i in eachindex(neighbors)
        for (k, nbr) in enumerate(neighbors[i])
            j = nbr.j
            col_offset = pair_to_col[i][k]

            species_pair = canonical_species_pair(
                crys.species[i],
                crys.species[j],
            )

            push!(
                entries,
                IFC2PairEntry(
                    col_offset,
                    species_pair,
                    nbr.r_ij,
                    T,
                ),
            )
        end
    end

    return entries
end


#####################


_tensor_coeff(::Type{LongitudinalChannel}, α::Int, β::Int, rr) = rr
_tensor_coeff(::Type{TransverseChannel}, α::Int, β::Int, rr) = ifelse(α == β, 1.0, 0.0) - rr
_tensor_coeff(ch, α::Int, β::Int, rr) = error("Unknown channel $ch")

@inline function tensor_coeff(ch::Type{<:ChannelIndex}, α::Int, β::Int, rhat::SVector{3,Float64})
    rr = rhat[α] * rhat[β]
    return _tensor_coeff(ch, α, β, rr)
end

function unit_vector(rij::SVector{3,Float64})
    r = norm(rij)
    r == 0.0 && throw(ArgumentError("rij has zero length"))
    return r, rij / r
end

function build_P(
    model::BSplineLTIFC2Model,
    entries::AbstractVector{<:IFC2PairEntry},
    n_ifc_unknowns::Int,
) 

    I = Int[]
    J = Int[]
    V = Float64[]

    # Each pair contributes up to:
    # 9 Cartesian components × 2 channels × radial_order × temperature_order.
    sizehint!(I, 9 * 2 * 16 * length(entries))
    sizehint!(J, 9 * 2 * 16 * length(entries))
    sizehint!(V, 9 * 2 * 16 * length(entries))

    channels = (LongitudinalChannel(), TransverseChannel())

    for e in entries
        #! DOUBLE CHECK r_ij is NEAREST IMAGE VECTOR
        r, rhat = unit_vector(e.rij)

        Br_vals = active_basis_values(model.radial_basis, r)
        BT_vals = active_basis_values(model.temperature_basis, e.T)

        for α in 1:3
            for β in 1:3
                row = ifc2_col(e.col_offset, α, β)

                for ch in channels
                    tc = tensor_coeff(ch, α, β, rhat)
                    tc == 0.0 && continue

                    for (ir, Br) in Br_vals
                        for (iT, BT) in BT_vals
                            col = theta_index(model, e.species_pair, ch, ir, iT)
                            val = tc * Br * BT

                            push!(I, row)
                            push!(J, col)
                            push!(V, val)
                        end
                    end
                end
            end
        end
    end

    return sparse(I, J, V, n_ifc_unknowns, ntheta(model))
end

struct FittedIFC2Function{M,V}
    model::M
    theta::V
end

function fit_ifc2_function(
    A::AbstractMatrix{Float64},
    b::AbstractVector{Float64},
    model::BSplineLTIFC2Model{S},
    entries::AbstractVector{IFC2PairEntry{S}};
    λ::Float64 = 0.0,
) where {S}

    P = build_P(model, entries, size(A, 2))
    Aeff = A * P

    nθ = size(P, 2)

    if λ == 0.0
        θ = Aeff \ b
    else
        R = sqrt(λ) * sparse(I, nθ, nθ)
        θ = [Aeff; R] \ [b; zeros(nθ)]
    end

    xfit = P * θ
    residual = Aeff * θ - b

    f = FittedIFC2Function(model, θ)

    return (; function_model = f, θ, xfit, P, Aeff, residual_norm = norm(residual))
end

function scalar_channel_value(
    model::BSplineLTIFC2Model{S},
    θ::AbstractVector{Float64},
    species_pair::Tuple{S,S},
    ch::Type{<:ChannelIndex},
    r::Float64,
    T::Float64,
) where {S}

    Br_vals = active_basis_values(model.radial_basis, r)
    BT_vals = active_basis_values(model.temperature_basis, T)

    out = 0.0

    for (ir, Br) in Br_vals
        for (iT, BT) in BT_vals
            idx = theta_index(model, species_pair, ch, ir, iT)
            out += θ[idx] * Br * BT
        end
    end

    return out
end

function (f::FittedIFC2Function{<:BSplineLTIFC2Model{S}})(
    rij::AbstractVector{<:Real},
    T::Real,
    species_pair::Tuple{S,S},
) where {S}

    rij_s = SVector{3,Float64}(rij)
    r, rhat = unit_vector(rij_s)

    model = f.model
    θ = f.theta

    kL = scalar_channel_value(model, θ, canonical_species_pair(species_pair), LongitudinalChannel, r, Float64(T))
    kT = scalar_channel_value(model, θ, canonical_species_pair(species_pair), TransverseChannel, r, Float64(T))

    rr = rhat * transpose(rhat)
    I3 = SMatrix{3,3,Float64}(I)

    return kL * rr + kT * (I3 - rr)
end



# API GOAL
# cutoffs = SpeciesCutoffs(
#     SpeciesPair(:Al, :Al, 9.0),
#     SpeciesPair(:Al, :Ni, 8.5),
#     SpeciesPair(:Ni, :Ni, 8.0),
# )

# model = BSplineLTIFC2Model(
#     cutoffs;
#     radial_nbreaks = 12,
#     temperature_breaks = range(300.0, 1200.0; length = 5),
#     radial_order = 4,
#     temperature_order = 3,
# )

# fit = fit_ifc2_function(
#     problem,
#     model;
#     λ = 1e-8,
# )

# Φ = fit.function_model

# rij = SVector{3,Float64}(1.2, 0.3, -0.1)
# Φij = Φ(rij, 800.0, (:Al, :Ni))
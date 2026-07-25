struct BSplineLTIFC2Model{BR, BT} <: AbstractIFC2FunctionModel
    species_pairs::Vector{SpeciesPair}
    species_pair_to_index::Dict{Tuple{Symbol, Symbol},Int}
    radial_basis::BR
    temperature_basis::BT
end

function BSplineLTIFC2Model(
    species_pairs::Vector{SpeciesPair},
    radial_breaks::AbstractVector{<:Real},
    temperature_breaks::AbstractVector{<:Real};
    radial_order::Int = 4,
    temperature_order::Int = 3,
) 
    spairs = unique(canonical_species_pair.(species_pairs))
    spair_to_index = Dict(sp => i for (i, sp) in enumerate(spairs))

    Br = BSplineBasis(BSplineOrder(radial_order), collect(Float64, radial_breaks))
    BT = BSplineBasis(BSplineOrder(temperature_order), collect(Float64, temperature_breaks))

    return BSplineLTIFC2Model(spairs, spair_to_index, Br, BT)
end

nr_basis(model::BSplineLTIFC2Model) = length(model.radial_basis)
nT_basis(model::BSplineLTIFC2Model) = length(model.temperature_basis)
n_species_pairs(model::BSplineLTIFC2Model) = length(model.species_pairs)
nchannels(::BSplineLTIFC2Model) = 2

function ntheta(model::BSplineLTIFC2Model)
    return n_species_pairs(model) *
           nchannels(model) *
           nr_basis(model) *
           nT_basis(model)
end


abstract type ChannelIndex end
struct LongitudinalChannel <: ChannelIndex end
struct TransverseChannel <: ChannelIndex end

channel_index(::Type{LongitudinalChannel}) = 1
channel_index(::Type{TransverseChannel}) = 2

function theta_index(
    model::BSplineLTIFC2Model,
    sp::SpeciesPair,
    ch::Type{<:ChannelIndex},
    ir::Int,
    iT::Int,
) where {S}

    isp = model.species_pair_to_index[canonical_species_pair(sp)]
    nr = nr_basis(model)
    nT = nT_basis(model)
    ich = channel_index(ch)

    return (((isp - 1) * nchannels(model) + (ich - 1)) * nr + (ir - 1)) * nT + iT
end

function active_basis_values(B, x::Real)
    xf = Float64(x)

    if !BSplineKit.isindomain(B, xf)
        xmin, xmax = BSplineKit.boundaries(B)
        throw(DomainError(xf, "value outside B-spline domain [$xmin, $xmax]"))
    end

    i, vals = B(xf)

    out = Tuple{Int,Float64}[]

    for q in eachindex(vals)
        ib = i - q + 1
        val = Float64(vals[q])

        if 1 <= ib <= length(B) && val != 0.0
            push!(out, (ib, val))
        end
    end

    return out
end
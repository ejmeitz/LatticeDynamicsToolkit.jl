
#! PARAMETRIC ON N OR NO?
struct SpeciesPair
    pair::Tuple{Symbol, Symbol}
    rc2::Float64
    rc3::Float64
    rc4::Float64
end

@inline canonical_species_pair(a::Symbol, b::Symbol) =
    string(a) <= string(b) ? (a, b) : (b, a)

function SpeciesPair(a::Symbol, b::Symbol, cutoffs::Vararg{<:Real,N}) where {N}
    0 < N <= 3 || throw(ArgumentError("SpeciesPair expects 1, 2, or 3 cutoffs. Got $N."))
    pair = canonical_species_pair(a,b)
    cutoffs = ntuple(i -> i <= N ? Float64(cutoffs[i]) : -1.0, 3)
    return SpeciesPair(pair, cutoffs...)
end

canonical_species_pair(sp::SpeciesPair) = sp.pair

# Does not actually check, just let it error
ifc2_cutoff(sp::SpeciesPair) = sp.rc2
ifc3_cutoff(sp::SpeciesPair) = sp.rc3
ifc4_cutoff(sp::SpeciesPair) = sp.rc4

struct SpeciesCutoffs
    pairs::Dict{Tuple{Symbol,Symbol},SpeciesPair}
end

function SpeciesCutoffs(pairs::SpeciesPair...)
    d = Dict{Tuple{Symbol,Symbol},SpeciesPair}()

    for sp in pairs
        key = canonical_species_pair(sp)
        haskey(d, key) && throw(ArgumentError("Duplicate SpeciesPair for $key"))
        d[key] = sp
    end

    return SpeciesCutoffs(d)
end

function get_species_pair(cutoffs::SpeciesCutoffs, a::Symbol, b::Symbol)
    key = canonical_species_pair(a, b)

    return get(cutoffs.pairs, key) do
        throw(KeyError(key))
    end
end

function ifc2_cutoff(cutoffs::SpeciesCutoffs, a::Symbol, b::Symbol)
    sp = get_species_pair(cutoffs, a, b)
    return ifc2_cutoff(sp)
end

maximum_ifc2_cutoff(cutoffs::SpeciesCutoffs) =
    maximum(ifc2_cutoff(sp) for sp in values(cutoffs.pairs))


struct NeighborAtom
    r_ij::SVector{3,Float64} # minimum image vector
    d2::Float64 # squared distance
    central_atom::Int # Central atom index
    j::Int # Neighbor index
end

function per_species_neighborlist(cutoffs::SpeciesCutoffs, crys::CrystalStructure)

    # Build neighborlist with maximum cutoff
    sys = CellListMap.ParticleSystem(
        xpositions = crys.x_cart,
        unitcell = crys.L,
        cutoff = maximum_ifc2_cutoff(cutoffs),
        output = 0.0,
        parallel = false #! current setup is not thread safe
    )

    neighbors = [NeighborPair[] for _ in 1:length(crys)]

    # Loops all unique (i,j) pairs (i.e., lower triangle)
    # Only keep pairs within their respective cutoffs
    CellListMap.pairwise!(sys) do pair, _
        species_i = crys.species[pair.i]
        species_j = crys.species[pair.j]
        sp = get_species_pair(cutoffs, species_i, species_j)
        rc_sq = ifc2_cutoff(sp)^2
        if pair.d2 <= rc_sq
            if pair.i < pair.j
                r_ij = SVector{3,Float64}(pair.y .- pair.x)
                nbr_atom = NeighborAtom(r_ij, pair.d2, pair.i, pair.j)
                push!(neighbors[pair.i], nbr_atom)
            else
                r_ji = SVector{3,Float64}(pair.x .- pair.y)
                nbr_atom = NeighborAtom(r_ji, pair.d2, pair.i, pair.j)
                push!(neighbors[pair.j], nbr_atom)
            end
        end
    end

    return neighbors
end
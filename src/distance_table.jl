function make_nl_bidirectional(nl::Vector{<:Tuple}, na::Int)
    nbrs = [Vector{Tuple{Int,Float64}}() for _ in 1:na]
    @inbounds for (i,j,d) in nl
        push!(nbrs[i], (j, d))
        push!(nbrs[j], (i, d)) 
    end
    return nbrs
end

function make_distance_table(
        crys::CrystalStructure{T},
        r_cut::T,
        tol::T = T(1e-4)
    ) where {T <: AbstractFloat}

    # Returns List of (idx1, idx2, dist) cannot assume ordering
    nl_oneway = CellListMap.neighborlist(crys.x_cart, r_cut + tol; unitcell = crys.L)
    nl = make_nl_bidirectional(nl_oneway, length(crys))

    dtas = Vector{DistanceTableAtom{T}}(undef, length(crys))
    Δf = MVector{3,T}(0,0,0)
    n = MVector{3,Int}(0,0,0)

    # Process into format similar to TDEP distance table
    for i in eachindex(crys.x_cart) # cannot parallelize naively

        nbrs = nl[i]
        n_neighbors = length(nbrs)

        vs = @SVector zeros(MVector{3,T}, n_neighbors)
        lvs = @SVector zeros(MVector{3,T}, n_neighbors)
        ns = @SVector zeros(MVector{3, Int16}, n_neighbors)

        # See lines 791 - 814 in lo_distancetable
        for (k,(j, _)) in enumerate(nbrs)
            Δf .= crys.x_frac[j] .- crys.x_frac[i]
            n   .= round.(Int16, Δf) # -1 or 0 or 1
            wrap = Δf .- n         # in [-0.5, 0.5)
            vs[k]  .= crys.L * wrap     # Cartesian minimum-image vector
            lvs[k] .= crys.L * (-n)   
            ns[k] .= n  
        end

        inds = SVector([nbrs[j][1] for j in 1:n_neighbors]...)
        dists = SVector([nbrs[j][2] for j in 1:n_neighbors]...)

        dtas[i] = DistanceTableAtom{T, n_neighbors}(
            i,
            SVector(vs),
            SVector(lvs),
            SVector(ns),
            inds,
            dists
        )

    end

    return DistanceTable{T}([dtas...])
end
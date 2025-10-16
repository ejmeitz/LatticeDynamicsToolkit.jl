function make_nl_bidirectional(nl::Vector{<:Tuple}, na::Int)
    nbrs = [Vector{Tuple{Int,Float64}}() for _ in 1:na]
    @inbounds for (i,j,d) in nl
        push!(nbrs[i], (j, d))
        push!(nbrs[j], (i, d)) 
    end
    return nbrs
end

function make_distance_table(
        crys::CrystalStructure,
        r_cut::Float64;
        include_self::Bool = false
    )

    # Returns List of (idx1, idx2, dist) cannot assume ordering
    nl_oneway = CellListMap.neighborlist(crys.x_cart, r_cut; unitcell = crys.L)
    nl = make_nl_bidirectional(nl_oneway, length(crys))

    dtas = Vector{DistanceTableAtom}(undef, length(crys))
    Δf = MVector{3,Float64}(0,0,0)
    n = MVector{3,Int}(0,0,0)

    # Process into format similar to TDEP distance table
    for i in eachindex(crys.x_cart) # cannot parallelize naively

        nbrs = nl[i]
        n_neighbors = length(nbrs) + (include_self ? 1 : 0)

        vs = @SVector zeros(MVector{3, Float64}, n_neighbors)
        lvs = @SVector zeros(MVector{3, Float64}, n_neighbors)
        ns = @SVector zeros(MVector{3, Int16}, n_neighbors)
        inds  = @MVector zeros(Int, n_neighbors)
        dists = @MVector zeros(Float64, n_neighbors)

        k = 1
        if include_self
            vs[k]    .= SVector{3, Float64}(0,0,0)
            lvs[k]   .= SVector{3, Float64}(0,0,0)
            ns[k]    .= SVector{3,Int16}(0,0,0)
            inds[k]  = i
            dists[k] = 0.0
            k += 1
        end

        # See lines 791 - 814 in lo_distancetable
        for (j, d) in nbrs
            Δf .= crys.x_frac[j] .- crys.x_frac[i]
            n   .= round.(Int16, Δf) # -1 or 0 or 1
            wrap = Δf .- n         # in [-0.5, 0.5)
            vs[k]  .= crys.L * wrap     # Cartesian minimum-image vector
            lvs[k] .= crys.L * (-n)   
            ns[k] .= n
            inds[k]  = j
            dists[k] = d  
            k += 1
        end


        dtas[i] = DistanceTableAtom{n_neighbors}(
            i,
            SVector(vs),
            SVector(lvs),
            SVector(ns),
            SVector(inds),
            SVector(dists)
        )

    end

    return DistanceTable([dtas...])
end
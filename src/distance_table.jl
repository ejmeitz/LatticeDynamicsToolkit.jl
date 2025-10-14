
function make_distance_table(
        x_cart::AbstractVector{SVector{3,T}},
        x_frac::AbstractVector{SVector{3,T}},
        L_ss::AbstractMatrix{T},
        r_cut::T,
        tol::T = T(1e-4)
    ) where {T <: AbstractFloat}

    # Returns List of (idx1, idx2, dist) cannot assume ordering
    nl = neighborlist(x_cart, r_cut + tol; unitcell = L_ss, showprogress = true)
    # Order by idx1
    sort!(nl; by = (x -> x[1]))

    left_nl_idx = 1
    right_nl_idx = 1
    N = length(n1)

    if length(x_cart) > max(nl[end]...)
        error("x_cart has $(length(x_cart)) atoms, but last entry in neighborlist has atom $(max(nl[end]...))")
    end

    dtas = Vector{DistanceTableAtom{T}}(undef, length(x_cart))
    Δf = MVector{3,T}(0,0,0)
    n = MVector{3,Int}(0,0,0)

    # Process into format similar to TDEP distance table
    for i in eachindex(x_cart) # cannot parallelize naively

        left_nl_idx = right_nl_idx

        #! can remove after testing
        if first(nl[right_nl_idx]) != i
            error("Ordering assumption wrong. On iter $i, but got atom $(nl[right_nl_idx])")
        end

        # Gather all neighbors of this atom
        while first(nl[right_nl_idx]) == i
            right_nl_idx += 1
            right_nl_idx == N && (right_nl_idx += 1; break)
        end

        nbrs = @views nl[left_nl_idx : right_nl_idx - 1]
        n_neighbors = length(nbrs)

        vs = @SVector zeros(SVector{3,T}, n_neighbors)
        lvs = @SVector zeros(SVector{3,T}, n_neighbors)
        ns = @SVector zeros(SVector{3, Int16}, n_neighbors)

        # See lines 791 - 814 in lo_distancetable
        for i in 1:n_neighbors
            idx1 = nbrs[i][1]; idx2 = nbrs[i][2]
            Δf .= x_frac[idx2] .- x_frac_i[idx1]
            n   .= round.(Int16, Δf) # -1 or 0 or 1
            wrap = Δf .- n         # in [-0.5, 0.5)
            vs[i]  += L * wrap     # Cartesian minimum-image vector
            lvs[i] += L * (-n)   
            ns[i] += n  
        end

        inds = SVector([nbrs[i][2] for i in 1:n_neighbors])
        dists = SVector([nbrs[i][3] for i in 1:n_neighbors])

        dtas[i] = DistaceTableAtom{T, n_neighbors}(
            i,
            vs,
            lvs,
            inds,
            dists
        )

    end

    return DistanceTable{T}([dtas...])
end
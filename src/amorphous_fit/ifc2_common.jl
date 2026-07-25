function _build_upper_neighbor_list(crystal::CrystalStructure, r_cut::Float64)
    na = length(crystal)
    nl_pairs = CellListMap.neighborlist(crystal.x_cart, r_cut; unitcell=crystal.L)
    neighbors = [Int[] for _ in 1:na]

    for (i, j, _) in nl_pairs
        if i < j
            push!(neighbors[i], j)
        else
            push!(neighbors[j], i)
        end
    end

    for i in 1:na
        sort!(neighbors[i])
    end
    return neighbors
end

function _build_pair_to_col(neighbors::Vector{Vector{Int}})
    na = length(neighbors)
    pair_to_col = Vector{Vector{Int}}(undef, na)
    col = 1
    for i in 1:na
        pair_to_col[i] = Vector{Int}(undef, length(neighbors[i]))
        for k in eachindex(neighbors[i])
            pair_to_col[i][k] = col
            col += 9
        end
    end
    return pair_to_col, col - 1
end

function _sample_position(problem::AmorphousFitProblem, sample::Int, atom::Int)
    col0 = 3 * (atom - 1)
    return SVector{3,Float64}(
        problem.positions[sample, col0 + 1],
        problem.positions[sample, col0 + 2],
        problem.positions[sample, col0 + 3],
    )
end

function _sample_force(problem::AmorphousFitProblem, sample::Int, atom::Int, component::Int)
    return problem.forces[sample, 3 * (atom - 1) + component]
end

function _displacement(crystal::CrystalStructure, x::SVector{3,Float64}, x0::SVector{3,Float64})
    Δx_frac = crystal.L_inv * (x - x0)
    return crystal.L * (Δx_frac - round.(Δx_frac))
end

function _pack_ifc2(
    φ::AbstractVector{Float64},
    r_cut::Float64,
    neighbors::Vector{Vector{Int}},
    pair_to_col::Vector{Vector{Int}},
)
    na = length(neighbors)
    blocks = Vector{Vector{SMatrix{3,3,Float64,9}}}(undef, na)
    @inbounds for i in 1:na
        blocks[i] = Vector{SMatrix{3,3,Float64,9}}(undef, length(neighbors[i]))
        for k in eachindex(neighbors[i])
            col0 = pair_to_col[i][k]
            blocks[i][k] = SMatrix{3,3,Float64,9}(φ[col0:col0+8])
        end
    end
    return AmorphousIFC2(na, r_cut, neighbors, blocks)
end

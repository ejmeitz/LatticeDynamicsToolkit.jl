export fit_ifc2

"""
    build_upper_neighbor_list(crystal::CrystalStructure, r_cut::Float64) 
        -> Vector{Vector{Int}}

Build neighbor list containing only upper-triangle pairs (i < j).
Uses CellListMap with periodic boundary conditions.

Returns `neighbors` where `neighbors[i]` contains indices `j > i` within `r_cut`.
"""
function build_upper_neighbor_list(crystal::CrystalStructure, r_cut::Float64)
    na = length(crystal)
    
    # CellListMap returns pairs (i, j, dist) - order is NOT guaranteed
    nl_pairs = CellListMap.neighborlist(crystal.x_cart, r_cut; unitcell = crystal.L)
    
    # Convert to adjacency list format, ensuring we store with smaller index first
    neighbors = [Int[] for _ in 1:na]
    
    for (i, j, _) in nl_pairs
        # Ensure we always store as (smaller, larger) for upper triangle
        if i < j
            push!(neighbors[i], j)
        else
            push!(neighbors[j], i)
        end
    end
    
    # Sort each neighbor list for consistent ordering
    for i in 1:na
        sort!(neighbors[i])
    end
    
    return neighbors
end

function build_pair_to_col(neighbors::Vector{Vector{Int}})
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

function compute_displacements!(
    u_t::Vector{SVector{3,Float64}},
    x_t::Vector{SVector{3,Float64}},
    crystal,
    positions::Matrix{Float64},
    x0::Vector{SVector{3,Float64}},
    t::Int,
    na::Int
)
    L    = crystal.L
    Linv = crystal.L_inv

    @inbounds for i in 1:na
        idx = (t - 1) * na + i

        xi = SVector{3,Float64}(positions[1, idx], positions[2, idx], positions[3, idx])
        x_t[i] = xi

        dx1 = xi[1] - x0[i][1]
        dx2 = xi[2] - x0[i][2]
        dx3 = xi[3] - x0[i][3]

        f1 = Linv[1,1]*dx1 + Linv[1,2]*dx2 + Linv[1,3]*dx3
        f2 = Linv[2,1]*dx1 + Linv[2,2]*dx2 + Linv[2,3]*dx3
        f3 = Linv[3,1]*dx1 + Linv[3,2]*dx2 + Linv[3,3]*dx3

        f1w = f1 - round(f1)
        f2w = f2 - round(f2)
        f3w = f3 - round(f3)

        u1 = L[1,1]*f1w + L[1,2]*f2w + L[1,3]*f3w
        u2 = L[2,1]*f1w + L[2,2]*f2w + L[2,3]*f3w
        u3 = L[3,1]*f1w + L[3,2]*f2w + L[3,3]*f3w

        u_t[i] = SVector{3,Float64}(u1, u2, u3)
    end

    return nothing
end

"""
Compute mean absolute antisymmetry of onsite ASR sum:
S_i = Σ_{j≠i} H_{ij} where H_{ij} is the block used in the assembled Hessian.
For your storage (i<j stores Φ_ij), this means:
- S_i gets +Φ_ij
- S_j gets +Φ_ij^T

Returns mean_i ||antisym(S_i)||_F and max_i ||antisym(S_i)||_F.
"""
function onsite_antisymmetry_stats(
    na::Int,
    neighbors::Vector{Vector{Int}},
    blocks::Vector{Vector{SMatrix{3,3,Float64,9}}}
)
    S = Vector{MMatrix{3,3,Float64,9}}(undef, na)
    @inbounds for i in 1:na
        S[i] = @MMatrix zeros(3,3)
    end

    @inbounds for i in 1:na
        for (k, j) in enumerate(neighbors[i])
            Φ = blocks[i][k]

            # S_i += Φ
            Si = S[i]
            Si[1,1] += Φ[1,1]; Si[1,2] += Φ[1,2]; Si[1,3] += Φ[1,3]
            Si[2,1] += Φ[2,1]; Si[2,2] += Φ[2,2]; Si[2,3] += Φ[2,3]
            Si[3,1] += Φ[3,1]; Si[3,2] += Φ[3,2]; Si[3,3] += Φ[3,3]
            S[i] = Si

            # S_j += Φ^T
            Sj = S[j]
            Sj[1,1] += Φ[1,1]; Sj[1,2] += Φ[2,1]; Sj[1,3] += Φ[3,1]
            Sj[2,1] += Φ[1,2]; Sj[2,2] += Φ[2,2]; Sj[2,3] += Φ[3,2]
            Sj[3,1] += Φ[1,3]; Sj[3,2] += Φ[2,3]; Sj[3,3] += Φ[3,3]
            S[j] = Sj
        end
    end

    sum_abs = 0.0
    max_abs = 0.0

    @inbounds for i in 1:na
        Si = S[i]
        a12 = 0.5 * (Si[1,2] - Si[2,1])
        a13 = 0.5 * (Si[1,3] - Si[3,1])
        a23 = 0.5 * (Si[2,3] - Si[3,2])
        absA = sqrt(2.0 * (a12*a12 + a13*a13 + a23*a23))
        sum_abs += absA
        max_abs = max(max_abs, absA)
    end

    return sum_abs / na, max_abs
end

# ----------------------------
# Augmented operator: [A; sqrtμ C]
# ----------------------------

struct IFC2AugOperator{TCrystal, TMat} <: AbstractMatrix{Float64}
    crystal::TCrystal
    positions::TMat
    neighbors::Vector{Vector{Int}}
    pair_to_col::Vector{Vector{Int}}
    x0::Vector{SVector{3,Float64}}

    na::Int
    n_samples::Int
    n_pairs::Int
    n_unknowns::Int
    n_equations::Int
    n_constraints::Int
    sqrtμ::Float64

    # scratch
    x_t::Vector{SVector{3,Float64}}
    u_t::Vector{SVector{3,Float64}}
end

Base.size(A::IFC2AugOperator) = (A.n_equations + A.n_constraints, A.n_unknowns)

# Indices in the 9-vector (column-major) for (α,β) and (β,α)
# Φ[1,2] index = 4, Φ[2,1] index = 2
const IDX_12 = 4; const IDX_21 = 2
const IDX_13 = 7; const IDX_31 = 3
const IDX_23 = 8; const IDX_32 = 6

# y := [A*x; sqrtμ*C*x]
function LinearAlgebra.mul!(y::AbstractVector{Float64}, A::IFC2AugOperator, x::AbstractVector{Float64})
    @assert length(x) == A.n_unknowns
    @assert length(y) == A.n_equations + A.n_constraints
    fill!(y, 0.0)

    na = A.na
    ns = A.n_samples
    neighbors = A.neighbors
    pair_to_col = A.pair_to_col

    # --- top part: A*x (force equations) ---
    @inbounds for t in 1:ns
        compute_displacements!(A.u_t, A.x_t, A.crystal, A.positions, A.x0, t, na)
        u_t = A.u_t

        for i in 1:na
            ui = u_t[i]
            row_i0 = (t - 1) * 3 * na + 3 * (i - 1)

            for (k, j) in enumerate(neighbors[i])
                uj = u_t[j]
                du1 = ui[1] - uj[1]
                du2 = ui[2] - uj[2]
                du3 = ui[3] - uj[3]

                col0 = pair_to_col[i][k]

                v1 = x[col0 + 0]
                v2 = x[col0 + 1]
                v3 = x[col0 + 2]
                v4 = x[col0 + 3]
                v5 = x[col0 + 4]
                v6 = x[col0 + 5]
                v7 = x[col0 + 6]
                v8 = x[col0 + 7]
                v9 = x[col0 + 8]

                row_j0 = (t - 1) * 3 * na + 3 * (j - 1)

                # F_i += Φ_ij * Δu
                y[row_i0 + 1] += v1*du1 + v4*du2 + v7*du3
                y[row_i0 + 2] += v2*du1 + v5*du2 + v8*du3
                y[row_i0 + 3] += v3*du1 + v6*du2 + v9*du3

                # F_j += -Φ_ij^T * Δu
                y[row_j0 + 1] += -(v1*du1 + v2*du2 + v3*du3)
                y[row_j0 + 2] += -(v4*du1 + v5*du2 + v6*du3)
                y[row_j0 + 3] += -(v7*du1 + v8*du2 + v9*du3)
            end
        end
    end

    # --- bottom part: sqrtμ * C*x (3 constraints per atom) ---
    off = A.n_equations
    sμ = A.sqrtμ

    @inbounds for i in 1:na
        for (k, j) in enumerate(neighbors[i])
            col0 = pair_to_col[i][k]

            # d_ab = Φ[α,β] - Φ[β,α]
            d12 = x[col0 + (IDX_12-1)] - x[col0 + (IDX_21-1)]
            d13 = x[col0 + (IDX_13-1)] - x[col0 + (IDX_31-1)]
            d23 = x[col0 + (IDX_23-1)] - x[col0 + (IDX_32-1)]

            ri = off + 3*(i-1)
            rj = off + 3*(j-1)

            # For atom i: +d ; for atom j: -d  (because S_j gets Φ^T)
            y[ri + 1] += sμ * d12
            y[ri + 2] += sμ * d13
            y[ri + 3] += sμ * d23

            y[rj + 1] += -sμ * d12
            y[rj + 2] += -sμ * d13
            y[rj + 3] += -sμ * d23
        end
    end

    return y
end

# g := Aaug' * y = A' * y_top + sqrtμ * C' * y_bottom
function LinearAlgebra.mul!(g::AbstractVector{Float64}, At::Adjoint{<:Any,<:IFC2AugOperator}, y::AbstractVector{Float64})
    A = At.parent
    @assert length(g) == A.n_unknowns
    @assert length(y) == A.n_equations + A.n_constraints
    fill!(g, 0.0)

    na = A.na
    ns = A.n_samples
    neighbors = A.neighbors
    pair_to_col = A.pair_to_col

    # --- top contribution: A' * y_top ---
    @inbounds for t in 1:ns
        compute_displacements!(A.u_t, A.x_t, A.crystal, A.positions, A.x0, t, na)
        u_t = A.u_t

        for i in 1:na
            ui = u_t[i]
            row_i0 = (t - 1) * 3 * na + 3 * (i - 1)

            yi1 = y[row_i0 + 1]
            yi2 = y[row_i0 + 2]
            yi3 = y[row_i0 + 3]

            for (k, j) in enumerate(neighbors[i])
                uj = u_t[j]
                du1 = ui[1] - uj[1]
                du2 = ui[2] - uj[2]
                du3 = ui[3] - uj[3]

                col0 = pair_to_col[i][k]

                row_j0 = (t - 1) * 3 * na + 3 * (j - 1)
                yj1 = y[row_j0 + 1]
                yj2 = y[row_j0 + 2]
                yj3 = y[row_j0 + 3]

                # From i-rows
                g[col0 + 0] += yi1 * du1
                g[col0 + 1] += yi2 * du1
                g[col0 + 2] += yi3 * du1

                g[col0 + 3] += yi1 * du2
                g[col0 + 4] += yi2 * du2
                g[col0 + 5] += yi3 * du2

                g[col0 + 6] += yi1 * du3
                g[col0 + 7] += yi2 * du3
                g[col0 + 8] += yi3 * du3

                # From j-rows (transpose mapping with minus sign)
                g[col0 + 0] += yj1 * (-du1)
                g[col0 + 1] += yj1 * (-du2)
                g[col0 + 2] += yj1 * (-du3)

                g[col0 + 3] += yj2 * (-du1)
                g[col0 + 4] += yj2 * (-du2)
                g[col0 + 5] += yj2 * (-du3)

                g[col0 + 6] += yj3 * (-du1)
                g[col0 + 7] += yj3 * (-du2)
                g[col0 + 8] += yj3 * (-du3)
            end
        end
    end

    # --- bottom contribution: sqrtμ * C' * y_bottom ---
    off = A.n_equations
    sμ = A.sqrtμ

    @inbounds for i in 1:na
        ri = off + 3*(i-1)
        yi12 = y[ri + 1]
        yi13 = y[ri + 2]
        yi23 = y[ri + 3]

        for (k, j) in enumerate(neighbors[i])
            col0 = pair_to_col[i][k]
            rj = off + 3*(j-1)

            # delta = y_i - y_j
            d12 = yi12 - y[rj + 1]
            d13 = yi13 - y[rj + 2]
            d23 = yi23 - y[rj + 3]

            # For each constraint component:
            # Φ[αβ] gets +delta, Φ[βα] gets -delta
            g[col0 + (IDX_12-1)] += sμ * d12
            g[col0 + (IDX_21-1)] += -sμ * d12

            g[col0 + (IDX_13-1)] += sμ * d13
            g[col0 + (IDX_31-1)] += -sμ * d13

            g[col0 + (IDX_23-1)] += sμ * d23
            g[col0 + (IDX_32-1)] += -sμ * d23
        end
    end

    return g
end

# 5-arg mul! overloads (avoid fallback paths)
function LinearAlgebra.mul!(
    y::AbstractVector{Float64},
    A::IFC2AugOperator,
    x::AbstractVector{Float64},
    α::Number,
    β::Number
)
    if β == 0
        mul!(y, A, x)
        if α != 1
            @inbounds @simd for i in eachindex(y)
                y[i] = α * y[i]
            end
        end
        return y
    else
        tmp = similar(y)
        mul!(tmp, A, x)
        @inbounds @simd for i in eachindex(y)
            y[i] = β * y[i] + α * tmp[i]
        end
        return y
    end
end

function LinearAlgebra.mul!(
    g::AbstractVector{Float64},
    At::Adjoint{<:Any,<:IFC2AugOperator},
    y::AbstractVector{Float64},
    α::Number,
    β::Number
)
    if β == 0
        mul!(g, At, y)
        if α != 1
            @inbounds @simd for i in eachindex(g)
                g[i] = α * g[i]
            end
        end
        return g
    else
        tmp = similar(g)
        mul!(tmp, At, y)
        @inbounds @simd for i in eachindex(g)
            g[i] = β * g[i] + α * tmp[i]
        end
        return g
    end
end

function make_ifc2_aug_operator(
    crystal::CrystalStructure,
    positions::Matrix{Float64};
    r_cut::Float64,
    μ::Float64
)
    na = length(crystal)
    n_samples = size(positions, 2) ÷ na
    @assert size(positions) == (3, na * n_samples)

    neighbors = build_upper_neighbor_list(crystal, r_cut)
    n_pairs = sum(length.(neighbors))

    pair_to_col, n_unknowns = build_pair_to_col(neighbors)
    n_equations = 3 * na * n_samples
    n_constraints = 3 * na

    x0 = Vector{SVector{3,Float64}}(undef, na)
    @inbounds for i in 1:na
        xi = crystal.x_cart[i]
        x0[i] = SVector{3,Float64}(xi[1], xi[2], xi[3])
    end

    x_t = Vector{SVector{3,Float64}}(undef, na)
    u_t = Vector{SVector{3,Float64}}(undef, na)

    return IFC2AugOperator(
        crystal, positions, neighbors, pair_to_col, x0,
        na, n_samples, n_pairs, n_unknowns, n_equations, n_constraints,
        sqrt(μ),
        x_t, u_t
    )
end

# ----------------------------
# Fit function (penalty constraints)
# ----------------------------

function fit_ifc2(
    crystal::CrystalStructure,
    positions::Matrix{Float64},
    forces::Matrix{Float64},
    r_cut::Float64;
    λ::Float64 = 1e-6,
    μ::Float64 = 1e2,  # larger enforces hermicity more strictly
    rtol::Float64 = 1e-4,
    atol::Float64 = 0.0,
    maxiter::Int = 2000,
    verbose::Bool = true,
    verbosity::Int = KrylovKit.STARTSTOP_LEVEL
)
    na = length(crystal)
    n_samples = size(positions, 2) ÷ na
    @assert size(positions) == (3, na * n_samples) "positions must be (3, na*n_samples)"
    @assert size(forces)    == (3, na * n_samples) "forces must be (3, na*n_samples)"

    Aop = make_ifc2_aug_operator(crystal, positions; r_cut=r_cut, μ=μ)

    F = vec(forces)
    b = zeros(Float64, Aop.n_equations + Aop.n_constraints)
    @inbounds b[1:Aop.n_equations] .= F

    verbose && @info "Fitting IFC2 with LSMR (matrix-free) + onsite-sum symmetry penalty" n_atoms=na n_samples=n_samples n_pairs=Aop.n_pairs n_unknowns=Aop.n_unknowns n_equations=Aop.n_equations n_constraints=Aop.n_constraints λ=λ μ=μ rtol=rtol

    φ, info = lssolve(Aop, b, λ; rtol=rtol, atol=atol, maxiter=maxiter, verbosity=verbosity)

    verbose && @info "LSMR done" converged=info.converged numiter=info.numiter numops=info.numops normres=info.normres

    # Fit quality on force equations only
    yfit = zeros(Float64, Aop.n_equations + Aop.n_constraints)
    mul!(yfit, Aop, φ)
    r_force = @views yfit[1:Aop.n_equations] .- F
    rmse = sqrt(mean(r_force .^ 2))
    max_err = maximum(abs.(r_force))
    verbose && @info "Fit quality (forces)" rmse=rmse max_error=max_err

    # Unpack φ into SMatrix blocks
    neighbors = Aop.neighbors
    pair_to_col = Aop.pair_to_col
    blocks = Vector{Vector{SMatrix{3,3,Float64,9}}}(undef, na)
    @inbounds for i in 1:na
        n_nbrs = length(neighbors[i])
        blocks[i] = Vector{SMatrix{3,3,Float64,9}}(undef, n_nbrs)
        for k in 1:n_nbrs
            col0 = pair_to_col[i][k]
            blocks[i][k] = SMatrix{3,3,Float64,9}(φ[col0:col0+8])
        end
    end

    # Diagnostic: antisymmetry of onsite ASR sum (this is exactly what μ penalizes)
    mean_abs, max_abs = onsite_antisymmetry_stats(na, neighbors, blocks)
    verbose && @info "Onsite-sum antisymmetry (ASR diagonal symmetry diagnostic)" mean_abs=mean_abs max_abs=max_abs

    return AmorphousIFC2(na, r_cut, neighbors, blocks), info
end

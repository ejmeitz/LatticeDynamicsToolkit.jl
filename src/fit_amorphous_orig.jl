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


"""
    fit_ifc2(crystal, positions, forces; r_cut, λ=1e-6, hard_constraints=true, μ=1e6) -> AmorphousIFC2

Fit second-order interatomic force constants for amorphous/disordered systems.

# Arguments
- `crystal::CrystalStructure`: Reference structure (equilibrium positions)
- `positions::Matrix{Float64}`: Cartesian positions `(3, na * n_timesteps)` in bohr
- `forces::Matrix{Float64}`: Cartesian forces `(3, na * n_timesteps)` in Hartree/bohr
- `r_cut::Float64`: Cutoff radius for IFC interactions (bohr)
- `λ::Float64`: Ridge regularization parameter (default 1e-6)
- `hard_constraints::Bool`: If true, use Lagrange multipliers for exact Hessian symmetry (default true).
                            If false, use penalty method with weight μ.
- `μ::Float64`: Hessian symmetry penalty weight, only used if hard_constraints=false (default 1e6).

# Returns
`AmorphousIFC2` with fitted force constants satisfying:
- Acoustic sum rule: Φ_ii = -Σ_{j≠i} Φ_ij (exact, by construction)
- Newton's 3rd law: Φ_ij = Φ_ji^T (exact, by construction)
- Hessian symmetry: Φ_ii = Φ_ii^T (exact with Lagrange multipliers, approximate with penalty)

# Method
Solves the linear system `F_i = Σ_{j≠i} Φ_ij · (u_i - u_j)` where:
- Only upper-triangle pairs (i < j) are independent unknowns
- Relative displacements `Δu = u_i - u_j` automatically enforce ASR
- Symmetric contribution to both F_i and F_j enforces Newton's 3rd law
- Hessian symmetry constraint: Σ_j Φ_ij^{αβ} = Σ_j Φ_ij^{βα} (symmetric diagonal blocks)
"""
function fit_ifc2(
    crystal::CrystalStructure,
    positions::Matrix{Float64},
    forces::Matrix{Float64};
    r_cut::Float64,
    λ::Float64 = 1e-6,
    hard_constraints::Bool = true,
    μ::Float64 = 1e6,
    verbose::Bool = true
)
    na = length(crystal)
    n_samples = size(positions, 2) ÷ na
    @assert size(positions) == (3, na * n_samples) "positions must be (3, na * n_timesteps)"
    @assert size(forces) == (3, na * n_samples) "forces must be (3, na * n_timesteps)"
    
    # Build neighbor list (upper triangle only)
    neighbors = build_upper_neighbor_list(crystal, r_cut)
    
    # Count total number of unique pairs
    n_pairs = sum(length.(neighbors))
    n_unknowns = 9 * n_pairs  # 9 elements per 3×3 block
    n_equations = 3 * na * n_samples
    
    verbose && @info "Fitting AmorphousIFC2" n_atoms=na n_pairs=n_pairs n_samples=n_samples n_unknowns=n_unknowns n_equations=n_equations
    
    # Build mapping from (i, k) -> column offset in design matrix
    # where k is the index into neighbors[i]
    pair_to_col = Vector{Vector{Int}}(undef, na)
    col = 1
    for i in 1:na
        pair_to_col[i] = Vector{Int}(undef, length(neighbors[i]))
        for k in eachindex(neighbors[i])
            pair_to_col[i][k] = col
            col += 9
        end
    end
    
    # Build design matrix in COO format
    I_coo = Int[]
    J_coo = Int[]
    V_coo = Float64[]
    
    # Pre-allocate estimate
    sizehint!(I_coo, 18 * n_pairs * n_samples)  # Each pair contributes to 2 atoms × 3 components × 9 unknowns
    sizehint!(J_coo, 18 * n_pairs * n_samples)
    sizehint!(V_coo, 18 * n_pairs * n_samples)
    
    # Equilibrium positions
    x0 = crystal.x_cart
    
    # Loop over timesteps and build design matrix
    for t in 1:n_samples
        # Extract positions for this timestep
        x_t = [SVector{3,Float64}(positions[:, (t-1)*na + i]) for i in 1:na]
        
        # Compute displacements u = x - x0 (with minimum image convention)
        u_t = Vector{SVector{3,Float64}}(undef, na)
        for i in 1:na
            # Displacement in fractional coords, then wrap to minimum image
            Δx_frac = crystal.L_inv * (x_t[i] - x0[i])
            Δx_frac_wrapped = Δx_frac - round.(Δx_frac)
            u_t[i] = crystal.L * Δx_frac_wrapped
        end
        
        # For each unique pair (i, j) with i < j
        for i in 1:na
            for (k, j) in enumerate(neighbors[i])
                # Relative displacement Δu = u_i - u_j
                Δu = u_t[i] - u_t[j]
                
                # Column offset for this pair's 9 unknowns
                col_offset = pair_to_col[i][k]
                
                # Row offsets for atoms i and j in this timestep
                row_i = (t-1) * 3 * na + 3*(i-1)
                row_j = (t-1) * 3 * na + 3*(j-1)
                
                # Contribution to F_i: (Δu^T ⊗ I₃) · φ_ij
                # For α ∈ {1,2,3}: F_i[α] = Σ_β Φ_ij[α,β] * Δu[β]
                # With column-major storage: φ_ij = [Φ[1,1], Φ[2,1], Φ[3,1], Φ[1,2], ...]
                for α in 1:3
                    row = row_i + α
                    for β in 1:3
                        col = col_offset + (β-1)*3 + α - 1  # column-major index
                        push!(I_coo, row)
                        push!(J_coo, col)
                        push!(V_coo, Δu[β])
                    end
                end
                
                # Contribution to F_j: -(I₃ ⊗ Δu^T) · φ_ij = -Φ_ij^T · Δu
                # For α ∈ {1,2,3}: F_j[α] = -Σ_β Φ_ij[β,α] * Δu[β] = -Σ_β Φ_ij^T[α,β] * Δu[β]
                # With column-major storage of Φ_ij: Φ_ij[β,α] is at index (α-1)*3 + β
                for α in 1:3
                    row = row_j + α
                    for β in 1:3
                        col = col_offset + (α-1)*3 + β - 1  # Φ_ij[β,α] in column-major
                        push!(I_coo, row)
                        push!(J_coo, col)
                        push!(V_coo, -Δu[β])
                    end
                end
            end
        end
    end
    
    # Build sparse design matrix
    A = sparse(I_coo, J_coo, V_coo, n_equations, n_unknowns)
    
    # Build force vector (flatten column-major)
    F = vec(forces)
    
    # Build Hessian symmetry constraint matrix C
    # Constraint: for each atom i, Σ_j Φ_ij must be symmetric
    # This gives 3 constraints per atom: (1,2), (1,3), (2,3) elements must equal their transposes
    # C · φ = 0 where each row of C encodes one constraint
    n_constraints = 3 * na
    
    # Column-major indices for 3×3 matrix elements (1-indexed, offset by col_offset-1)
    # Φ[α,β] is at index (β-1)*3 + α
    constraint_pairs = [(1,2), (1,3), (2,3)]  # (α, β) pairs
    
    C_I = Int[]
    C_J = Int[]
    C_V = Float64[]
    
    for i in 1:na
        for (c_idx, (α, β)) in enumerate(constraint_pairs)
            constraint_row = (i-1)*3 + c_idx
            
            # Index of Φ[α,β] and Φ[β,α] within a 9-element column-major block
            idx_αβ = (β-1)*3 + α  # Φ[α,β]
            idx_βα = (α-1)*3 + β  # Φ[β,α]
            
            # Contribution from j > i (stored in neighbors[i])
            # Φ_ij[α,β] - Φ_ij[β,α] with coefficient +1
            for (k, j) in enumerate(neighbors[i])
                col_offset = pair_to_col[i][k]
                push!(C_I, constraint_row)
                push!(C_J, col_offset + idx_αβ - 1)
                push!(C_V, 1.0)
                push!(C_I, constraint_row)
                push!(C_J, col_offset + idx_βα - 1)
                push!(C_V, -1.0)
            end
            
            # Contribution from j < i (stored in neighbors[j] as Φ_ji)
            # Φ_ij = Φ_ji^T, so Φ_ij[α,β] = Φ_ji[β,α] and Φ_ij[β,α] = Φ_ji[α,β]
            # Contribution: Φ_ji[β,α] - Φ_ji[α,β] with coefficient +1
            for j in 1:(i-1)
                k = findfirst(==(i), neighbors[j])
                if !isnothing(k)
                    col_offset = pair_to_col[j][k]
                    push!(C_I, constraint_row)
                    push!(C_J, col_offset + idx_βα - 1)  # Φ_ji[β,α]
                    push!(C_V, 1.0)
                    push!(C_I, constraint_row)
                    push!(C_J, col_offset + idx_αβ - 1)  # Φ_ji[α,β]
                    push!(C_V, -1.0)
                end
            end
        end
    end
    
    C_mat = sparse(C_I, C_J, C_V, n_constraints, n_unknowns)
    
    # Solve with Hessian symmetry constraint
    AtF = A' * F  # This is just a vector, always cheap
    
    if hard_constraints
        # Lagrange multiplier method via Schur complement (exact constraints)
        # Requires forming AtA explicitly - may OOM for large/dense systems
        verbose && @info "Solving with Lagrange multipliers (Schur complement)..."
        
        AtA = A' * A
        H = AtA + λ * sparse(I, n_unknowns, n_unknowns)
        b = AtF
        
        # Factor H (positive definite)
        H_chol = cholesky(Symmetric(H))
        
        # Compute H⁻¹ b and H⁻¹ C' efficiently
        H_inv_b = H_chol \ b
        H_inv_Ct = H_chol \ Matrix(C_mat')  # Convert to dense for the solve
        
        # Schur complement: S = C H⁻¹ C' (small: n_constraints × n_constraints)
        S = C_mat * H_inv_Ct
        
        # Solve for Lagrange multipliers: S ν = C H⁻¹ b
        rhs_schur = C_mat * H_inv_b
        ν = S \ rhs_schur
        
        # Recover primal: φ = H⁻¹ (b - C' ν)
        φ = H_inv_b - H_inv_Ct * ν
        
        # Report constraint violation (should be ~machine epsilon)
        constraint_violation = norm(C_mat * φ)
        verbose && @info "Constraint violation (should be ~0)" hessian_symmetry_norm=constraint_violation
    else
        # Penalty method (soft constraints) with direct Cholesky solver
        # minimize ||Aφ - F||² + λ||φ||² + μ||Cφ||²
        # (A'A + λI + μC'C) φ = A'F
        verbose && @info "Solving with penalty method (μ=$μ) using direct Cholesky..."
        
        verbose && @info "Computing A'A (may be slow/memory intensive)..."
        AtA = A' * A
        CtC = C_mat' * C_mat
        AtA_reg = AtA + λ * sparse(I, n_unknowns, n_unknowns) + μ * CtC
        
        verbose && @info "Factorizing..." nnz_AtA_reg=nnz(AtA_reg)
        chol = cholesky(Symmetric(AtA_reg))
        φ = chol \ AtF
        
        # Report constraint violation
        constraint_violation = norm(C_mat * φ)
        verbose && @info "Constraint violation" hessian_symmetry_norm=constraint_violation
    end
    
    # Compute fit quality
    residual = A * φ - F
    rmse = sqrt(mean(residual.^2))
    max_err = maximum(abs.(residual))
    verbose && @info "Fit complete" rmse_Hartree_bohr=rmse max_error_Hartree_bohr=max_err
    
    # Unpack φ into SMatrix blocks
    blocks = Vector{Vector{SMatrix{3,3,Float64,9}}}(undef, na)
    for i in 1:na
        n_nbrs = length(neighbors[i])
        blocks[i] = Vector{SMatrix{3,3,Float64,9}}(undef, n_nbrs)
        for k in 1:n_nbrs
            col_offset = pair_to_col[i][k]
            # Extract 9 elements in column-major order
            φ_ij = SMatrix{3,3,Float64,9}(φ[col_offset:col_offset+8])
            blocks[i][k] = φ_ij
        end
    end
    
    return AmorphousIFC2(na, r_cut, neighbors, blocks)
end


"""
    fit_ifc2(crystal, positions_path, forces_path, n_timesteps; kwargs...) -> AmorphousIFC2

Convenience method that reads positions and forces from files.
"""
function fit_ifc2(
    crystal::CrystalStructure,
    positions_path::String,
    forces_path::String,
    n_timesteps::Int;
    kwargs...
)
    na = length(crystal)
    positions = read_positions(positions_path, crystal, n_timesteps)
    forces = read_forces(forces_path, na, n_timesteps)
    return fit_ifc2(crystal, positions, forces; kwargs...)
end

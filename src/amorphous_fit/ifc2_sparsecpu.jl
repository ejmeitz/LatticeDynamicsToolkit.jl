function _assemble_ifc2_design(problem::AmorphousFitProblem)
    crystal = problem.crystal
    na = length(crystal)
    ns = n_samples(problem)
    neighbors = _build_upper_neighbor_list(crystal, problem.r_cut)
    pair_to_col, n_unknowns = _build_pair_to_col(neighbors)
    n_pairs = n_unknowns ÷ 9
    n_equations = 3 * na * ns

    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, 18 * n_pairs * ns)
    sizehint!(J, 18 * n_pairs * ns)
    sizehint!(V, 18 * n_pairs * ns)

    F = Vector{Float64}(undef, n_equations)
    x0 = crystal.x_cart
    u_t = Vector{SVector{3,Float64}}(undef, na)

    @inbounds for s in 1:ns
        for a in 1:na
            u_t[a] = _displacement(crystal, _sample_position(problem, s, a), x0[a])
            for α in 1:3
                row = (s - 1) * 3 * na + 3 * (a - 1) + α
                F[row] = _sample_force(problem, s, a, α)
            end
        end

        for i in 1:na
            for (k, j) in enumerate(neighbors[i])
                Δu = u_t[i] - u_t[j]
                col_offset = pair_to_col[i][k]
                row_i = (s - 1) * 3 * na + 3 * (i - 1)
                row_j = (s - 1) * 3 * na + 3 * (j - 1)

                for α in 1:3
                    row = row_i + α
                    for β in 1:3
                        col = col_offset + (β - 1) * 3 + α - 1
                        push!(I, row)
                        push!(J, col)
                        push!(V, Δu[β])
                    end
                end

                for α in 1:3
                    row = row_j + α
                    for β in 1:3
                        col = col_offset + (α - 1) * 3 + β - 1
                        push!(I, row)
                        push!(J, col)
                        push!(V, -Δu[β])
                    end
                end
            end
        end
    end

    A = sparse(I, J, V, n_equations, n_unknowns)
    return A, F, neighbors, pair_to_col
end

function _assemble_ifc2_constraints(neighbors::Vector{Vector{Int}}, pair_to_col::Vector{Vector{Int}}, n_unknowns::Int)
    na = length(neighbors)
    n_constraints = 3 * na
    constraint_pairs = ((1, 2), (1, 3), (2, 3))
    I = Int[]
    J = Int[]
    V = Float64[]

    @inbounds for i in 1:na
        for (c_idx, (α, β)) in enumerate(constraint_pairs)
            row = (i - 1) * 3 + c_idx
            idx_αβ = (β - 1) * 3 + α
            idx_βα = (α - 1) * 3 + β

            for (k, _) in enumerate(neighbors[i])
                col_offset = pair_to_col[i][k]
                push!(I, row); push!(J, col_offset + idx_αβ - 1); push!(V, 1.0)
                push!(I, row); push!(J, col_offset + idx_βα - 1); push!(V, -1.0)
            end

            for j in 1:(i - 1)
                k = findfirst(==(i), neighbors[j])
                if k !== nothing
                    col_offset = pair_to_col[j][k]
                    push!(I, row); push!(J, col_offset + idx_βα - 1); push!(V, 1.0)
                    push!(I, row); push!(J, col_offset + idx_αβ - 1); push!(V, -1.0)
                end
            end
        end
    end

    return sparse(I, J, V, n_constraints, n_unknowns)
end

function _solve_ifc2(::LSMRAlgorithm, A, b, options::AmorphousFitOptions)
    tol = max(options.atol, options.rtol * norm(b))
    alg = LSMR(; maxiter=options.maxiter, tol=tol, verbosity=options.verbosity)
    return lssolve(A, b, alg, options.regularization)
end


_solve_ifc2(::KrylovAuto, A, b, options::AmorphousFitOptions) =
    _solve_ifc2(LSMRAlgorithm(), A, b, options)

function _fit_ifc2_stage(
    problem::AmorphousFitProblem,
    ::SparseCPU,
    constraints::ApproximateConstraints,
    ::ZeroInitialGuess,
    ::LeastSquaresFormulation,
    algorithm::Union{KrylovAuto,LSMRAlgorithm},
    options::AmorphousFitOptions,
)
    A, F, neighbors, pair_to_col = _assemble_ifc2_design(problem)
    C = _assemble_ifc2_constraints(neighbors, pair_to_col, size(A, 2))
    A_aug = vcat(A, sqrt(constraints.weight) * C)
    b_aug = vcat(F, zeros(Float64, size(C, 1)))

    φ, stats = if norm(b_aug) == 0
        zeros(Float64, size(A_aug, 2)), nothing
    else
        _solve_ifc2(algorithm, A_aug, b_aug, options)
    end
    ifc2 = _pack_ifc2(φ, problem.r_cut, neighbors, pair_to_col)

    residual = A * φ - F
    constraint_residual = C * φ
    messages = if stats === nothing
        ["Zero right-hand side: returned zero IFCs without Krylov iteration."]
    else
        [
            "KrylovKit converged: $(stats.converged == 1)",
            "KrylovKit iterations: $(stats.numiter)",
            "KrylovKit operations: $(stats.numops)",
            "KrylovKit normal residual: $(stats.normres)",
        ]
    end
    diagnostics = FitDiagnostics(
        force_rmse=sqrt(mean(residual .^ 2)),
        force_max_abs=maximum(abs.(residual)),
        constraint_norm=norm(constraint_residual),
        messages=messages,
    )

    return IFC2FitResult(ifc2, diagnostics)
end

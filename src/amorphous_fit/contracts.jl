export
    validate_problem,
    n_samples,
    n_atoms

n_atoms(problem::AmorphousFitProblem) = length(problem.crystal)
n_samples(problem::AmorphousFitProblem) = size(problem.positions, 1)

function validate_problem(problem::AmorphousFitProblem)
    na = n_atoms(problem)
    ndof = 3 * na

    size(problem.positions, 2) == ndof || throw(ArgumentError(
        "positions must have shape (n_samples, 3 * na); got $(size(problem.positions)) for na=$na."
    ))
    size(problem.forces) == size(problem.positions) || throw(ArgumentError(
        "forces must have the same shape as positions; got forces=$(size(problem.forces)), positions=$(size(problem.positions))."
    ))
    problem.r_cut > 0 || throw(ArgumentError("r_cut must be positive, got $(problem.r_cut)."))

    return problem
end

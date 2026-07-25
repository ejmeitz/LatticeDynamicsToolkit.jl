export
    fit_amorphous,
    fit_ifc2_stage,
    fit_ifc3_residual_stage

function _stub_error(stage::AbstractString)
    throw(ArgumentError(
        "$stage is scaffolded but not implemented yet. " *
        "This pass only defines the amorphous fitting contracts and dispatch surface."
    ))
end

function fit_amorphous(problem::AmorphousFitProblem, options::AmorphousFitOptions = AmorphousFitOptions())
    validate_problem(problem)
    ifc2_result = fit_ifc2_stage(problem, options)
    return (; ifc2 = ifc2_result, ifc3 = nothing)
end

function fit_ifc2_stage(problem::AmorphousFitProblem, options::AmorphousFitOptions)
    validate_problem(problem)
    validate_options(options)
    return _fit_ifc2_stage(problem, options.storage, options.constraints,
        options.initial_guess, options.formulation, options.algorithm, options)
end

function fit_ifc3_residual_stage(
    problem::AmorphousFitProblem,
    ifc2_result::IFC2FitResult,
    options::AmorphousFitOptions,
)
    validate_problem(problem)
    validate_options(options)
    return _fit_ifc3_residual_stage(problem, ifc2_result, options)
end

function _fit_ifc2_stage(
    problem::AmorphousFitProblem,
    storage::DesignMatrixStorage,
    constraints::ConstraintEnforcement,
    initial_guess::InitialGuessStrategy,
    formulation::KrylovProblemFormulation,
    algorithm::KrylovAlgorithmChoice,
    options::AmorphousFitOptions,
)
    _stub_error("IFC2 fitting with $(typeof(storage)), $(typeof(constraints)), $(typeof(initial_guess)), $(typeof(formulation)), and $(typeof(algorithm))")
end

function _fit_ifc3_residual_stage(
    problem::AmorphousFitProblem,
    ifc2_result::IFC2FitResult,
    options::AmorphousFitOptions,
)
    _stub_error("Residual IFC3 fitting")
end

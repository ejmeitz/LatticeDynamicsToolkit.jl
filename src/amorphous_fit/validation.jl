export validate_options

function validate_options(options::AmorphousFitOptions)
    options.regularization >= 0 || throw(ArgumentError(
        "regularization must be non-negative, got $(options.regularization)."
    ))
    if options.constraints isa ApproximateConstraints
        options.constraints.weight >= 0 || throw(ArgumentError(
            "Approximate constraint weight must be non-negative, got $(options.constraints.weight)."
        ))
    end
    options.maxiter > 0 || throw(ArgumentError("maxiter must be positive, got $(options.maxiter)."))
    options.rtol >= 0 || throw(ArgumentError("rtol must be non-negative, got $(options.rtol)."))
    options.atol >= 0 || throw(ArgumentError("atol must be non-negative, got $(options.atol)."))
    return options
end

function validate_options(options::AmorphousFitOptions{S,C,G,PriorRegularizedFormulation,A}) where {S,C,G,A}
    validate_options(AmorphousFitOptions(
        options.storage,
        options.constraints,
        options.initial_guess,
        LeastSquaresFormulation(),
        options.algorithm,
        options.regularization,
        options.maxiter,
        options.rtol,
        options.atol,
        options.verbosity,
    ))
    options.initial_guess isa ProvidedIFC2Guess || throw(ArgumentError(
        "PriorRegularizedFormulation requires ProvidedIFC2Guess so the prior is explicit."
    ))
    return options
end

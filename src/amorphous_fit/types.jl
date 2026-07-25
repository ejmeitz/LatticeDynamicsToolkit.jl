export
    AmorphousFitProblem,
    AmorphousFitOptions,
    FitDiagnostics,
    IFC2FitResult,
    IFC3ResidualFitResult,
    DesignMatrixStorage,
    MatrixFree,
    SparseCPU,
    SparseGPU,
    ConstraintEnforcement,
    ApproximateConstraints,
    ExactConstraints,
    InitialGuessStrategy,
    ZeroInitialGuess,
    ProvidedIFC2Guess,
    KrylovProblemFormulation,
    LeastSquaresFormulation,
    PriorRegularizedFormulation,
    KrylovAlgorithmChoice,
    KrylovAuto,
    LSQRAlgorithm,
    LSMRAlgorithm

abstract type DesignMatrixStorage end
struct MatrixFree <: DesignMatrixStorage end
struct SparseCPU <: DesignMatrixStorage end
struct SparseGPU <: DesignMatrixStorage end

abstract type ConstraintEnforcement end
Base.@kwdef struct ApproximateConstraints <: ConstraintEnforcement
    weight::Float64 = 1e2
end
struct ExactConstraints <: ConstraintEnforcement end

abstract type InitialGuessStrategy end
struct ZeroInitialGuess <: InitialGuessStrategy end
struct ProvidedIFC2Guess <: InitialGuessStrategy
    ifc2::AmorphousIFC2
end

abstract type KrylovProblemFormulation end
struct LeastSquaresFormulation <: KrylovProblemFormulation end
struct PriorRegularizedFormulation <: KrylovProblemFormulation end

abstract type KrylovAlgorithmChoice end
struct KrylovAuto <: KrylovAlgorithmChoice end
struct LSQRAlgorithm <: KrylovAlgorithmChoice end
struct LSMRAlgorithm <: KrylovAlgorithmChoice end

struct AmorphousFitProblem{TP<:AbstractMatrix{Float64},TF<:AbstractMatrix{Float64}}
    crystal::CrystalStructure
    positions::TP
    forces::TF
    r_cut::Float64
end

struct AmorphousFitOptions{
    S<:DesignMatrixStorage,
    C<:ConstraintEnforcement,
    G<:InitialGuessStrategy,
    F<:KrylovProblemFormulation,
    A<:KrylovAlgorithmChoice,
}
    storage::S
    constraints::C
    initial_guess::G
    formulation::F
    algorithm::A
    regularization::Float64
    maxiter::Int
    rtol::Float64
    atol::Float64
    verbosity::Int
end

function AmorphousFitOptions(;
    storage::DesignMatrixStorage = SparseCPU(),
    constraints::ConstraintEnforcement = ApproximateConstraints(),
    initial_guess::InitialGuessStrategy = ZeroInitialGuess(),
    formulation::KrylovProblemFormulation = LeastSquaresFormulation(),
    algorithm::KrylovAlgorithmChoice = LSMRAlgorithm(),
    regularization::Float64 = 1e-6,
    maxiter::Int = 2000,
    rtol::Float64 = 1e-4,
    atol::Float64 = 0.0,
    verbosity::Int = KrylovKit.STARTSTOP_LEVEL,
)
    return AmorphousFitOptions(
        storage,
        constraints,
        initial_guess,
        formulation,
        algorithm,
        regularization,
        maxiter,
        rtol,
        atol,
        verbosity,
    )
end

Base.@kwdef struct FitDiagnostics
    force_rmse::Union{Nothing,Float64} = nothing
    force_max_abs::Union{Nothing,Float64} = nothing
    constraint_norm::Union{Nothing,Float64} = nothing
    messages::Vector{String} = String[]
end

struct IFC2FitResult
    ifc2::Union{Nothing,AmorphousIFC2}
    diagnostics::FitDiagnostics
end

struct IFC3ResidualFitResult
    ifc3::Union{Nothing,IFC3}
    residual_forces::Union{Nothing,Matrix{Float64}}
    diagnostics::FitDiagnostics
end

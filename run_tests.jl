# Convenience script to run tests: julia run_tests.jl
# Alternative: julia -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")'
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
include(joinpath(@__DIR__, "test", "runtests.jl"))

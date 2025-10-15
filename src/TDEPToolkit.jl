module TDEPToolkit

using LinearAlgebra
using StaticArrays
using CellListMap
using OhMyThreads
using PrecompileTools

const lo_tol = 1e-5
const lo_sqtol = lo_tol^2

include("util.jl")
include("types.jl")
include("distance_table.jl")
include("io.jl")
include("remap.jl")
include("tep.jl")


##! ADD PRECOMPILATION OF IFC LOADING AND REMAPPING

end # module TDEPToolkit
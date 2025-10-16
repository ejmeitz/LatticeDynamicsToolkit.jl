module TDEPToolkit

using LinearAlgebra
using StaticArrays
using CellListMap
using OhMyThreads

const lo_tol = 1e-5
const lo_sqtol = lo_tol^2

include("util.jl")
include("types.jl")
include("distance_table.jl")
include("io.jl")
include("remap.jl")
include("tep.jl")


precompile(read_ifc2, (String, String))
precompile(read_ifc3, (String, String))
precompile(read_ifc4, (String, String))

precompile(CrystalStructure, (String,))

#! precompile energies and remapping?

end # module TDEPToolkit
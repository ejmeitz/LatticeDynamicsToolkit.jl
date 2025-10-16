module TDEPToolkit

using LinearAlgebra
using StaticArrays
using CellListMap
using OhMyThreads
import PeriodicTable

const periodic_table = PeriodicTable.elements

const lo_tol = 1e-5
const lo_sqtol = lo_tol^2
const kB = 8.617333262e-5 # eV / K
const hbar = 6.582119569e-16 # eV * s

include("util.jl")
include("types.jl")
include("modes.jl")
include("distance_table.jl")
include("io.jl")
include("remap.jl")
include("epot.jl")


precompile(read_ifc2, (String, String))
precompile(read_ifc3, (String, String))
precompile(read_ifc4, (String, String))

precompile(CrystalStructure, (String,))

#! precompile energies and remapping?

end # module TDEPToolkit
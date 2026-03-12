module LatticeDynamicsToolkit

using AtomsBase
using AtomsCalculators
using LinearAlgebra
using SparseArrays
using StaticArrays
using CellListMap
using OhMyThreads
import Spglib
import Unitful: ustrip
import PeriodicTable
using ProgressMeter
using Statistics: mean
import Random: randn!
using HDF5
using KrylovKit
using Roots
import SpecialFunctions: erfc
using TDEPWrapper

include("constants.jl")
include("util.jl")

include("io/io_utils.jl")
include("io/poscar.jl")
include("types/types.jl")
include("types/stubs.jl")
include("types/meshes.jl")
include("types/ifcs.jl")
include("io/forceconstants.jl")
include("io/tdep_io.jl")
include("io/hdf5_io.jl")


include("ewald.jl")
include("harmonic/dynmat_longrange.jl")
include("harmonic/dynmat.jl")
include("harmonic/thermo.jl")
include("harmonic/dispersion.jl")
include("harmonic/canonical_configs.jl")
include("anharmonic/free_energy.jl")

include("distance_table.jl")
include("remap.jl")
include("epot.jl")
include("fit_amorphous.jl")


precompile(read_ifc2, (String, String))
precompile(read_ifc3, (String, String))
precompile(read_ifc4, (String, String))

precompile(CrystalStructure, (String,))

#! precompile energies and remapping?

end # module LatticeDynamicsToolkit
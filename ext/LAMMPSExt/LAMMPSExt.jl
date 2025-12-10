module LAMMPSExt

using LAMMPS
using LatticeDynamicsToolkit
import LatticeDynamicsToolkit: bohr_to_A, A_to_bohr, Hartree_to_eV, emu_to_amu, lo_tol, forceconstant_2nd_HartreeBohr_to_eVA
using AtomsBase
using AtomsCalculators
using LinearAlgebra
using Unitful
using OhMyThreads
using ProgressMeter
using CellListMap
using StaticArrays
import Random: randn!

include("lammps_calc.jl")
include("energy_dataset.jl")

include("gauss_legendre.jl")
include("TI.jl")

end
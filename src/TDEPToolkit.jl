module TDEPToolkit

using AtomsBase
using AtomsCalculators
using LinearAlgebra
using StaticArrays
using CellListMap
using OhMyThreads
import Unitful: ustrip
import PeriodicTable
using ProgressMeter
import Random: randn!

const periodic_table = PeriodicTable.elements

const lo_tol = 1e-5
const lo_sqtol = lo_tol^2

# Hartrees to eV
const Hartree_to_eV = 27.21138602
const eV_to_Hartree = 1.0 / Hartree_to_eV

# Bohr radius to angstrom
const bohr_to_A = 0.52917721067
const A_to_bohr = 1.0 / bohr_to_A

# Atomic and electron mass unit
const amu_to_kg = 1.660539040E-27
const emu_to_kg = 9.10938356E-31
const amu_to_emu = amu_to_kg/emu_to_kg
const emu_to_amu = 1.0/amu_to_emu

const kB_eV = 8.6173303E-5
const kB_Hartree = kB_eV*eV_to_Hartree

const hbar_eV = 6.582119514E-16
const hbar_Hartree = hbar_eV*eV_to_Hartree

const forceconstant_2nd_eVA_to_HartreeBohr = eV_to_Hartree / (A_to_bohr^2)
const forceconstant_3rd_eVA_to_HartreeBohr = eV_to_Hartree / (A_to_bohr^3)
const forceconstant_4th_eVA_to_HartreeBohr = eV_to_Hartree / (A_to_bohr^4)

const forceconstant_2nd_HartreeBohr_to_eVA = 1.0 / forceconstant_2nd_eVA_to_HartreeBohr
const forceconstant_3rd_HartreeBohr_to_eVA = 1.0 / forceconstant_3rd_eVA_to_HartreeBohr
const forceconstant_4th_HartreeBohr_to_eVA = 1.0 / forceconstant_4th_eVA_to_HartreeBohr

const frequency_Hartree_to_THz=1e-12/(2*pi)/hbar_Hartree
const frequency_THz_to_Hartree=1.0/frequency_Hartree_to_THz



include("util.jl")
include("types.jl")
include("modes.jl")
include("distance_table.jl")
include("io.jl")
include("remap.jl")
include("canonical_configs.jl")
include("epot.jl")


precompile(read_ifc2, (String, String))
precompile(read_ifc3, (String, String))
precompile(read_ifc4, (String, String))

precompile(CrystalStructure, (String,))

#! precompile energies and remapping?

end # module TDEPToolkit
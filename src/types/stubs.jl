export  LAMMPSCalculator, ThermodynmicIntegration, TISettings

"""
    LAMMPSCalculator(
        sys::CrystalStructure,
        potential_definition::Union{String, Array{String}};
        label_type_map::Dict{Symbol, Int} = Dict{Symbol, Int}(),
        logfile_path::String = "none",
    )

Defines a general interaction that will call LAMMPS to calculate forces and energies. Forces
and energies are calculated on a single thread. You must call LAMMPS.MPI.Init() for LAMMPS.jl
to load the LAMMPS executable on systems where MPI is available. To speed-up single point
calculations the neighbor list skin distance is set to the nearest-neighbor distance in the system
and never rebuilt.

The LAMMPS potential files can be found at:
`abspath(dirname(LAMMPS.locate()), "..", "share", "lammps", "potentials")`

Restrictions:
-------------
- CPU only
- Floats promote to Float64
- No triclinic boundary
- 3D systems only
- Fully periodic systems only
- Expects 'metal' unit system
- Only crystals

Arguments:
----------
- `sys::CrystalStructure`: The system object this interaction will be applied to. 
- `potential_definition::Union{String, Array{String}}` : Commands passed to lammps which define your interaction.
    For example, to define LJ you pass:
    `lj_cmds = ["pair_style lj/cut 8.5", "pair_coeff * * 0.0104 3.4", "pair_modify shift yes"]`
- `label_type_map::Dict{Symbol, Int} = Dict{Symbol, Int}()` : By default atom types are assigned in the 
    order they appear in the system. This can make defining the potential for multi-atomic systems
    difficult. By providing this dictionary you can overide the type label assigned to each unique species. 
- `logfile_path::String = "none"` : Path where LAMMPS logfile is written. Defaults to no log file. 
"""
mutable struct LAMMPSCalculator{T, S}
    lmp::T # T will be LMP but that is not available here
    pot_cmds::S # S will be String or Vector{String}
end


"""
TODO
"""
struct TISettings 
    T::Float64
    pot_cmds::Union{String, Vector{String}}
    nsteps::Integer
    nsteps_equil::Integer
    n_lambda::Integer
    dt_ps::Float64
    T_damp::Float64
end

function TISettings(T, pot_cmds, nsteps, nsteps_equil; n_lambda = 9, dt_ps = 1e-3, T_damp = 100*dt_ps)
    return TISettings(T, pot_cmds, nsteps, nsteps_equil, n_lambda, dt_ps, T_damp)
end

"""
TODO
"""
function ThermodynmicIntegration end
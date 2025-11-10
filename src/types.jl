export 
    IFC2, IFC3, IFC4,
    CrystalStructure, ConfigSettings,
    ClassicalConfigSettings, QuantumConfigSettings,
    LAMMPSCalculator, ThermodynmicIntegration, TISettings,
    Quantum, Classical, Limit

abstract type IFCs end

# mimics the lo_fc2_pair type
struct FC2Data 
    idxs::SVector{2, Int} # indices in unitcell
    lvs::SVector{2, SVector{3, Float64}} # lattice vectors for the unit cell these atoms belong to
    r::SVector{3, Float64} # vector between atom 1 and 2
    n2::SVector{3,Int16} # image flags for each atom, n1 is always zero vector
    ifcs::SMatrix{3, 3, Float64, 9}
end 

# mimics the lo_fc3_triplet type
struct FC3Data
    idxs::SVector{3, Int} # indices in unitcell
    lvs::SVector{3, SVector{3, Float64}}
    rv1::SVector{3, Float64}
    rv2::SVector{3, Float64}
    rv3::SVector{3, Float64}
    n2::SVector{3,Int16}
    n3::SVector{3,Int16}
    ifcs::SArray{Tuple{3,3,3}, Float64, 3, 27}
end

struct FC4Data
    idxs::SVector{4, Int} # indices in unitcell
    lvs::SVector{4, SVector{3, Float64}}
    rv1::SVector{3, Float64}
    rv2::SVector{3, Float64}
    rv3::SVector{3, Float64}
    rv4::SVector{3, Float64}
    n2::SVector{3,Int16}
    n3::SVector{3,Int16}
    n4::SVector{3,Int16}
    ifcs::SArray{Tuple{3,3,3,3}, Float64, 4, 81}
end

struct IFC2 <: IFCs
    na::Int # number of atoms in the cell
    r_cut::Float64
    # all interactions for each atom in the unitcell, has length na_uc
    # by translational symmetry we can re-build all other IFCs
    atoms::Vector{Vector{FC2Data}}
end

struct IFC3 <: IFCs
    na::Int
    r_cut::Float64
    atoms::Vector{Vector{FC3Data}}
end

struct IFC4 <: IFCs
    na::Int
    r_cut::Float64
    atoms::Vector{Vector{FC4Data}}
end

get_interactions(data::I, i::Int) where {I <: IFCs} = data.atoms[i]

get_kwarg(::IFC2) = :ifc2
get_kwarg(::IFC3) = :ifc3
get_kwarg(::IFC4) = :ifc4

build_kwargs(ifcs::IFCs...) = NamedTuple{get_kwarg.(ifcs)}(ifcs)

Base.show(io::IO, ifc::IFC2) =
    print(io, "2nd Order IFCs, cutoff = $(round(ifc.r_cut*bohr_to_A, digits = 5)) Ang from $(ifc.na) atom unit-cell")

Base.show(io::IO, ifc::IFC3) =
    print(io, "3rd Order IFCs, cutoff = $(round(ifc.r_cut*bohr_to_A, digits = 5)) Ang from $(ifc.na) atom unit-cell")

Base.show(io::IO, ifc::IFC4) =
    print(io, "4th Order IFCs, cutoff = $(round(ifc.r_cut*bohr_to_A, digits = 5)) Ang from $(ifc.na) atom unit-cell")

#######################

struct DistanceTableAtom{N}
    central_atom::Int
    # Vector from current center atom to other particle
    vs::SVector{N,SVector{3,Float64}}
    # Vector from current cetner atom to other particle with pbc accounted for
    lvs::SVector{N,SVector{3,Float64}}
    # image flags
    ns::SVector{N, SVector{3, Int16}}
    inds::SVector{N, Int}
    dists::SVector{N, Float64}
end

n_neighbors(::DistanceTableAtom{N}) where N = N

struct DistanceTable{D <: DistanceTableAtom}
    atoms::Vector{D}
end


##############################

struct CrystalStructure <: AbstractSystem{3}
    x_frac::Vector{SVector{3,Float64}}
    x_cart::Vector{SVector{3,Float64}} # bohr
    species::Vector{Symbol}
    m::Vector{Float64} # in emu
    invsqrtm::Vector{Float64} # 1/sqrt(m)
    L::SMatrix{3,3,Float64,9} # bohr
    L_inv::SMatrix{3,3,Float64,9}
    sym_data::Spglib.Dataset
end


"""
    CrystalStructure(poscar_path::String)

- `poscar_path::String` : Path to POSCAR file to parse, only fractional coords supported
"""
function CrystalStructure(poscar_path::String; symprec::Float64 = 1e-5) 

    species, x_frac, cell = read_poscar_data(poscar_path)

    # convert lattice vectors to bohr
    # will make all coordinates in bohr as well
    cell *= A_to_bohr
    cell_inv = inv(cell)

    ls = length(species); lc = length(x_frac)
    if length(species) != length(x_frac)
        error("Length of species vectors built from POSCAR ($(ls)) does not match number of atoms in POSCAR ($(lc)). Weird.")
    end

    x_cart = to_cart_coords.(Ref(cell), x_frac)

    # convert mass to electron mass units (emu)
    m = ustrip.([periodic_table[s].atomic_mass for s in species]) .* amu_to_emu 
    invsqrtm = 1.0 ./ sqrt.(m)

    spglib_cell = Spglib.SpglibCell(
        cell,
        x_frac,
        atomic_number.(species)
    )

    sym_data = Spglib.get_dataset(spglib_cell, symprec)

    return CrystalStructure(x_frac, x_cart, species, m, invsqrtm, cell, cell_inv, sym_data)
end


volume(cs::CrystalStructure) = det(cs.L)
primitive_volume(cs::CrystalStructure) = det(cs.sym_data.primitive_lattice)
n_atoms_primitive(cs::CrystalStructure) = length(unique(cs.sym_data.mapping_to_primitive))

Base.length(sys::CrystalStructure) = length(sys.x_frac)
Base.size(sys::CrystalStructure) = size(sys.x_frac)

Base.getindex(sys::CrystalStructure, i::Integer) = AtomView(sys, i)

AtomsBase.cell_vectors(sys::CrystalStructure) = AtomsBase._auto_cell_vectors(tuple(eachcol(sys.L)...))
AtomsBase.periodicity(::CrystalStructure) = (true, true, true)

function Base.getindex(system::CrystalStructure, x::Symbol)
    if x === :cell_vectors
        cell_vectors(system)
    elseif x === :periodicity
        periodicity(system)
    else
        throw(KeyError(x))
    end
end
Base.haskey(::CrystalStructure, x::Symbol) = x in (:cell_vectors, :periodicity)
Base.keys(::CrystalStructure) = (:cell_vectors, :periodicity)

# Atom and atom property access
AtomsBase.atomkeys(::CrystalStructure) = (:x_frac, :x_cart, :m, :invsqrtm)
AtomsBase.cell(sys::CrystalStructure) = AtomsBase.PeriodicCell(cell_vectors(sys), periodicity(sys))

AtomsBase.hasatomkey(system::CrystalStructure, x::Symbol) = x in atomkeys(system)

function Base.getindex(system::CrystalStructure, i::Union{Integer,AbstractVector}, x::Symbol)
    getfield(system, x)[i]
end

Base.getindex(system::CrystalStructure, ::Colon, x::Symbol) = getfield(system, x)

AtomsBase.position(s::CrystalStructure, ::Colon) = s.x_cart
AtomsBase.position(sys::CrystalStructure, i::Union{Integer, AbstractVector}) = sys.x_cart[i]

AtomsBase.mass(s::CrystalStructure, ::Colon) = s.m
AtomsBase.mass(sys::CrystalStructure, i::Union{Integer, AbstractVector}) = sys.m[i]

AtomsBase.species(s::CrystalStructure, ::Colon) = AtomsBase.ChemicalSpecies.(s.species)
AtomsBase.species(sys::CrystalStructure, i::Union{Integer, AbstractVector}) = AtomsBase.ChemicalSpecies(sys.species[i])

AtomsBase.atomic_symbol(s::CrystalStructure, ::Colon) = s.species
AtomsBase.atomic_symbol(s::CrystalStructure, i::Union{Integer, AbstractVector}) = s.species[i]

AtomsBase.atomic_number(s::CrystalStructure, ::Colon) = atomic_number.(s.species)
AtomsBase.atomic_number(s::CrystalStructure, i::Union{Integer, AbstractVector}) = atomic_number(s.species[i])

######################

# Gamma centered IBZ mesh
struct IBZMesh{I <: Integer}
    mesh::SVector{3, I}
    k_ibz::Vector{SVector{3, Float64}}
    weights::Vector{Float64}
    radius::Float64 # no integration weight included
    n_atoms_prim::Integer
end


n_full_q_point(ibz::IBZMesh) = prod(ibz.mesh)
Base.length(ibz::IBZMesh) = length(ibz.k_ibz)
Base.eachindex(ibz::IBZMesh) = eachindex(ibz.weights)

# Add IBZ constructor from Spglib using the types defined here
function Spglib.BrillouinZoneMesh(uc::CrystalStructure, mesh; symprec = 1e-5)

    cell = Spglib.SpglibCell(
        uc.L,
        uc.x_frac,
        Int.(AtomsBase.atomic_number(uc, :))
    )

    return Spglib.get_ir_reciprocal_mesh(cell, mesh, symprec)
end

function IBZMesh(uc::CrystalStructure, mesh; symprec = 1e-5)

    bzm = Spglib.BrillouinZoneMesh(uc, mesh; symprec = symprec)
    
    N = length(bzm.grid_address)

    # Indices of irreducible k-points (sorted for stable downstream use)
    ir_idx = sort!(unique(bzm.ir_mapping_table))

    if length(ir_idx) > 0.5 * length(bzm.grid_address)
        @warn "The number of irreducible k-points ($(length(ir_idx))) is more than half the total number of k-points ($(length(bzm.grid_address)))." *
              " This may indicate that the symmetry detection failed. Consider increasing symprec (currently $(symprec))."
    end

    counts = zeros(Int, N)
    @inbounds for i in 1:N
        j = bzm.ir_mapping_table[i]
        counts[j] += 1
    end
    multiplicity = counts[ir_idx]

    # Get list of fractional IBZ points
    invmesh = 1.0 ./ Float64.(bzm.mesh)
    k_ibz = Vector{SVector{3,Float64}}(undef, length(ir_idx))
    @inbounds for (p, j) in pairs(ir_idx)
        k_ibz[p] = (Float64.(bzm.grid_address[j])) .* invmesh
    end

    weights = multiplicity ./ sum(multiplicity)

    @assert length(k_ibz) == length(weights) "something wrong in IBZMesh construction"

    # cbrt(3 / (4*pi*V*Nk))
    radius = cbrt(3.0/primitive_volume(uc)/prod(mesh)/4.0/pi)

    N_prim = n_atoms_primitive(uc)

    return IBZMesh(bzm.mesh, k_ibz, weights, radius, N_prim)
end

######################

struct DispersionData{N} # N = N branch
    ibz::IBZMesh
    freqs::Vector{SVector{N, Float64}}  # each entry sorted smallest --> largest
    vels::Vector{SMatrix{3, N, Float64}} # Length is Nq
end 

# Constructor in dispersion.jl

######################

abstract type Limit end
struct Quantum <: Limit end
struct Classical <: Limit end

struct ConfigSettings{L <: Limit}
    n_configs::Int
    temperature::Float64
end

function ConfigSettings(n_configs::Int, temperature::Float64, ::Type{L}) where {L <: Limit}
    return ConfigSettings{L}(n_configs, temperature)
end

# Type aliases
const QuantumConfigSettings = ConfigSettings{Quantum}
const ClassicalConfigSettings = ConfigSettings{Classical}

######################

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
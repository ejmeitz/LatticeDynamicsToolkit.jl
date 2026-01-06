export 
    CrystalStructure, ConfigSettings,
    ClassicalConfigSettings, QuantumConfigSettings,
    Quantum, Classical, Limit
    

#######################
# Distance Table Types
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
# Crystal Structure
##############################

struct CrystalStructure <: AbstractSystem{3}
    x_frac::Vector{SVector{3,Float64}}
    x_cart::Vector{SVector{3,Float64}} # bohr
    species::Vector{Symbol}
    m::Vector{Float64} # in emu
    invsqrtm::Vector{Float64} # 1/sqrt(m)
    L::SMatrix{3,3,Float64,9} # bohr
    L_inv::SMatrix{3,3,Float64,9}
    sym_data::Union{Spglib.Dataset, Nothing}
end


"""
    CrystalStructure(poscar_path::String; symprec=1e-5, compute_symmetry=true)

- `poscar_path::String` : Path to POSCAR file to parse, only fractional coords supported
- `symprec::Float64` : Symmetry precision for Spglib (default 1e-5)
- `compute_symmetry::Bool` : Whether to compute space group symmetry (default true).
   Set to `false` for amorphous systems to skip Spglib analysis.
"""
function CrystalStructure(poscar_path::String; symprec::Float64 = 1e-5, compute_symmetry::Bool = true) 

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

    sym_data = nothing
    if compute_symmetry
        spglib_cell = Spglib.SpglibCell(
            cell,
            x_frac,
            atomic_number.(species)
        )
        sym_data = Spglib.get_dataset(spglib_cell, symprec)
    end

    return CrystalStructure(x_frac, x_cart, species, m, invsqrtm, cell, cell_inv, sym_data)
end


volume(cs::CrystalStructure) = det(cs.L)

function primitive_volume(cs::CrystalStructure)
    isnothing(cs.sym_data) && error("Symmetry data not available. Set compute_symmetry=true when constructing CrystalStructure.")
    return det(cs.sym_data.primitive_lattice)
end

function n_atoms_primitive(cs::CrystalStructure)
    isnothing(cs.sym_data) && error("Symmetry data not available. Set compute_symmetry=true when constructing CrystalStructure.")
    return length(unique(cs.sym_data.mapping_to_primitive))
end

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

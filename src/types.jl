export IFC2, IFC3, IFC4, CrystalStructure, ClassicalConfigSettings, QuantumConfigSettings

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

# Could make this an AtomsBase compatible type
struct CrystalStructure
    x_frac::Vector{SVector{3,Float64}}
    x_cart::Vector{SVector{3,Float64}} # bohr
    species::Vector{Symbol}
    m::Vector{Float64} # in emu
    invsqrtm::Vector{Float64} # 1/sqrt(m)
    L::SMatrix{3,3,Float64,9} # bohr
    L_inv::SMatrix{3,3,Float64,9} 
end

Base.length(cs::CrystalStructure) = length(cs.x_frac)

"""

Parameters:
-----------
- `poscar_path::String` : Path to POSCAR file to parse, only fractional coords supported
"""
function CrystalStructure(poscar_path::String) 

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

    return CrystalStructure(x_frac, x_cart, species, m, invsqrtm, cell, cell_inv)
end


###################

abstract type ConfigSettings end

struct QuantumConfigSettings <: ConfigSettings
    n_configs::Int
    temperature::Float64
end

struct ClassicalConfigSettings <: ConfigSettings
    n_configs::Int
    temperature::Float64
end
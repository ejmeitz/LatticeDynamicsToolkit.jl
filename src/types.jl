export IFCs, CrystalStructure

abstract type FCData{O,T} end
abstract type AtomFC{O,T,N} end

# mimics the lo_fc2_pair type
struct FC2Data{T} <: FCData{2,T}
    idxs::SVector{2, Int} # indices in unitcell
    lvs::SVector{2, SVector{3, T}} # lattice vectors for the unit cell these atoms belong to
    r::SVector{3, T} # vector between atom 1 and 2
    ifcs::SMatrix{3, 3, T}
end 

# mimics the lo_fc3_triplet type
struct FC3Data{T} <: FCData{3,T}
    idxs::SVector{3, Int} # indices in unitcell
    lvs::SVector{3, SVector{3, T}}
    rv1::SVector{3, T}
    rv2::SVector{3, T}
    rv3::SVector{3, T}
    ifcs::SArray{Tuple{3,3,3}, T}
end

struct FC4Data{T} <: FCData{4,T}
    idxs::SVector{4, Int} # indices in unitcell
    lvs::SVector{4, SVector{3, T}}
    rv1::SVector{3, T}
    rv2::SVector{3, T}
    rv3::SVector{3, T}
    rv4::SVector{3, T}
    ifcs::SArray{Tuple{3,3,3,3}, T}
end

# mimics the lo_fc2_atom type
struct AtomFC2{T,N} <: AtomFC{2,T,N}
    pairs::SVector{N, FC2Data{T}}
end

# mimics the lo_fc3_atom type
struct AtomFC3{T,N} <: AtomFC{3,T,N}
    triplets::SVector{N, FC3Data{T}}
end

# mimics the lo_fc4_atom type
struct AtomFC4{T,N} <: AtomFC{4,T,N}
    quartets::SVector{N, FC4Data{T}}
end

n_neighbors(::AtomFC{O,T,N}) where {O,T,N} = N

struct IFCs{O, T}
    na::Int # number of atoms in the cell
    r_cut::T
    # all interactions for each atom in the unitcell, has length na_uc
    # by translational symmetry we can re-build all other IFCs
    atoms::AbstractVector{<:AtomFC{O,T}}
end

get_interactions(data::IFCs{2}, i::Int) = data.atoms[i].pairs 
get_interactions(data::IFCs{3}, i::Int) = data.atoms[i].triplets 
get_interactions(data::IFCs{4}, i::Int) = data.atoms[i].quartets 


Base.show(io::IO, ifc::IFCs{O,T}) where {O,T} =
    print(io, "Order $(O) IFCs, cutoff = $(round(ifc.r_cut, digits = 6)) Ang from $(ifc.na) aotm unit-cell")


#######################

struct DistanceTableAtom{T,N}
    central_atom::Int
    # Vector from current center atom to other particle
    vs::SVector{N,SVector{3,T}}
    # Vector from current cetner atom to other particle with pbc accounted for
    lvs::SVector{N,SVector{3,T}}
    # image flags
    ns::SVector{N, SVector{3, Int16}}
    inds::SVector{N, Int}
    dists::SVector{N, T}
end

n_neighbors(::DistanceTableAtom{<:Any,N}) where N = N

struct DistanceTable{T}
    atoms::AbstractVector{<:DistanceTableAtom{T}}
end


##############################

struct CrystalStructure{T}
    x_frac::AbstractVector{SVector{3,T}}
    x_cart::AbstractVector{SVector{3,T}}
    species::AbstractVector{Symbol}
    L::AbstractMatrix{T}
end

Base.length(cs::CrystalStructure) = length(cs.x_frac)

"""

Parameters:
-----------
- `poscar_path::String` : Path to file to parse
- `poscar_is_frac::Bool = true` : Whether or not the coordinates in the file are fractional.
- `FT::Type{FLOAT_TYPE} = Float64` : Which type to parse the coords and cell data as
"""
function CrystalStructure(
        poscar_path::String; 
        poscar_is_frac::Bool = true,
        FT::Type{FLOAT_TYPE} = Float64
    ) where {FLOAT_TYPE <: AbstractFloat}

    species, x_frac, cell = read_poscar_data(
        poscar_path;
        ssposcar_is_frac = poscar_is_frac,
        store_frac_coords = true, 
        FT = FLOAT_TYPE
    )

    ls = length(species); lc = length(x_frac)
    if length(species) != length(x_frac)
        error("Length of species vectors built from POSCAR ($(ls)) does not match number of atoms in POSCAR ($(lc)). Weird.")
    end

    x_cart = to_cart_coords.(Ref(cell), x_frac)

    return CrystalStructure{FLOAT_TYPE}(x_frac, x_cart, species, cell)

end
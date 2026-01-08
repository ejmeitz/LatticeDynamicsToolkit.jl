export 
    IFC2, IFC3, IFC4,
    AmorphousIFC2

abstract type IFCs end

# mimics the lo_fc2_pair type
struct FC2Data 
    idxs::SVector{2, Int} # indices in unitcell
    lvs::SVector{2, SVector{3, Float64}} # cartesian lattice vectors for the unit cell these atoms belong to
    r::SVector{3, Float64} # cartesian vector between atom 1 and 2
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
# Amorphous IFC Types #
#######################

"""
    AmorphousIFC2

Second-order interatomic force constants for amorphous/disordered systems.
Unlike `IFC2`, this type stores force constants for ALL atom pairs (no translational symmetry).

# Fields
- `na::Int`: Number of atoms
- `r_cut::Float64`: Cutoff radius (bohr)
- `neighbors::Vector{Vector{Int}}`: `neighbors[i]` contains indices `j > i` within cutoff
- `blocks::Vector{Vector{SMatrix{3,3,Float64,9}}}`: `blocks[i][k]` is Φ_ij for `j = neighbors[i][k]`

# Notes
- Only upper-triangle pairs (i < j) are stored explicitly
- Φ_ji = Φ_ij^T by Newton's 3rd law (enforced during fitting)
- Diagonal blocks Φ_ii are computed from ASR: Φ_ii = -Σ_{j≠i} Φ_ij
"""
struct AmorphousIFC2
    na::Int
    r_cut::Float64
    neighbors::Vector{Vector{Int}}
    blocks::Vector{Vector{SMatrix{3,3,Float64,9}}}
end

Base.show(io::IO, ifc::AmorphousIFC2) =
    print(io, "Amorphous 2nd Order IFCs, cutoff = $(round(ifc.r_cut*bohr_to_A, digits = 5)) Ang, $(ifc.na) atoms")

get_kwarg(::AmorphousIFC2) = :ifc2

####################################
# Dense Matrix Conversion Methods  #
####################################

"""
    Matrix(ifc2::IFC2, sc::CrystalStructure) -> Matrix{Float64}

Convert crystalline IFC2 to a dense 3N×3N force constant matrix.
Assumes `ifc2` has been remapped to the supercell `sc`.

The resulting matrix enforces:
- Exact symmetry by construction: Φ_ij = Φ_ji^T
- Acoustic sum rule: Σ_j Φ_ij = 0
"""
function (::Type{<:Matrix})(ifc2::IFC2, sc::CrystalStructure)

    ifc2.na != length(sc) && error(ArgumentError("Cannot convert IFC to Matrix. IFC are built on $(ifc2.na) cell, but supercell has $(length(sc)). Try remapping IFCs to the supercell."))

    n_dof = 3*ifc2.na
    out = zeros(Float64, n_dof, n_dof)
    
    @inbounds for a1 in 1:ifc2.na
        r1 = 3*(a1-1)
        for pair in get_interactions(ifc2, a1)
            a2 = pair.idxs[2]
            r2 = 3*(a2-1)
            out[r1+1:r1+3, r2+1:r2+3] .= pair.ifcs 
        end
    end

    # enforce exact symmetry
    @inbounds for j in 1:n_dof, i in j+1:n_dof
        s = 0.5 * (out[i,j] + out[j,i])
        out[i,j] = s; out[j,i] = s
    end

    # enforce acoustic sum rule
    for i in 1:ifc2.na # index of block matrix
        for α in 1:3
            for β in 1:3
                ii = 3*(i-1) + α
                jj = 3*(i-1) + β # i == j because we're on diagonal
                out[ii,jj] = 0.0 # remove existing value
                out[ii,jj] = @views -1*sum(out[ii, β:3:end])
            end
        end
    end

    return out
end

"""
    Matrix(ifc2::AmorphousIFC2) -> Matrix{Float64}

Convert AmorphousIFC2 to a dense 3N×3N force constant matrix.

The matrix is constructed with:
- Symmetry Φ_ij = Φ_ji^T (by construction from fitting)
- ASR via row sums: Φ_ii[α,β] = -Σ_j Φ_ij[α,β]
- Exact symmetry enforced manually via symmetrization
"""
function (::Type{<:Matrix})(ifc2::AmorphousIFC2)
    na = ifc2.na
    n_dof = 3 * na
    out = zeros(Float64, n_dof, n_dof)
    
    # Fill off-diagonal blocks (upper triangle stored, lower from symmetry)
    # Diagonal blocks remain zero for now
    @inbounds for i in 1:na
        ri = 3*(i-1)
        for (k, j) in enumerate(ifc2.neighbors[i])
            rj = 3*(j-1)
            Φ_ij = ifc2.blocks[i][k]
            # Upper triangle: i < j
            out[ri+1:ri+3, rj+1:rj+3] .= Φ_ij
            # Lower triangle: Φ_ji = Φ_ij^T
            out[rj+1:rj+3, ri+1:ri+3] .= transpose(Φ_ij)
        end
    end
    
    # Apply ASR via row sums: diagonal block (i,i)[α,β] = -sum of row α, every 3rd column starting at β
    @inbounds for i in 1:na
        for α in 1:3
            for β in 1:3
                ii = 3*(i-1) + α
                jj = 3*(i-1) + β
                out[ii, jj] = -sum(@view out[ii, β:3:n_dof])
            end
        end
    end

    #Enforce exact symmetry
    @inbounds for j in 1:n_dof, i in j+1:n_dof
        s = 0.5 * (out[i,j] + out[j,i])
        out[i,j] = s; out[j,i] = s
    end
    
    return out
end

export remap

function to_frac_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L
    return mod.(cell \ position, L(1.0))
end

# Assumes coords are all in box already
function remap(
        r_frac_ss::AbstractVector{SVector{3,T}},
        atom_types_ss, 
        A::AbstractMatrix{T},
        ifcs::Vararg{<:IFCs, N}
    ) where {T,N}

    # Check input
    n_ss = length(r_frac_ss)
    if n_ss != length(atom_types_ss)
        throw(DimensionMismatch("the number of atomic positions in the supercell $(n_ss) do not match the number of atom types ($(length(atom_types_ss)))"))
    end

    n_uc_ifc = getproperty.(ifcs, :na)
    r_cut_ifc = getproperty.(ifcs, :r_cut)

    if !allequal(n_uc_ifc)
        throw(ArgumentError("You passed IFCs built from different size cells $(n_uc_ifc). You must call remap on each separately."))
    end
    n_uc = n_uc_ifc[1]

    # Build mapping from supercell atom --> unit cell atom (Spglib.jl)
    lattice = Cell(A, r_frac_ss, atom_types_ss)
    dset = get_Dataset(cell, 1e-5)
    ss_to_uc_map = #TODO

    if lattice.n_std_atoms != n_uc
        throw(ArgumentError("Spglib identified your supercell as having a $(lattice.n_std_atoms) atom unit cell. But you gave me $(n_uc) atoms in the unitcell."))
    end

    #! CHECK THAT std_positions match n_uc


    # Build neighborlist for super cell (CellListMap.jl)
    r_cart_ss = to_frac_coords.(Ref(A), r_frac_ss)

    # Only build neighbor lists for unique cutoff distances
    unique_rcs = unique(r_cut_ifc)
    nl_map = [findfirst(rc -> rc == ifc.r_cut, unique_rcs) for ifc in ifcs]
    #! CHECK IT WANTS UNITCELL AS ROWS
    nls = [neighborlist(r_cart_ss, rc + T(1e-4); unitcell = A', showprogress = true) for rc in unique_rcs]

    return [remap(ifcs[i], nls[nl_map[i]], ss_to_uc_map) for i in eachindex(ifcs)]
end

function remap(ifc2::IFCs{2,T}, nl, ss_to_uc_map) where T

end

function remap(ifc3::IFCs{3,T}, nl, ss_to_uc_map) where T

end

function remap(ifc4::IFCs{4,T}, nl, ss_to_uc_map) where T

end
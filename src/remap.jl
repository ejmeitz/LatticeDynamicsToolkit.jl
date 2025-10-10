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


"""
    map_super_to_unitcell(
        L_uc, L_sc,
        x_frac_uc::AbstractVector{<:SVector{3,T}},
        x_frac_sc::AbstractVector{<:SVector{3,T}};
        species_uc=nothing, species_sc=nothing,
        tol::Real=1e-3,
    ) -> (s2u::Vector{Int}, translations::Vector{SVector{3,Int}})

Map each supercell atom → its unit-cell atom (index) and unit-cell translation (i,j,k).
Distances are computed in Å using the unit-cell lattice `L_uc`.
- `L_uc`, `L_sc`: 3×3 with **columns** a,b,c (Cartesian, Å)
- `x_frac_uc`, `x_frac_sc`: vectors of SVector{3} fractional coords
- `tol` is in **Å** (e.g., 1e-3–1e-2 Å for slightly relaxed structures)
"""
function map_super_to_unitcell(
        L_uc::AbstractMatrix{T},
        L_sc::AbstractMatrix{T},
        x_frac_uc::AbstractVector{<:SVector{3,T}},
        x_frac_sc::AbstractVector{<:SVector{3,T}};
        species_uc::Union{Nothing,AbstractVector}=nothing,
        species_sc::Union{Nothing,AbstractVector}=nothing,
        tol::Real=1e-3,
    ) where {T <: AbstractFloat}

    @assert size(L_uc) == (3,3) && size(L_sc) == (3,3)
    n_uc = length(x_frac_uc)
    n_sc = length(x_frac_sc)
    if (species_uc !== nothing) || (species_sc !== nothing)
        @assert species_uc !== nothing && species_sc !== nothing "Provide both species arrays."
        @assert length(species_uc) == n_uc && length(species_sc) == n_sc
    end

    # PBC wrap to [-0.5, 0.5) per component in fractional space
    pbc_wrap(d::SVector{3,T}) = d .- round.(d)

    s2u = zeros(Int, n_sc)
    # translations = Vector{SVector{3,Int}}(undef, n_sc)

    f_uc_raw = MVector{3, T}(0,0,0)
    f_red = MVector{3, T}(0,0,0)
    d_frac =  MVector{3, T}(0,0,0)
    d_cart = MVector{3, T}(0,0,0)

    M = inv(L_uc) * L_sc
    tol_sq = tol*tol

    @inbounds for j in 1:n_sc
        # Cartesian → uc frac (not reduced)
        mul!(f_uc_raw, M, x_frac_sc[j])
        # reduce to [0,1)
        f_red .= f_uc_raw .- floor.(f_uc_raw)

        #Δ     = f_uc_raw - f_red
        # t_ijk = SVector{3,Int}(round.(Int, Δ))  # translation in uc basis

        # nearest match in Å using uc lattice: distance = ‖L_uc * wrap(frac_diff)‖2
        best_i, best_d_sq = 0, T(Inf)
        for i in 1:n_uc
            if species_uc !== nothing
                species_uc[i] == species_sc[j] || continue
            end
            d_frac .= pbc_wrap(x_frac_uc[i] - f_red)
            mul!(d_cart, L_uc, d_frac)
            d_sq = dot(d_cart, d_cart)
            if d_sq < best_d_sq
                best_d_sq = d_sq
                best_i = i
            end
        end

        if !(best_i != 0 && best_d_sq ≤ tol_sq)
            throw(error("No unit-cell match within tol for supercell atom $j (best_d=$(sqrt(best_d)) Å)"))
        end
        
        s2u[j] = best_i
        # translations[j] = t_ijk
    end

    return s2u#, translations
end


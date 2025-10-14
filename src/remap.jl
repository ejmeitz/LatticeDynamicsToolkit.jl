export remap


# Assumes coords are all in box already
function remap(
        x_frac_sc::AbstractVector{SVector{3,T}},
        atom_types_sc, 
        L_sc::AbstractMatrix{T},
        x_frac_uc::AbstractVector{SVector{3,T}},
        atom_types_uc,
        L_uc::AbstractMatrix{T},
        ifcs::Vararg{<:IFCs, N}
    ) where {T,N}

    # Check input
    n_uc = length(x_frac_uc)
    n_sc = length(x_frac_sc)
    if n_sc != length(atom_types_sc)
        throw(DimensionMismatch("the number of atomic positions in the super cell $(n_sc) do not match the number of atom types ($(length(atom_types_sc)))"))
    end
    if n_uc != length(atom_types_uc)
        throw(DimensionMismatch("the number of atomic positions in the unitcell $(n_uc) do not match the number of atom types ($(length(atom_types_uc)))"))
    end

    n_uc_ifc = getproperty.(ifcs, :na)
    r_cut_ifc = getproperty.(ifcs, :r_cut)

    if !allequal(n_uc_ifc)
        throw(ArgumentError("You passed IFCs built from different size unitcells $(n_uc_ifc). You must call remap on each separately."))
    end
    if !all(n_uc .== n_uc_ifc)
        throw(ArgumentError("Your force constants are calcualted on a unitcells with $(n_uc_ifc) atoms, but you passed a unitcell with $(n_uc) atoms."))
    end

    ss_to_uc_map = map_super_to_unitcell(
        L_uc,
        L_sc,
        x_frac_uc,
        x_frac_sc;
        species_uc=atom_types_uc,
        species_sc=atom_types_sc
    ) 

    #######################
    # Build Neighborlists #
    #######################

    r_cart_ss = to_cart_coords.(Ref(L_sc), x_frac_sc)

    # Only build neighbor lists for unique cutoff distances
    unique_rcs = unique(r_cut_ifc)
    nl_map = [findfirst(rc -> rc == ifc.r_cut, unique_rcs) for ifc in ifcs]
    # TDEP adds a little tol to r_cut
    nls = [neighborlist(r_cart_ss, rc + T(1e-4); unitcell = L_sc, showprogress = true) for rc in unique_rcs]

    return [remap(ifcs[i], nls[nl_map[i]], ss_to_uc_map) for i in eachindex(ifcs)]
end

function remap(ifc2::IFCs{2,T}, nl, ss_to_uc_map) where T

end

function remap(ifc3::IFCs{3,T}, nl, ss_to_uc_map) where T

end

function remap(ifc4::IFCs{4,T}, nl, ss_to_uc_map) where T

end



@inline sqnorm(v::SVector{3,<:Real}) = dot(v, v)

"""
Quantize a triple of displacement vectors into a tolerant Dict key.
Use a quantum q ~ sqrt(lo_sqtol) from your code or whatever tolerance you prefer.
"""
@inline function key9(rv2::SVector{3,T}, rv3::SVector{3,T}, rv4::SVector{3,T}, q::T) where {T<:AbstractFloat}
    return (
        round(Int, rv2[1]/q), round(Int, rv2[2]/q), round(Int, rv2[3]/q),
        round(Int, rv3[1]/q), round(Int, rv3[2]/q), round(Int, rv3[3]/q),
        round(Int, rv4[1]/q), round(Int, rv4[2]/q), round(Int, rv4[3]/q),
    )
end

"""
Build a fast lookup for a unit-cell atom uca: (rv2,rv3,rv4) -> (m, mwm, idx)
Assumes fc.atom[uca].quartet[ii] has fields rv2, rv3, rv4, m, mwm.
"""
function build_uc_quartet_lookup(fc, uca; q::Float64=1e-8)
    lut = Dict{NTuple{9,Int}, Tuple{Float64,Float64,Int}}()
    qs = q
    for ii in 1:fc.atom[uca].n
        qrt = fc.atom[uca].quartet[ii]
        k = key9(qrt.rv2, qrt.rv3, qrt.rv4, qs)
        lut[k] = (qrt.m, ii)
    end
    return lut
end

# ---- neighbor table interface you said you have ---------------------------
# Expected minimal API of `dt` (distance table / neighbor list):
#   dt.particle[a1].n                :: Int
#   dt.particle[a1].ind[i]           :: Int              (neighbor index)
#   dt.particle[a1].v[i]             :: SVector{3,T}     (Cartesian displacement)
#   dt.particle[a1].lv[i]            :: SVector{3,Int}   (image lattice vector)
#
# And a constructor like:
#   build_dist_table!(dt, ss_positions::Vector{SVector{3,T}}, L_sc::SMatrix{3,3,T}, cutoff::T)
#
# If your names differ, just adapt in the spots below.

# ---- main remap -----------------------------------------------------------

"""
    remap_fc4!(fcss, fc, uc, ss, dt, s2u; tol_key=1e-8, pad=10.0*sqrt_tol)

- `fc`  : unit-cell force constants (4th order)
- `uc`  : unit-cell structure (only used for fc scan; keep for parity with Fortran)
- `ss`  : supercell structure (must expose `na` and `positions` + `L_sc` you use to build dt)
- `dt`  : neighbor table built over the supercell with cutoff ≥ rc + pad
- `s2u` : mapping supercell atom index → unit-cell atom index (your existing map)
- `tol_key` : quantization (in Å) used for tolerant matching of (rv2,rv3,rv4)
- `pad` : extra margin added to `rc` when building the neighbor list

After the call, `fcss` mirrors your Fortran: one atom entry per supercell atom,
same number of quartets as its corresponding unit-cell atom, and (m,mwm) copied.
"""
function remap_fc4(fc::IFCs{4,T}, dt::DistanceTable{T}, s2u::AbstractVector{<:Integer};
                    tol_key::Float64 = 1e-8) where T

    # 1) Actual cutoff rc from the unit-cell FC (max norm of rv2/rv3/rv4)
    rc = 0.0
    for uca in 1:fc.na
        for qidx in 1:fc.atom[uca].n
            qrt = fc.atom[uca].quartet[qidx]
            rc = max(rc, sqrt(sqnorm(qrt.rv2)))
            rc = max(rc, sqrt(sqnorm(qrt.rv3)))
            rc = max(rc, sqrt(sqnorm(qrt.rv4)))
        end
    end
    rc_sq = rc*rc + 10*tol_key  # mirrors your rc_sq = rc**2 + 10*lo_sqtol


    # Pre-build LUTs per unique unit-cell atom species/index to avoid O(n^4) lookups
    # Many supercell atoms map to the same unit-cell index; cache those LUTs.
    lut_cache = Dict{Int, Dict{NTuple{9,Int},Tuple{Float64,Float64,Int}}}()

    data = AtomFC4{T}[]

    # 4) For every supercell atom, generate quartets by neighbor triples
    for a1 in 1:na_ss
        uca = s2u[a1]  
        n_uc_quartets = fc.atom[uca].n

        # LUT for this unit-cell prototype (rv2,rv3,rv4 -> (m,mwm,idx))
        lut = get!(lut_cache, uca) do
            build_uc_quartet_lookup(fc, uca; q=tol_key)
        end

        l = 0
        quartets = Vector{FC4Data{T}}(undef, n_uc_quartets)

        @inbounds for i in 1:n_neighbors
            vi = dt.atoms[a1].vs[i]
            lvi = dt.atoms[a1].lvs[i]
            ii  = dt.atoms[a1].inds[i]
            for j in 1:n_neighbors
                vj  = dt.atoms[a1].vs[j]
                if sqnorm(vi - vj) >= rc_sq; continue; end
                lvj = dt.atoms[a1].lvs[j]
                jj  = dt.atoms[a1].inds[j]
                for k in 1:n_neighbors
                    vk = dt.atoms[a1].vs[k]
                    if sqnorm(vj - vk) >= rc_sq; continue; end
                    if sqnorm(vi - vk) >= rc_sq; continue; end

                    l += 1

                    # Find matching UC quartet (by tolerant key of rv2/rv3/rv4)
                    key = key9(vi, vj, vk, tol_key)
                    ifcs, _ = get(lut, key, nothing)
                    if ifcs === nothing
                        error("Could not locate quartet for a1=$a1; this should be impossible.")
                    end

                    quartets[l] = FC4Data{T}(
                        SVector(a1, ii, jj, dt.atoms[a1].inds[k]),
                        SVector(SVector{3,Int}(0,0,0), lvi, lvj, dt.atoms[a1].lvs[k]),
                        SVector{3,T}(0,0,0),
                        vi,
                        vj,
                        vk,
                        ifcs
                    )

                end
            end
        end

        # Count sanity check (like your Fortran)
        if l != n_uc_quartets
            error("Inconsistent number of quartets for atom $a1 (got $l, expected $n_uc_quartets).")
        end

        push!(data, AtomFC4{T, n_uc_quartets}(SVector{n_uc_quartets}(quartets)))

    end

    return IFCs{4, T}(na, sqrt(rc_sq), data)
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


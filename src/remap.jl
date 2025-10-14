export remap


# Assumes coords are all in box already
function remap(
        new_sc::CrystalStructure{T},
        uc::CrystalStructure{T},
        ifcs_old::Vararg{<:IFCs, N}
    ) where {T <: AbstractFloat, N}

    # Check input
    n_uc = length(uc)
    n_uc_ifc = getproperty.(ifcs_old, :na)
    r_cut_ifc = getproperty.(ifcs_old, :r_cut)

    if !allequal(n_uc_ifc)
        throw(ArgumentError("You passed ifcs_old built from different size unitcells $(n_uc_ifc). You must call remap on each separately."))
    end
    if !all(n_uc .== n_uc_ifc)
        throw(ArgumentError("Your force constants are calcualted on a unitcells with $(n_uc_ifc) atoms, but you passed a unitcell with $(n_uc) atoms."))
    end

    ss_to_uc_map = map_super_to_unitcell(uc, new_sc) 

    #######################
    # Build Neighborlists #
    #######################

    # Only build neighbor lists for unique cutoff distances
    unique_rcs = unique(r_cut_ifc)
    nl_map = [findfirst(rc -> rc == ifc.r_cut, unique_rcs) for ifc in ifcs_old]
    # TDEP adds a little tol to r_cut
    dts = [make_distance_table(new_sc, rc + T(1e-4)) for rc in unique_rcs]

    return [remap(ifcs_old[i], dts[nl_map[i]], ss_to_uc_map) for i in eachindex(ifcs_old)]
end

function remap(ifc2::IFCs{2,T}, nl, ss_to_uc_map) where T

end

function remap(ifc3::IFCs{3,T}, nl, ss_to_uc_map) where T

end

function remap(ifc4::IFCs{4,T}, nl, ss_to_uc_map) where T

end



@inline sqnorm(v::SVector{3,<:Real}) = dot(v, v)


# If your UC quartets only store lv (Å) and not integer flags, convert once:
@inline function lv_cart_to_n(L::AbstractMatrix{T}, lv::SVector{3,T}) where {T<:AbstractFloat}
    n_real = -(L \ lv)
    return SVector{3,Int}(round.(n_real; r=RoundNearestTiesAway))
end

"""
Build UC lookup for atom `uca`, keyed by **(i2,n2,i3,n3,i4,n4)**, all integers.
Value is `(m, idx)`
"""
function build_uc_quartet_lookup_int(fc4::IFCs{4, T}, uca, L_uc::AbstractMatrix{T}) where T
    KEY_TYPE = Tuple{Int, Int, Int, SVector{3, T}, SVector{3, T}, SVector{3,T}}
    VALUE_TYPE = SArray{Tuple{3,3,3,3}, T}
    lut = Dict{KEY_TYPE, VALUE_TYPE}()

    for q in get_interactions(fc4, uca)
        q = fc.atom[uca].quartet[ii]

        n2 = lv_cart_to_n(L_uc, q.lv2)
        n3 = lv_cart_to_n(L_uc, q.lv3)
        n4 = lv_cart_to_n(L_uc, q.lv4)

        k = (q.idxs..., n2, n3, n4)
        lut[k] = q.ifcs 
    end
    return lut
end

# ---------- main remap (uses integer flags from dt) ------------------------

"""
    remap_fc4(fc::IFCs{4,T}, dt::DistanceTable{T}, s2u::AbstractVector{<:Integer};
              L_uc::Union{Nothing,SMatrix}=nothing)

- Uses the stored integer image flags `dt.atoms[a].ns[i]`.
- Keys the UC lookup with integers (no float tolerance).
- Expects your UC quartets either already have integer flags (n2,n3,n4), or pass `L_uc`
  so we can convert `lv`→integer once in the lookup builder.

Returns your IFCs with per-atom quartets filled; only `m` is copied (no mwm).
"""
function remap_fc4(fc::IFCs{4,T}, dt::DistanceTable{T}, s2u::AbstractVector{<:Integer},
                   L_uc::AbstractMatrix{T}) where {T}

    rc = zero(Float64)
    for uca in 1:fc.na, q in fc.atom[uca].quartet
        rc = max(rc, sqrt(sqnorm(q.rv2)))
        rc = max(rc, sqrt(sqnorm(q.rv3)))
        rc = max(rc, sqrt(sqnorm(q.rv4)))
    end
    rc_sq = rc^2

    KEY_TYPE = Tuple{Int, Int, Int, SVector{3, T}, SVector{3, T}, SVector{3,T}}
    VALUE_TYPE = SArray{Tuple{3,3,3,3}, T}
    all_luts = Dict{Int, Dict{KEY_TYPE, VALUE_TYPE}}()

    data  = AtomFC4{T}[]
    na_ss = length(dt.atoms)

    # 3) per supercell atom
    for a1 in 1:na_ss
        uca = s2u[a1]
        n_uc_quartets = fc.atom[uca].n

        lut = get!(all_luts, uca) do
            build_uc_quartet_lookup_int(fc, uca, L_uc)
        end

        l = 0
        quartets = Vector{FC4Data{T}}(undef, n_uc_quartets)

        inds = dt.atoms[a1].inds
        vs   = dt.atoms[a1].vs     # Cartesian min-image vectors (Å)
        ns   = dt.atoms[a1].ns     # integer image flags (same basis as UC LUT)
        nnb  = length(inds)

        @inbounds for i in 1:nnb
            vi  = vs[i]
            ii  = inds[i]
            n2  = ns[i]

            for j in 1:nnb
                vj  = vs[j]
                (sqnorm(vi - vj) <= rc_sq) || continue
                jj  = inds[j]
                n3  = ns[j]

                for k in 1:nnb
                    vk  = vs[k]
                    (sqnorm(vj - vk) <= rc_sq) || continue
                    (sqnorm(vi - vk) <= rc_sq) || continue
                    kk  = inds[k]
                    n4  = ns[k]

                    # integer key: (unit-cell neighbor index, integer image flags)
                    key = (s2u[ii], s2u[jj], s2u[kk], n2, n3, n4)

                    m = get(lut, key, nothing)
                    m === nothing && error("No UC quartet for a1=$a1 (key=$(key))")

                    l += 1

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

        (l == n_uc_quartets) || error("Inconsistent number of quartets for atom $a1: got $l, expected $n_uc_quartets")
        push!(data, AtomFC4{T, n_uc_quartets}(SVector{n_uc_quartets}(quartets)))
    end

    return IFCs{4, T}(na_ss, sqrt(rc_sq), data)
end

"""
    map_super_to_unitcell(
        uc::CrystalStructure{T},
        sc::CrystalStructure{T},
        tol::Real=1e-3,
    ) -> (s2u::Vector{Int}, translations::Vector{SVector{3,Int}})

Map each supercell atom → its unit-cell atom (index) and unit-cell translation (i,j,k).
Distances are computed in Å using the unit-cell lattice `L_uc`.
- `tol` is in **Å** (e.g., 1e-3–1e-2 Å for slightly relaxed structures)
"""
function map_super_to_unitcell(
        uc::CrystalStructure{T},
        sc::CrystalStructure{T},
        tol::Real=1e-3,
    ) where {T <: AbstractFloat}

    n_uc = length(uc)
    n_sc = length(sc)

    # PBC wrap to [-0.5, 0.5) per component in fractional space
    pbc_wrap(d::SVector{3,T}) = d .- round.(d)

    s2u = zeros(Int, n_sc)
    # translations = Vector{SVector{3,Int}}(undef, n_sc)

    f_uc_raw = MVector{3, T}(0,0,0)
    f_red = MVector{3, T}(0,0,0)
    d_frac =  MVector{3, T}(0,0,0)
    d_cart = MVector{3, T}(0,0,0)

    M = inv(uc.L) * sc.L
    tol_sq = tol*tol

    @inbounds for j in 1:n_sc
        # Cartesian → uc frac (not reduced)
        mul!(f_uc_raw, M, sc.x_frac[j])
        # reduce to [0,1)
        f_red .= f_uc_raw .- floor.(f_uc_raw)

        #Δ     = f_uc_raw - f_red
        # t_ijk = SVector{3,Int}(round.(Int, Δ))  # translation in uc basis

        # nearest match in Å using uc lattice: distance = ‖L_uc * wrap(frac_diff)‖2
        best_i, best_d_sq = 0, T(Inf)
        for i in 1:n_uc
            uc.species[i] == sc.species[j] || continue
            
            d_frac .= pbc_wrap(uc.x_frac[i] - f_red)
            mul!(d_cart, uc.L, d_frac)
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


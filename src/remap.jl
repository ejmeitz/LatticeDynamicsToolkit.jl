export remap


function remap(new_sc::CrystalStructure, uc::CrystalStructure, ifcs::IFCs...)

    remap_checks(uc, ifcs...)

    ss_to_uc_map = map_super_to_unitcell(uc, new_sc) 

    return remap.(ifcs, Ref(new_sc), Ref(ss_to_uc_map))

end

function remap_checks(
        uc::CrystalStructure,
        ifcs::IFCs...
    )

    # Check input
    n_uc = length(uc)
    n_uc_ifc = [ifc.na for ifc in ifcs]

    if !allequal(n_uc_ifc)
        throw(ArgumentError("You passed ifcs_old built from different size unitcells $(n_uc_ifc)."))
    end
    if !all(n_uc .== n_uc_ifc)
        throw(ArgumentError("Your force constants are calcualted on a cell with $(n_uc_ifc) atoms, but you passed a cell with $(n_uc) atoms."))
    end
end

function remap(
    fc::IFC2,
    sc::CrystalStructure,
    s2u::AbstractVector{Int};
) 

    na_sc = length(sc)
    data  = Vector{Vector{FC2Data}}(undef, na_sc)

    z3 = SVector{3,Float64}(0,0,0)

    @inbounds for a1 in 1:na_sc

        uca = s2u[a1]
        r0 = sc.x_cart[a1]
        
        pairs = get_interactions(fc, uca)
        pair_data = Vector{FC2Data}(undef, length(pairs))

        for (i,p) in enumerate(pairs)
            # Target absolute position in SC: r1 = r0 + r
            r1_cart = p.r + r0
            r1_frac = clean_fractional_coordinates.(sc.L_inv * r1_cart)

            # Find matching atom a2 in SC by fractional coordinates
            i2 = -1
            for a2 in 1:na_sc
                if sum(abs.(sc.x_frac[a2] - r1_frac)) < lo_tol
                    i2 = a2
                    break
                end
            end
            
            if i2 == 0
                error("Failed mapping second-order IFCs for atom $a1 (pair $i): no match for target position.")
            end

            # Fix the lattice image for atom 2 so that:
            r2_cart   = p.r - sc.x_cart[i2] + sc.x_cart[a1]
            r2_lat    = round.(Int, sc.L_inv * r2_cart) 
            lv2_cart  = sc.L * r2_lat 
            
            n2 = SVector{3,Int16}(r2_lat)

            # Store as FC2Data: (i1,i2), (lv1=0, lv2), (rv1=0, rv2=r), and the 3×3 ifcs
            pair_data[i] = FC2Data(
                SVector(a1, i2),
                SVector(z3, lv2_cart),
                p.r,
                n2,
                p.ifcs
            )
        end

        data[a1] = pair_data
    end

    return IFC2(na_sc, fc.r_cut, data)   # or IFCs{2,T}(na_sc, fc.r_cut, data) if you store the cutoff in fc
end


function remap(
        fc::IFC3,
        sc::CrystalStructure,
        s2u::AbstractVector{Int};
        n_threads::Integer = Threads.nthreads()
    )

    rc = 0.0
    for uca in 1:fc.na, t in get_interactions(fc, uca)
        rc = max(rc, sqrt(sqnorm(t.rv2)))
        rc = max(rc, sqrt(sqnorm(t.rv3)))
    end
    rc += lo_tol
    rc_sq  = rc * rc

    dt = make_distance_table(sc, rc; include_self = true)

    na_sc = length(dt.atoms)

    # One entry per supercell atom
    data = Vector{Vector{FC3Data}}(undef, na_sc)

    @inbounds @tasks for a1 in 1:na_sc
        @set ntasks=n_threads

        uca = s2u[a1]
        l = 0

        # Expect same # of triplets as its corresponding UC atom
        triplets = get_interactions(fc, uca)
        n_uc_triplets = length(triplets)
        triplet_data = Vector{FC3Data}(undef, n_uc_triplets)

        inds = dt.atoms[a1].inds
        vs   = dt.atoms[a1].vs
        lvs  = dt.atoms[a1].lvs
        nnb  = length(inds)

        for i in 1:nnb
            vi = vs[i]
            for j in 1:nnb
                vj = vs[j]
                # Only consider neighbor pairs within cutoff
                (sqnorm(vi - vj) < rc_sq) || continue
                l += 1

                # Find matching UC triplet by relative vectors (within tol)
                ifc_idx = -1
                for (idx, t) in enumerate(triplets)
                    (sqnorm(t.rv2 - vi) < lo_sqtol) || continue
                    (sqnorm(t.rv3 - vj) < lo_sqtol) || continue
                    ifc_idx = idx
                    break
                end

                if ifc_idx == -1
                    error("Could not locate triplet; cells/mapping likely inconsistent.")
                end

                triplet_data[l] = FC3Data(
                    SVector(a1, inds[i], inds[j]),
                    SVector(SVector{3,Float64}(0,0,0), lvs[i], lvs[j]),
                    SVector{3,Float64}(0,0,0),
                    vi,
                    vj,
                    dt.atoms[a1].ns[i],
                    dt.atoms[a1].ns[j],
                    fc.atoms[uca][ifc_idx].ifcs
                )
            end
        end

        (l == n_uc_triplets) || error("Inconsistent number of triplets for atom $a1: got $l, expected $n_uc_triplets")
        data[a1] = triplet_data
    end

    return IFC3(na_sc, fc.r_cut, data)
end


function remap(
        fc::IFC4,
        sc::CrystalStructure,
        s2u::AbstractVector{Int};
        n_threads::Integer = Threads.nthreads()
    ) 

    rc = 0.0
    for uca in 1:fc.na, q in get_interactions(fc, uca)
        rc = max(rc, sqrt(sqnorm(q.rv2)))
        rc = max(rc, sqrt(sqnorm(q.rv3)))
        rc = max(rc, sqrt(sqnorm(q.rv4)))
    end

    rc_sq = rc*rc + 10*lo_sqtol

    dt = make_distance_table(sc, rc + 10*lo_tol; include_self = true)

    na_sc = length(dt.atoms)
    data  = Vector{Vector{FC4Data}}(undef, na_sc)

    @tasks for a1 in 1:na_sc
        @set ntasks=n_threads

        uca = s2u[a1]

        quartets = get_interactions(fc, uca)
        # The corresponding atom in the unitcell has this
        # many interactions, so we expect atom `a1` to have
        # the same number of interactions
        n_uc_quartets = length(quartets)

        l = 0
        quartet_data = Vector{FC4Data}(undef, n_uc_quartets)

        inds = dt.atoms[a1].inds
        vs   = dt.atoms[a1].vs
        lvs = dt.atoms[a1].lvs
        nnb  = length(inds)

        @inbounds for i in 1:nnb
            vi  = vs[i]
            for j in 1:nnb
                vj  = vs[j]
                (sqnorm(vi - vj) < rc_sq) || continue
                for k in 1:nnb
                    vk  = vs[k]
                    (sqnorm(vj - vk) < rc_sq) || continue
                    (sqnorm(vi - vk) < rc_sq) || continue
                    l += 1

                    ifc_idx = -1
                    for (idx, q) in enumerate(quartets)
                        (sqnorm(q.rv2 - vi) <= lo_sqtol) || continue
                        (sqnorm(q.rv3 - vj) <= lo_sqtol) || continue
                        (sqnorm(q.rv4 - vk) <= lo_sqtol) || continue
                        ifc_idx = idx
                        break
                    end

                    if ifc_idx == -1
                        error("Could not locate quartet, should be impossible")
                    end

                    quartet_data[l] = FC4Data( 
                        SVector(a1, inds[i], inds[j], inds[k]), 
                        SVector(SVector{3,Float64}(0,0,0), lvs[i], lvs[j], lvs[k]),
                        SVector{3,Float64}(0,0,0), 
                        vi, vj, vk,
                        dt.atoms[a1].ns[i],
                        dt.atoms[a1].ns[j],
                        dt.atoms[a1].ns[k],
                        fc.atoms[uca][ifc_idx].ifcs
                    )

                end
            end
        end

        (l == n_uc_quartets) || error("Inconsistent number of quartets for atom $a1: got $l, expected $n_uc_quartets")
        data[a1] = quartet_data
    end

    return IFC4(na_sc, sqrt(rc_sq), data)
end

"""
    map_super_to_unitcell(
        uc::CrystalStructure,
        sc::CrystalStructure,
        tol::Real=1e-3,
    ) -> (s2u::Vector{Int}, translations::Vector{SVector{3,Int}})

Map each supercell atom → its unit-cell atom (index) and unit-cell translation (i,j,k).
Distances are computed in Å using the unit-cell lattice `L_uc`.
- `tol` is in **Å** (e.g., 1e-3–1e-2 Å for slightly relaxed structures)
"""
function map_super_to_unitcell(
        uc::CrystalStructure,
        sc::CrystalStructure,
        tol::Real=1e-3,
    )

    n_uc = length(uc)
    n_sc = length(sc)

    # PBC wrap to [-0.5, 0.5) per component in fractional space
    pbc_wrap(d) = d .- round.(d)

    s2u = zeros(Int, n_sc)
    # translations = Vector{SVector{3,Int}}(undef, n_sc)

    f_uc_raw = MVector{3, Float64}(0,0,0)
    f_red = MVector{3, Float64}(0,0,0)
    d_frac =  MVector{3, Float64}(0,0,0)
    d_cart = MVector{3, Float64}(0,0,0)

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
        best_i, best_d_sq = 0, Inf
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


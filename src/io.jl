export read_ifc2, read_ifc3, read_ifc4


"""
    read_ifc2(ifc2_path, ucposcar_path)
    read_ifc2(path, r_frac_uc, L_uc) -> IFCs{2, T}

Read a 2nd-order outfile.forceconstant from TDEP. This follows the logic implemented in TDEP.
Does NOT handle polar force constants.

Inputs
- `path::AbstractString`: file path
- `r_frac_uc::AbstractVector{SVector{3,T}}`: 3×na fractional coords of atoms in the unitcell
- `L_uc::AbstractMatrix{T}`: 3×3 lattice (columns are lattice vectors), Cartesian = A * fractional

"""
function read_ifc2(ifc2_path::AbstractString, ucposcar_path::AbstractString)
    _, x_frac_uc, L_uc = read_poscar_data(ucposcar_path)
    L_uc *= A_to_bohr
    return read_ifc2(ifc2_path, x_frac_uc, L_uc)
end

function read_ifc2(
        path::AbstractString,
        r_frac_uc::AbstractVector{SVector{3,Float64}},
        L_uc::AbstractMatrix{Float64};
        chop_tol::Float64 = 1e-13
    )

    @assert size(L_uc) == (3,3) "Lattice A must be 3×3."

    max_rcut = 0.0
    lv1_frac = SVector{3, Float64}(0,0,0)

    data = Vector{FC2Data}[]

    open(path, "r") do io
        na = readline_skip_text!(io, Int) # number of atoms in unit cell
        resize!(data, na)

        # cross-check given positions
        @assert length(r_frac_uc) == na "r_frac_uc has length $(length(r_frac_uc)) but file says there are $(na) atoms in the unitcell."

        # recomputed based on max from whole file
        _cutoff_ignored = readline_skip_text!(io, Float64)

        for a1 in 1:na
            n_neighbors = readline_skip_text!(io, Int)
            pair_data = Vector{FC2Data}(undef, n_neighbors) 
            for i in 1:n_neighbors
                a2 = readline_skip_text!(io, Int) # unitcell index of neighbor
                lv2_frac = read_svec3!(io, Float64)
                ifcs = read_mat3_rows!(io, Float64; conv = forceconstant_2nd_eVA_to_HartreeBohr)

                # frac positions of the two atoms in their image cells
                v1_frac = r_frac_uc[a1] # + lv1_frac 
                v2_frac = lv2_frac + SVector{3,Float64}(r_frac_uc[a2])

                # Convert to Cartesian
                lv1_cart = SVector{3,Float64}(L_uc * lv1_frac)
                lv2_cart = SVector{3,Float64}(L_uc * lv2_frac)
                r_cart   = SVector{3,Float64}(L_uc * (v2_frac - v1_frac))

                lv1_cart = chop3(lv1_cart, chop_tol)
                lv2_cart = chop3(lv2_cart, chop_tol)
                r_cart   = chop3(r_cart, chop_tol)

                max_rcut = max(max_rcut, norm(r_cart))

                # Calculate image flags
                n2 = SVector{3,Int16}(round.(Int16, v2_frac - v1_frac))

                pair_data[i] = FC2Data(
                    SVector{2,Int}(a1, a2),
                    SVector{2,SVector{3,Float64}}(lv1_cart, lv2_cart),
                    r_cart,
                    n2,
                    ifcs
                )
            end
            
            data[a1] = pair_data
        end

        # technically theres more polar stuff, but ignore that for now

        IFC2(na, max_rcut + lo_sqtol, data)
    end
end

"""
    read_ifc3(ifc3_path, ucposcar_path) -> IFCs{3,T}
    read_ifc3(path, r_frac_uc, L_uc) -> IFCs{3,T}

Read a 3rd-order `outfile.forceconstant_thirdorder` (TDEP-style).

Inputs
- `path::AbstractString`: file path
- `r_frac_uc::AbstractVector{SVector{3,T}}`: fractional positions of atoms in the unit cell (length = na)
- `L_uc::AbstractMatrix{T}`: 3×3 lattice with **columns as lattice vectors** so `cart = A * frac`

Returns
- `IFCs{3,T}` with `na_uc`, `r_cut`, and `atoms::Vector{AtomFC3{T,N}}`
"""
function read_ifc3(ifc3_path::AbstractString, ucposcar_path::AbstractString)
    _, x_frac, L_uc = read_poscar_data(ucposcar_path)
    L_uc *= A_to_bohr
    return read_ifc3(ifc3_path, x_frac, L_uc)
end

function read_ifc3(
        path::AbstractString,
        r_frac_uc::AbstractVector{SVector{3,Float64}},
        L_uc::AbstractMatrix{Float64};
        chop_tol::Float64 = 1e-13
    )

    @assert size(L_uc) == (3,3) "Lattice must be 3×3."

    data = Vector{FC3Data}[]
    max_rcut = 0.0

    L = SMatrix{3,3}(L_uc)

    open(path, "r") do io
        na = readline_skip_text!(io, Int)
        @assert length(r_frac_uc) == na "r_frac_uc length $(length(r_frac_uc)) but file says na=$na."
        _cutoff_ignored = readline_skip_text!(io, Float64)  # recomputed below

        resize!(data, na)

        for a1 in 1:na
            n_tr = readline_skip_text!(io, Int)
            trip_data = Vector{FC3Data}(undef, n_tr)

            for i in 1:n_tr
                i1 = readline_skip_text!(io, Int)
                i2 = readline_skip_text!(io, Int)
                i3 = readline_skip_text!(io, Int)

                # Fortran stores i1 explicitly; in the routine it equals a1. Keep and sanity-check.
                @assert i1 == a1 "Triplet central index mismatch: got i1=$i1 for a1=$a1."

                lv1_frac = read_svec3!(io, Float64)
                lv2_frac = read_svec3!(io, Float64)
                lv3_frac = read_svec3!(io, Float64)

                ifcs = read_tensor3!(io, Float64; conv = forceconstant_3rd_eVA_to_HartreeBohr)

                # Fractional positions of the three atoms in their image cells
                v1_frac = lv1_frac + r_frac_uc[a1]
                v2_frac = lv2_frac + r_frac_uc[i2]
                v3_frac = lv3_frac + r_frac_uc[i3]

                # Convert to Cartesian
                lv1_cart = L * lv1_frac
                lv2_cart = L * lv2_frac
                lv3_cart = L * lv3_frac

                v1 = L * v1_frac
                v2 = L * v2_frac
                v3 = L * v3_frac

                # Relative vectors (Fortran: rv1=0, rv2=v2-v1, rv3=v3-v1)
                rv1 = SVector{3, Float64}(0, 0, 0)
                rv2 = chop3(v2 - v1, chop_tol)
                rv3 = chop3(v3 - v1, chop_tol)

                lv1_cart = chop3(lv1_cart, chop_tol)
                lv2_cart = chop3(lv2_cart, chop_tol)
                lv3_cart = chop3(lv3_cart, chop_tol)

                # Image flags
                n2 = SVector{3,Int16}(round.(Int16, v2_frac - v1_frac))
                n3 = SVector{3,Int16}(round.(Int16, v3_frac - v1_frac))

                # Cutoff per Fortran: max(||rv2||, ||rv3||, ||rv3 - rv2||)
                max_rcut = max(max_rcut, norm(rv2))
                max_rcut = max(max_rcut, norm(rv3))
                max_rcut = max(max_rcut, norm(rv3 - rv2))

                trip_data[i] = FC3Data(
                    SVector{3,Int}(i1, i2, i3),
                    SVector{3,SVector{3,Float64}}(lv1_cart, lv2_cart, lv3_cart),
                    rv1, rv2, rv3,
                    n2, n3,
                    ifcs
                )
            end

            data[a1] = trip_data
        end
        IFC3(na, max_rcut + lo_sqtol, data)
    end
end

"""
    read_ifc4(ifc4_path, ucposcar_path) -> IFCs{4,T}
    read_ifc4(path, r_frac_uc, A) -> IFCs{4,T}

Read a 4th-order TDEP-style force constant file.

Inputs
- `path::AbstractString`
- `r_frac_uc::AbstractVector{SVector{3,T}}`: fractional positions for the unit cell (length = na)
- `L_uc::AbstractMatrix{T}` (3×3): lattice with **columns** as lattice vectors ⇒ `cart = A * frac`

Returns
- `IFCs{4,T}` with `na_uc`, `r_cut`, and `atoms::Vector{AtomFC4{T,N}}`
"""
function read_ifc4(ifc4_path::AbstractString, ucposcar_path::AbstractString)
    _, x_frac, L_uc = read_poscar_data(ucposcar_path)
    L_uc *= A_to_bohr
    return read_ifc4(ifc4_path, x_frac, L_uc)
end

function read_ifc4(path::AbstractString,
                    r_frac_uc::AbstractVector{SVector{3,Float64}},
                    L_uc::AbstractMatrix{Float64};
                    chop_tol = 1e-13)

    @assert size(L_uc) == (3,3) "Lattice A must be 3×3."

    data = Vector{FC4Data}[]
    max_rc_sq = 0.0

    L = SMatrix{3,3}(L_uc)

    open(path, "r") do io
        na = readline_skip_text!(io, Int)
        resize!(data, na)
        @assert length(r_frac_uc) == na "r_frac_uc length $(length(r_frac_uc)) but file says na=$na."
        _cutoff_ignored = readline_skip_text!(io, Float64)  # recomputed below

        for a1 in 1:na
            n_quart = readline_skip_text!(io, Int)
            quartets = Vector{FC4Data}(undef, n_quart)

            @sync for i in 1:n_quart
                i1 = readline_skip_text!(io, Int)
                i2 = readline_skip_text!(io, Int)
                i3 = readline_skip_text!(io, Int)
                i4 = readline_skip_text!(io, Int)
                @assert i1 == a1 "Quartet central index mismatch: got i1=$i1 for a1=$a1."

                lv1_frac = read_svec3!(io, Float64)
                lv2_frac = read_svec3!(io, Float64)
                lv3_frac = read_svec3!(io, Float64)
                lv4_frac = read_svec3!(io, Float64)

                ifcs = read_tensor4!(io, Float64; conv = forceconstant_4th_eVA_to_HartreeBohr)

                # fractional positions including image shifts
                v1_frac = lv1_frac + r_frac_uc[a1]
                v2_frac = lv2_frac + r_frac_uc[i2]
                v3_frac = lv3_frac + r_frac_uc[i3]
                v4_frac = lv4_frac + r_frac_uc[i4]

                # Cartesian shifts and positions
                lv1_cart = L * lv1_frac
                lv2_cart = L * lv2_frac
                lv3_cart = L * lv3_frac
                lv4_cart = L * lv4_frac

                v1 = L * v1_frac
                v2 = L * v2_frac
                v3 = L * v3_frac
                v4 = L * v4_frac

                # relative vectors (Fortran: rv1=0)
                rv1 = SVector{3, Float64}(0,0,0)
                rv2 = chop3(v2 - v1, chop_tol)
                rv3 = chop3(v3 - v1, chop_tol)
                rv4 = chop3(v4 - v1, chop_tol)

                lv1_cart = chop3(lv1_cart, chop_tol)
                lv2_cart = chop3(lv2_cart, chop_tol)
                lv3_cart = chop3(lv3_cart, chop_tol)
                lv4_cart = chop3(lv4_cart, chop_tol)

                # Image flags
                n2 = round.(Int16, v2_frac .- v1_frac)
                n3 = round.(Int16, v3_frac .- v1_frac)
                n4 = round.(Int16, v4_frac .- v1_frac)

                # cutoff per Fortran
                max_rc_sq = max(max_rc_sq, sqnorm(rv2))
                max_rc_sq = max(max_rc_sq, sqnorm(rv3))
                max_rc_sq = max(max_rc_sq, sqnorm(rv4))
                max_rc_sq = max(max_rc_sq, sqnorm(rv2 - rv3))
                max_rc_sq = max(max_rc_sq, sqnorm(rv2 - rv4))
                max_rc_sq = max(max_rc_sq, sqnorm(rv3 - rv4))

                quartets[i] = FC4Data(
                    SVector{4,Int}(i1, i2, i3, i4),
                    SVector{4,SVector{3,Float64}}(lv1_cart, lv2_cart, lv3_cart, lv4_cart),
                    rv1, rv2, rv3, rv4,
                    n2, n3, n4,
                    ifcs
                )
            end

            data[a1] = quartets
        end
        IFC4(na, sqrt(max_rc_sq) + lo_sqtol, data)
    end
end

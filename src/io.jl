export read_fc2

# Helper functions
readline_skip_text!(io, T) = parse(T, first(split(strip(readline(io)))))

function read_vec3!(io, T)
    xs = split(strip(readline(io)))
    @assert length(xs) == 3 "Expected 3 components for a 3-vector."
    SVector{3,T}(parse.(T, xs))
end

"""
    read_ifc2(path, R_frac, A) -> IFCs{2, T}

Read a 2nd-order outfile.forceconstant from TDEP. This follows the logic implemented in TDEP.
Does NOT handle polar force constants.

Inputs
- `path::AbstractString`: file path
- `r_frac_uc::AbstractVector{SVector{3,T}}`: 3×na fractional coords of atoms in the unitcell
- `A::AbstractMatrix{T}`: 3×3 lattice (columns are lattice vectors), Cartesian = A * fractional

"""
function read_ifc2(path::AbstractString,
                        r_frac_uc::AbstractVector{SVector{3,T}},
                        A::AbstractMatrix{T};
                        chop_tol::T = T(1e-13)) where {T<:Real}

    
    function read_mat3_rows!(io) 
        r1 = read_vec3!(io, T)
        r2 = read_vec3!(io, T)
        r3 = read_vec3!(io, T)
        M = hcat(r1, r2, r3)'   # r1 is first row, etc.
        return SMatrix{3,3,T}(M)     
    end

    chop3(v::SVector{3,T}) = SVector{3,T}(ntuple(i->(abs(v[i]) < chop_tol ? zero(T) : v[i]), 3))

    @assert size(A) == (3,3) "Lattice A must be 3×3."

    max_rcut = zero(T)

    data = AtomFC2{T}[]

    open(path, "r") do io
        na = readline_skip_text!(io, Int) # number of atoms in unit cell

        # cross-check given positions
        @assert length(r_frac_uc) == na "r_frac_uc has length $(length(r_frac_uc)) but file says there are $(na) atoms in the unitcell."

        # recomputed based on max from whole file
        _cutoff_ignored = readline_skip_text!(io, Float64)

        for a1 in 1:na
            n_neighbors = readline_skip_text!(io, Int)
            pair_data = Vector{FC2Data{T}}(undef, n_neighbors) 
            for i in 1:n_neighbors
                a2 = readline_skip_text!(io, Int) # unitcell index of neighbor
                lv2_frac = read_vec3!(io, T)
                ifcs = read_mat3_rows!(io)

                # Fortran routine uses lv1 == 0
                lv1_frac = SVector{3,T}(0,0,0)

                # frac positions of the two atoms in their image cells
                v1_frac = lv1_frac + SVector{3,T}(r_frac_uc[a1])
                v2_frac = lv2_frac + SVector{3,T}(r_frac_uc[a2])

                # Convert to Cartesian
                lv1_cart = SVector{3,T}(A * lv1_frac)
                lv2_cart = SVector{3,T}(A * lv2_frac)
                r_cart   = SVector{3,T}(A * (v2_frac - v1_frac))

                lv1_cart = chop3(lv1_cart)
                lv2_cart = chop3(lv2_cart)
                r_cart   = chop3(r_cart)

                max_rcut = max(max_rcut, norm(r_cart))

                pair_data[i] = FC2Data{T}(
                    SVector{2,Int}(a1, a2),
                    SVector{2,SVector{3,T}}(lv1_cart, lv2_cart),
                    r_cart,
                    ifcs
                )
            end
            
            # Would be nice to not have to copy into SVector here
            push!(data, AtomFC2{T, n_neighbors}(SVector{n_neighbors}(pair_data)))
        end

        # technically theres more polar stuff, but ignore that for now

        r_cut = max_rcut + sqrt(eps(T))
        println(typeof(data))
        return IFCs{2, T}(na, r_cut, data)

    end
end

"""
    read_ifc3(path, r_frac_uc, A) -> IFCs{3,T}

Read a 3rd-order `outfile.forceconstant_thirdorder` (TDEP-style).

Inputs
- `path::AbstractString`: file path
- `r_frac_uc::AbstractVector{SVector{3,T}}`: fractional positions of atoms in the unit cell (length = na)
- `A::AbstractMatrix{T}`: 3×3 lattice with **columns as lattice vectors** so `cart = A * frac`

Returns
- `IFCs{3,T}` with `na_uc`, `r_cut`, and `atoms::Vector{AtomFC3{T,N}}`
"""
function read_ifc3(path::AbstractString,
                           r_frac_uc::AbstractVector{SVector{3,T}},
                           A::AbstractMatrix{T};
                           chop_tol::T = T(1e-13)) where {T<:Real}

    @assert size(A) == (3,3) "Lattice A must be 3×3."

    # local helpers (reuse your global readline_skip_text! / read_vec3!)
    chop3(v::SVector{3,T}) = SVector{3,T}(ntuple(i -> (abs(v[i]) < chop_tol ? zero(T) : v[i]), 3))

    # Read 3×3×3 tensor in the same order as Fortran:
    # do ii=1:3; do jj=1:3; read v(1:3); m(ii,jj,:) = v
    function read_tensor3!(io)::SArray{Tuple{3,3,3},T}
        M = Array{T,3}(undef, 3, 3, 3)
        for ii in 1:3, jj in 1:3
            v = read_vec3!(io, T)
            @inbounds M[ii, jj, 1] = v[1]
            @inbounds M[ii, jj, 2] = v[2]
            @inbounds M[ii, jj, 3] = v[3]
        end
        return SArray{Tuple{3,3,3},T}(M)
    end

    data = AtomFC3{T}[]
    max_rcut = zero(T)

    open(path, "r") do io
        na = readline_skip_text!(io, Int)
        @assert length(r_frac_uc) == na "r_frac_uc length $(length(r_frac_uc)) but file says na=$na."
        _cutoff_ignored = readline_skip_text!(io, T)  # recomputed below

        for a1 in 1:na
            n_tr = readline_skip_text!(io, Int)
            trip_data = Vector{FC3Data{T}}(undef, n_tr)

            for i in 1:n_tr
                i1 = readline_skip_text!(io, Int)
                i2 = readline_skip_text!(io, Int)
                i3 = readline_skip_text!(io, Int)

                # Fortran stores i1 explicitly; in the routine it equals a1. Keep and sanity-check.
                @assert i1 == a1 "Triplet central index mismatch: got i1=$i1 for a1=$a1."

                lv1_frac = read_vec3!(io, T)
                lv2_frac = read_vec3!(io, T)
                lv3_frac = read_vec3!(io, T)

                ifcs = read_tensor3!(io)  # SArray{Tuple{3,3,3},T}

                # Fractional positions of the three atoms in their image cells
                v1_frac = lv1_frac + r_frac_uc[a1]
                v2_frac = lv2_frac + r_frac_uc[i2]
                v3_frac = lv3_frac + r_frac_uc[i3]

                # Convert to Cartesian
                lv1_cart = SVector{3,T}(A * lv1_frac)
                lv2_cart = SVector{3,T}(A * lv2_frac)
                lv3_cart = SVector{3,T}(A * lv3_frac)

                v1 = SVector{3,T}(A * v1_frac)
                v2 = SVector{3,T}(A * v2_frac)
                v3 = SVector{3,T}(A * v3_frac)

                # Relative vectors (Fortran: rv1=0, rv2=v2-v1, rv3=v3-v1)
                rv1 = SVector{3,T}(0, 0, 0)
                rv2 = chop3(v2 - v1)
                rv3 = chop3(v3 - v1)

                lv1_cart = chop3(lv1_cart)
                lv2_cart = chop3(lv2_cart)
                lv3_cart = chop3(lv3_cart)

                # Cutoff per Fortran: max(||rv2||, ||rv3||, ||rv3 - rv2||)
                max_rcut = max(max_rcut, norm(rv2))
                max_rcut = max(max_rcut, norm(rv3))
                max_rcut = max(max_rcut, norm(rv3 - rv2))

                trip_data[i] = FC3Data{T}(
                    SVector{3,Int}(i1, i2, i3),
                    SVector{3,SVector{3,T}}(lv1_cart, lv2_cart, lv3_cart),
                    rv1, rv2, rv3,
                    ifcs
                )
            end

            push!(data, AtomFC3{T, n_tr}(SVector{n_tr}(trip_data)))
        end

        r_cut = max_rcut + sqrt(eps(T))
        return IFCs{3, T}(na, r_cut, data)
    end
end

"""
    read_ifc4(path, r_frac_uc, A) -> IFCs{4,T}

Read a 4th-order TDEP-style force constant file.

Inputs
- `path::AbstractString`
- `r_frac_uc::AbstractVector{SVector{3,T}}`: fractional positions for the unit cell (length = na)
- `A::AbstractMatrix{T}` (3×3): lattice with **columns** as lattice vectors ⇒ `cart = A * frac`

Returns
- `IFCs{4,T}` with `na_uc`, `r_cut`, and `atoms::Vector{AtomFC4{T,N}}`
"""
function read_ifc4(path::AbstractString,
                           r_frac_uc::AbstractVector{SVector{3,T}},
                           A::AbstractMatrix{T};
                           chop_tol::T = T(1e-13)) where {T<:Real}

    @assert size(A) == (3,3) "Lattice A must be 3×3."

    chop3(v::SVector{3,T}) = SVector{3,T}(ntuple(i -> (abs(v[i]) < chop_tol ? zero(T) : v[i]), 3))

    # Fortran loop order:
    # do ii=1:3; do jj=1:3; do kk=1:3; read v0(1:3); m(ii,jj,kk,:) = v0
    function read_tensor4!(io)::SArray{Tuple{3,3,3,3},T}
        M = Array{T,4}(undef, 3,3,3,3)
        for ii in 1:3, jj in 1:3, kk in 1:3
            v0 = read_vec3!(io, T)
            @inbounds M[ii,jj,kk,1] = v0[1]
            @inbounds M[ii,jj,kk,2] = v0[2]
            @inbounds M[ii,jj,kk,3] = v0[3]
        end
        return SArray{Tuple{3,3,3,3},T}(M)
    end

    data = AtomFC4{T}[]
    max_rcut = zero(T)

    open(path, "r") do io
        na = readline_skip_text!(io, Int)
        @assert length(r_frac_uc) == na "r_frac_uc length $(length(r_frac_uc)) but file says na=$na."
        _cutoff_ignored = readline_skip_text!(io, T)  # recomputed below

        for a1 in 1:na
            n_quart = readline_skip_text!(io, Int)
            quartets = Vector{FC4Data{T}}(undef, n_quart)

            for i in 1:n_quart
                i1 = readline_skip_text!(io, Int)
                i2 = readline_skip_text!(io, Int)
                i3 = readline_skip_text!(io, Int)
                i4 = readline_skip_text!(io, Int)
                @assert i1 == a1 "Quartet central index mismatch: got i1=$i1 for a1=$a1."

                lv1_frac = read_vec3!(io, T)
                lv2_frac = read_vec3!(io, T)
                lv3_frac = read_vec3!(io, T)
                lv4_frac = read_vec3!(io, T)

                ifcs = read_tensor4!(io)  # SArray{Tuple{3,3,3,3},T}

                # fractional positions including image shifts
                v1_frac = lv1_frac + r_frac_uc[a1]
                v2_frac = lv2_frac + r_frac_uc[i2]
                v3_frac = lv3_frac + r_frac_uc[i3]
                v4_frac = lv4_frac + r_frac_uc[i4]

                # Cartesian shifts and positions
                lv1_cart = SVector{3,T}(A * lv1_frac)
                lv2_cart = SVector{3,T}(A * lv2_frac)
                lv3_cart = SVector{3,T}(A * lv3_frac)
                lv4_cart = SVector{3,T}(A * lv4_frac)

                v1 = SVector{3,T}(A * v1_frac)
                v2 = SVector{3,T}(A * v2_frac)
                v3 = SVector{3,T}(A * v3_frac)
                v4 = SVector{3,T}(A * v4_frac)

                # relative vectors (Fortran: rv1=0)
                rv1 = SVector{3,T}(0,0,0)
                rv2 = chop3(v2 - v1)
                rv3 = chop3(v3 - v1)
                rv4_3 = chop3(v4 - v1)

                lv1_cart = chop3(lv1_cart)
                lv2_cart = chop3(lv2_cart)
                lv3_cart = chop3(lv3_cart)
                lv4_cart = chop3(lv4_cart)

                # cutoff per Fortran
                max_rcut = max(max_rcut, norm(rv2))
                max_rcut = max(max_rcut, norm(rv3))
                max_rcut = max(max_rcut, norm(rv4_3))
                max_rcut = max(max_rcut, norm(rv2 - rv3))
                max_rcut = max(max_rcut, norm(rv2 - rv4_3))
                max_rcut = max(max_rcut, norm(rv3 - rv4_3))

                # pack rv4 as SVector{4,T} to match your FC4Data field
                rv4 = SVector{4,T}(rv4_3[1], rv4_3[2], rv4_3[3], zero(T))

                quartets[i] = FC4Data{T}(
                    SVector{4,Int}(i1, i2, i3, i4),
                    SVector{4,SVector{3,T}}(lv1_cart, lv2_cart, lv3_cart, lv4_cart),
                    rv1, rv2, rv3, rv4,
                    ifcs
                )
            end

            push!(data, AtomFC4{T, n_quart}(SVector{n_quart}(quartets)))
        end

        r_cut = max_rcut + sqrt(eps(T))
        return IFCs{4, T}(na, r_cut, data)
    end
end
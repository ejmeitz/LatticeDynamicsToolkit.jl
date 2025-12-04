# lattice vecs are columns of cell matrix
to_frac_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L = mod.(cell \ position, L(1.0))

# lattice vecs are columns of cell matrix
to_cart_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L = cell * position

sqnorm(v::AbstractVector) = dot(v, v)

negsqrt(x::Real) = sign(x) * sqrt(abs(x))

is_gamma(q) = sqnorm(q) < lo_sqtol


 # Bose-Einstein distribution: n(ω,T) = 1/(exp(ℏω/k_B T) - 1)
function planck(temperature::Float64, omega::Float64)
   
    if temperature < lo_temperaturetol
        return 0.0
    end
    
    if omega < lo_freqtol
        return 0.0
    end
    
    x = omega / (kB_Hartree * temperature)
    
    if x > 1e2
        # Very large x → n ≈ 0 (avoid overflow)
        return 0.0
    else
        return 1.0 / (exp(x) - 1.0)
    end
end


# Temperature derivative of Bose-Einstein distribution: ∂n/∂T
function planck_deriv(temperature::Float64, omega::Float64)
    
    if temperature < lo_temperaturetol
        return 0.0
    end
    
    if omega < lo_freqtol
        return 0.0
    end
    
    x = omega / (kB_Hartree * temperature)
    n = planck(temperature, omega)
    
    return n * n * exp(x) * x / temperature
end


# Second temperature derivative of Bose-Einstein distribution: ∂²n/∂T²
function planck_secondderiv(temperature::Float64, omega::Float64)

    if temperature < lo_temperaturetol
        return 0.0
    end
    
    if omega < lo_freqtol
        return 0.0
    end
    
    x = omega / (kB_Hartree * temperature)
    ex = exp(x)
    
    return (ex * x * (2 + ex * (-2 + x) + x)) / ((ex - 1)^3 * temperature^2)
end

"""
Zeros out small complex parts, converts large complex parts to negative real number
"""
function clean_eigenvalue(freq_sq)
    if abs(imag(freq_sq)) > lo_sqtol
        return -abs(freq_sq)  # Mark as negative (imaginary freq)
    else
        return real(freq_sq)  # Take real part
    end
end

@inline function chop(x::Float64, tol::Float64)
    @inbounds for v in _WELLDEFINED_SMALL_64
        if abs(x - v) < tol
            return v
        end
    end
    return x
end

@inline function chop(x::ComplexF64, tol::Float64)
    re = real(x); im = imag(x)
    if abs(re) < tol; re = 0.0; end
    if abs(im) < tol; im = 0.0; end
    return ComplexF64(re, im)
end

function chop!(x::AbstractArray{T}, tol::Float64) where {T <: Union{Float64, ComplexF64}}
    @inbounds for i in eachindex(x)
        x[i] = chop(x[i], tol)
    end
    return x
end

chop3(v::SVector{3,T}, chop_tol::T) where T = SVector{3,T}(chop.(v, Ref(chop_tol)))


##########################################
# Stuff for reading forceconstant files #
##########################################

readline_skip_text!(io, T) = parse(T, first(split(strip(readline(io)))))


@inline function read_svec3!(io, ::Type{T}; conv = T(1.0)) where T
    xs = split(strip(readline(io)))
    return SVector{3,T}(conv .* parse.(T, xs))
end

@inline  function read_vec3!(io, out::AbstractVector{T}; conv = T(1.0)) where T
    out .= conv .* parse.(T, split(strip(readline(io))))
    return out
end

@inline function read_mat3_rows!(io, ::Type{T}; conv = T(1.0)) where T
    r1 = read_svec3!(io, T; conv = conv)
    r2 = read_svec3!(io, T; conv = conv)
    r3 = read_svec3!(io, T; conv = conv)
    M = hcat(r1, r2, r3)'   # r1 is first row, etc.
    return SMatrix{3, 3, T, 9}(M)     
end

@inline  function read_tensor3!(io, ::Type{T}; conv = T(1.0)) where T
    M = @MArray zeros(T, 3, 3, 3)
    for ii in 1:3, jj in 1:3
        @views read_vec3!(io, M[ii, jj, :]; conv = conv)
    end
    return SArray{Tuple{3,3,3}, T, 3, 27}(M)
end

@inline  function read_tensor4!(io, ::Type{T}; conv = T(1.0)) where T
    M = @MArray zeros(T, 3, 3, 3, 3)
    for ii in 1:3, jj in 1:3, kk in 1:3
        @views read_vec3!(io, M[ii, jj, kk, :]; conv = conv)
    end
    return SArray{Tuple{3,3,3,3}, T, 4, 81}(M)
end


##################################
# Stuff for reading POSCAR FILES #
##################################

function read_poscar_symbol_block(path::String)

    return open(path, "r") do f
        for _ in 1:5 
            readline(f)
        end
        species_line = readline(f)
        count_line = readline(f)

        symbols = Symbol.(split(strip(species_line)))
        counts = parse.(Int, split(strip(count_line)))
        (symbols, counts)
    end

end

function read_poscar_cell(path::String)

    return open(path, "r") do f
        readline(f)
        scale = parse(Float64, readline(f))
        lv1 = scale .* parse.(Float64, split(strip(readline(f))))
        lv2 = scale .* parse.(Float64, split(strip(readline(f))))
        lv3 = scale .* parse.(Float64, split(strip(readline(f))))

        cell = SMatrix{3,3,Float64}(hcat(lv1, lv2, lv3)) # cell vecs as columns

        readline(f) # skip species line

        n_atoms = sum(parse.(Int, split(strip(readline(f)))))

        (cell, n_atoms)
    end

end

function read_poscar_positions(path::String)

    cell, n_atoms = read_poscar_cell(path)

    x_frac = zeros(SVector{3, Float64}, n_atoms)

    # K = (convert_to_cart ? cell : Matrix{Float64}(I, 3, 3))
    parse_line = (line) -> SVector(parse.(Float64, split(strip(line))[1:3])...)

    open(path, "r") do f
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f) # skip species line
        readline(f) # natoms line
        coord_type = strip(readline(f)) # "direct coordinates" line
        
        if startswith(coord_type, 'C') || startswith(coord_type, "c")
            error("Cartesian POSCAR are not supported, please convert $(path) to fractional.")
        end

        for i in 1:n_atoms
            x_frac[i] = parse_line(readline(f))
        end
    end

    return x_frac, cell

end

function read_poscar_data(path::String)

    symbols, counts = read_poscar_symbol_block(path)
    species = reduce(vcat, [fill(s, c) for (s,c) in zip(symbols, counts)])
    positions, cell = read_poscar_positions(path)
    
    return species, positions, cell

end




"""
    clean_fractional_coordinate(x::T; tol::T = sqrt(eps(T))) -> T

Return a "clean" fractional coordinate `y` equivalent to `x` but snapped/wrapped into
the half-open interval `[0, 1)` using tolerance `tol`.

Behavior mirrors the Fortran `lo_clean_fractional_coordinates`:

- If `|y| < tol`, set `y = 0`.
- If `|y - 1| < tol`, set `y = 0`.
- If `y ≥ 1`, repeatedly subtract 1.
- If `y < 0`, repeatedly add 1.
- Stop when `-tol < y < 1`.
"""
@inline function clean_fractional_coordinates(x::T; tol::T = lo_sqtol) where {T<:AbstractFloat}
    y = x

    while true
        # snap very close to 0 and 1 to 0
        if abs(y) < tol
            y = zero(T)
        end
        if abs(y - one(T)) < tol
            y = zero(T)
        end

        # wrap into [0,1) by integer shifts of 1
        if y ≥ one(T)
            y -= one(T)
            continue
        end
        if y < zero(T)
            y += one(T)
            continue
        end

        # inside final window (-tol, 1) → done
        if (y > -tol) && (y < one(T))
            return y
        end
    end
end


"""
    gram_schmidt!(X)

Apply Gram-Schmidt orthogonalization to columns of complex matrix X.
"""
function gram_schmidt!(X::AbstractMatrix{ComplexF64})
    nr, nc = size(X)
    Q = copy(X)
    
    for k in 1:nc
        for i in 1:(k-1)
            Q_i = view(Q, :, i)
            Q_k = view(Q, :, k)
            R_ik = dot(Q_i, Q_k)  # conjugate dot product
            Q[:, k] .-= R_ik .* Q_i
        end
        R_kk = @views sqrt(real(dot(Q[:, k], Q[:, k])))
        if abs(R_kk) > lo_sqtol
            Q[:, k] ./= R_kk
        else
            Q[:, k] .= 0.0
        end
    end
    
    # Chop small values (matching lo_chop behavior)
    chop!(Q, 1e-13)
    X .= Q

    return X
end

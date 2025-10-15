# lattice vecs are columns of cell matrix
to_frac_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L = mod.(cell \ position, L(1.0))

# lattice vecs are columns of cell matrix
to_cart_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L = cell * position

sqnorm(v::AbstractVector) = dot(v, v)


##########################################
# Stuff for reading forceconstant files #
##########################################

readline_skip_text!(io, T) = parse(T, first(split(strip(readline(io)))))

chop3(v::SVector{3,T}, chop_tol::T) where T = SVector{3,T}(ntuple(i -> (abs(v[i]) < chop_tol ? zero(T) : v[i]), 3))


@inline function read_svec3!(io, ::Type{T}) where T
    xs = split(strip(readline(io)))
    return SVector{3,T}(parse.(T, xs))
end

@inline  function read_vec3!(io, out::AbstractVector{T}) where T
    out .= parse.(T, split(strip(readline(io))))
    return out
end

@inline  function read_tensor3!(io, ::Type{T}) where T
    M = @MArray zeros(T, 3, 3, 3)
    for ii in 1:3, jj in 1:3
        @views read_vec3!(io, M[ii, jj, :])
    end
    return SArray{Tuple{3,3,3},T}(M)
end

@inline  function read_tensor4!(io, ::Type{T}) where T
    M = @MArray zeros(T, 3, 3, 3, 3)
    for ii in 1:3, jj in 1:3, kk in 1:3
        @views read_vec3!(io, M[ii, jj, kk, :])
    end
    return SArray{Tuple{3,3,3,3},T}(M)
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

function read_poscar_cell(path::String, ::Type{T} = Float64) where T

    return open(path, "r") do f
        readline(f)
        scale = parse(T, readline(f))
        lv1 = scale .* parse.(T, split(strip(readline(f))))
        lv2 = scale .* parse.(T, split(strip(readline(f))))
        lv3 = scale .* parse.(T, split(strip(readline(f))))

        cell = SMatrix{3,3,T}(hcat(lv1, lv2, lv3)) # cell vecs as columns

        readline(f) # skip species line

        n_atoms = sum(parse.(Int, split(strip(readline(f)))))

        (cell, n_atoms)
    end

end

function read_poscar_positions(
        path,
        ::Type{T} = Float64;
        ssposcar_is_frac::Bool = true,
        store_frac_coords::Bool = false,
    ) where {T <: AbstractFloat}

    cell, n_atoms = read_poscar_cell(path, T)

    positions = zeros(SVector{3, T}, n_atoms)

    convert_to_cart = (!store_frac_coords && ssposcar_is_frac)

    K = (convert_to_cart ? cell : Matrix{T}(I, 3, 3))
    parse_line = (line) -> SVector(K * parse.(T, split(strip(line))[1:3])...)



    open(path, "r") do f
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f) # skip species line
        readline(f) # natoms line
        readline(f) # skip "direct coordinates" line

        for i in 1:n_atoms
            positions[i] = parse_line(readline(f))
        end
    end

    return positions, cell

end

function read_poscar_data(
        path, 
        ::Type{FLOAT_TYPE} = Float64; 
        ssposcar_is_frac::Bool = true,
        store_frac_coords::Bool = true
    ) where {FLOAT_TYPE <: AbstractFloat}

    symbols, counts = read_poscar_symbol_block(path)
    species = reduce(vcat, [fill(s, c) for (s,c) in zip(symbols, counts)])

    positions, cell = read_poscar_positions(path, FLOAT_TYPE;
                                            ssposcar_is_frac = ssposcar_is_frac,
                                            store_frac_coords = store_frac_coords)
    
    cell = FLOAT_TYPE.(cell)

    return species, positions, cell

end


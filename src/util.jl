# lattice vecs are columns of cell matrix
to_frac_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L = mod.(cell \ position, L(1.0))

# lattice vecs are columns of cell matrix
to_cart_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L = cell * position

# Helper functions
readline_skip_text!(io, T) = parse(T, first(split(strip(readline(io)))))

function read_vec3!(io, T)
    xs = split(strip(readline(io)))
    @assert length(xs) == 3 "Expected 3 components for a 3-vector."
    SVector{3,T}(parse.(T, xs))
end

function read_poscar_symbol_block(path::String)
    species_line = ""
    count_line = ""
    open(path, "r") do f
        for _ in 1:5 
            readline(f)
        end
        species_line = readline(f)
        count_line = readline(f)
    end
    symbols = Symbol.(split(strip(species_line)))
    counts = parse.(Int, split(strip(count_line)))

    return symbols, counts
end

function read_poscar_cell(path::String; n_atoms = nothing)

    cell = zeros(Float64, 3, 3)

    open(path, "r") do f
        readline(f)
        scale = parse(Float64, readline(f))
        lv1 = scale .* parse.(Float64, split(strip(readline(f))))
        lv2 = scale .* parse.(Float64, split(strip(readline(f))))
        lv3 = scale .* parse.(Float64, split(strip(readline(f))))

        cell .= hcat(lv1, lv2, lv3) # cell vecs as columns

        readline(f) # skip species line

        natoms_file = sum(parse.(Int, split(strip(readline(f)))))
        if !isnothing(n_atoms) && natoms_file != n_atoms
            error(ArgumentError("Poscar has $(natoms_file) but you told me it would have $(natoms)"))
        end
        n_atoms = natoms_file

    end

    return cell, n_atoms

end

function read_poscar_positions(
        path;
        n_atoms = nothing, 
        ssposcar_is_frac::Bool = true,
        store_frac_coords::Bool = false,
        FT::Type{T} = Float64
    ) where {T <: AbstractFloat}

    cell, n_atoms = read_poscar_cell(path; n_atoms = n_atoms)

    positions = zeros(SVector{3, T}, n_atoms)

    convert_to_cart = (!store_frac_coords && ssposcar_is_frac)

    if convert_to_cart
        parse_line = (line) -> SVector(cell * parse.(T, split(strip(line))[1:3])...)
    else
        parse_line = (line) -> SVector(parse.(T, split(strip(line))[1:3])...)
    end

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
        path; 
        n_atoms = nothing,
        ssposcar_is_frac::Bool = true,
        store_frac_coords::Bool = true,
        FT::Type{FLOAT_TYPE} = Float64
    ) where {FLOAT_TYPE <: AbstractFloat}

    symbols, counts = read_poscar_symbol_block(path)
    species = reduce(vcat, [fill(s, c) for (s,c) in zip(symbols, counts)])

    positions, cell = read_poscar_positions(path;
                                            n_atoms = n_atoms, 
                                            ssposcar_is_frac = ssposcar_is_frac,
                                            store_frac_coords = store_frac_coords,
                                            FT = FLOAT_TYPE)
    
    cell = FLOAT_TYPE.(cell)

    return species, positions, cell

end


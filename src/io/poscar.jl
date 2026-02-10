# POSCAR / VASP structure format parsing

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
    species = reduce(vcat, [fill(s, c) for (s, c) in zip(symbols, counts)])
    positions, cell = read_poscar_positions(path)

    return species, positions, cell
end

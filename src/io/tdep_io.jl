# TDEP-format positions and forces (infile.positions, infile.forces)

export read_positions, read_forces

"""
    read_positions(path::String, crystal::CrystalStructure, n_timesteps::Int) -> Matrix{Float64}

Read TDEP-format positions file (fractional coordinates, 3 columns, all timesteps concatenated).

Returns Cartesian positions in **bohr** as a Matrix of size `(3, na * n_timesteps)`.
Requires `crystal` to convert fractional → Cartesian.
"""
function read_positions(path::String, crystal::CrystalStructure, n_timesteps::Int)
    na = length(crystal)
    n_total = na * n_timesteps

    positions = Matrix{Float64}(undef, 3, n_total)

    open(path, "r") do io
        for i in 1:n_total
            line = readline(io)
            frac = parse.(Float64, split(strip(line)))
            # Convert fractional to Cartesian (in bohr, since crystal.L is in bohr)
            positions[:, i] = crystal.L * SVector{3,Float64}(frac)
        end
    end

    return positions
end

"""
    read_forces(path::String, na::Int, n_timesteps::Int) -> Matrix{Float64}

Read TDEP-format forces file (3 columns in eV/Å, all timesteps concatenated).

Returns forces in **Hartree/bohr** as a Matrix of size `(3, na * n_timesteps)`.
"""
function read_forces(path::String, na::Int, n_timesteps::Int)
    n_total = na * n_timesteps

    # Conversion: eV/Å → Hartree/bohr
    conv = eV_to_Hartree / A_to_bohr

    forces = Matrix{Float64}(undef, 3, n_total)

    open(path, "r") do io
        for i in 1:n_total
            line = readline(io)
            f_eVA = parse.(Float64, split(strip(line)))
            forces[:, i] = conv .* f_eVA
        end
    end

    return forces
end

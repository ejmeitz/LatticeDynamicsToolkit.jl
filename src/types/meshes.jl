export 
    AbstractQMesh, IBZMesh, FFTMesh, QPointFull,
    fft_third_grid_index, fft_fourth_grid_index,
    n_full_q_point, n_irr_point

# Abstract q-mesh type
abstract type AbstractQMesh{I <: Integer} end

# Lightweight IBZ mesh - just irreducible points
# Good for: band structures, DOS, properties that only need IBZ
struct IBZMesh{I <: Integer} <: AbstractQMesh{I}
    mesh::SVector{3, I}
    k_ibz::Vector{SVector{3, Float64}}
    weights::Vector{Float64}
    radius::Float64 # no integration weight included
    n_atoms_prim::Integer
end

# Full q-point data for a single point on the mesh
struct QPointFull
    r::SVector{3, Float64}  # fractional coordinates
    weight::Float64  # integration weight
    irreducible_index::Int  # which IBZ point this came from
end

# FFT mesh with full BZ data - needed for anharmonic calculations
# Has both irreducible points AND full mesh for momentum conservation
struct FFTMesh{I <: Integer} <: AbstractQMesh{I}
    mesh::SVector{3, I}  # grid dimensions (Nx, Ny, Nz)
    k_ibz::Vector{SVector{3, Float64}}  # irreducible points
    weights_ibz::Vector{Float64}  # integration weights for IBZ points
    full_index_ibz::Vector{Int}  # full mesh index for each IBZ point
    k_full::Vector{QPointFull}  # ALL q-points on the grid
    radius::Float64
    n_atoms_prim::Int
end

# Common interface for both mesh types
n_full_q_point(mesh::AbstractQMesh) = prod(mesh.mesh)
n_irr_point(mesh::AbstractQMesh) = length(mesh.k_ibz)
Base.length(mesh::AbstractQMesh) = length(mesh.k_ibz)
Base.eachindex(mesh::IBZMesh) = eachindex(mesh.weights)
Base.eachindex(mesh::FFTMesh) = eachindex(mesh.weights_ibz)

# Add IBZ constructor from Spglib using the types defined here
function Spglib.BrillouinZoneMesh(uc::CrystalStructure, mesh; symprec = 1e-5)

    cell = Spglib.SpglibCell(
        uc.L,
        uc.x_frac,
        Int.(AtomsBase.atomic_number(uc, :))
    )

    return Spglib.get_ir_reciprocal_mesh(cell, mesh, symprec)
end

function _build_ibz_mesh(uc::CrystalStructure, mesh; symprec = 1e-5)
     bzm = Spglib.BrillouinZoneMesh(uc, mesh; symprec = symprec)
    
    N = length(bzm.grid_address)

    # Indices of irreducible k-points (sorted for stable downstream use)
    ir_idx = sort!(unique(bzm.ir_mapping_table))

    if length(ir_idx) > 0.5 * length(bzm.grid_address)
        @warn "The number of irreducible k-points ($(length(ir_idx))) is more than half the total number of k-points ($(length(bzm.grid_address)))." *
              " This may indicate that the symmetry detection failed. Consider increasing symprec (currently $(symprec))."
    end

    counts = zeros(Int, N)
    @inbounds for i in 1:N
        j = bzm.ir_mapping_table[i]
        counts[j] += 1
    end
    multiplicity = counts[ir_idx]

    # Get list of fractional IBZ points
    invmesh = 1.0 ./ Float64.(bzm.mesh)
    k_ibz = Vector{SVector{3,Float64}}(undef, length(ir_idx))
    @inbounds for (p, j) in pairs(ir_idx)
        k_ibz[p] = (Float64.(bzm.grid_address[j])) .* invmesh
    end

    weights = multiplicity ./ sum(multiplicity)

    @assert length(k_ibz) == length(weights) "something wrong in IBZMesh construction"

    # cbrt(3 / (4*pi*V*Nk))
    radius = cbrt(3.0/primitive_volume(uc)/prod(mesh)/4.0/pi)

    N_prim = n_atoms_primitive(uc)

    return bzm, k_ibz, weights, radius, N_prim
end

function IBZMesh(uc::CrystalStructure, mesh; symprec = 1e-5)
    bzm, k_ibz, weights, radius, N_prim = _build_ibz_mesh(uc, mesh; symprec = symprec)
    return IBZMesh(bzm.mesh, k_ibz, weights, radius, N_prim)
end

"""
    FFTMesh(ibz::IBZMesh, bzm::Spglib.BrillouinZoneMesh)

Construct an FFTMesh from an IBZMesh by expanding to the full Brillouin zone.
Requires the original Spglib.BrillouinZoneMesh object for the full mesh mapping.
"""
function FFTMesh(uc::CrystalStructure, mesh; symprec = 1e-5)

    bzm, k_ibz, weights, radius, N_prim = _build_ibz_mesh(uc, mesh; symprec = symprec)
    ibz = IBZMesh(bzm.mesh, k_ibz, weights, radius, N_prim)

    Nx, Ny, Nz = ibz.mesh
    N_full = Nx * Ny * Nz
    
    # Generate all q-points on the full grid
    k_full = Vector{QPointFull}(undef, N_full)
    invmesh = 1.0 ./ Float64.(ibz.mesh)
    
    # Build map from irreducible index to full index
    ir_to_full = Dict{Int, Int}()
    for (full_idx, grid_addr) in enumerate(bzm.grid_address)
        ir_idx = bzm.ir_mapping_table[full_idx]
        if !haskey(ir_to_full, ir_idx)
            ir_to_full[ir_idx] = full_idx
        end
    end
    
    # Construct full q-point array with mapping to IBZ
    ir_sorted = sort(unique(bzm.ir_mapping_table))
    ir_to_ibz = Dict(ir_sorted[i] => i for i in 1:length(ir_sorted))
    
    for (full_idx, grid_addr) in enumerate(bzm.grid_address)
        k = SVector{3, Float64}(Float64.(grid_addr) .* invmesh)
        ir_idx = bzm.ir_mapping_table[full_idx]
        ibz_idx = ir_to_ibz[ir_idx]
        weight = 1.0 / N_full  # uniform weight for full mesh points
        k_full[full_idx] = QPointFull(k, weight, ibz_idx)
    end
    
    # Find which full index each IBZ point corresponds to
    full_index_ibz = [ir_to_full[ir_sorted[i]] for i in 1:length(ir_sorted)]
    
    return FFTMesh(
        ibz.mesh, 
        ibz.k_ibz, 
        ibz.weights, 
        full_index_ibz,
        k_full, 
        ibz.radius, 
        ibz.n_atoms_prim
    )
end

"""
    singlet_to_triplet(l, ny, nz)

Convert linear index (1-based) to 3D grid indices (i,j,k).
"""
function singlet_to_triplet(l::Int, ny::Int, nz::Int)
    k = mod(l - 1, nz) + 1
    j = div(mod(l - 1, ny * nz), nz) + 1
    i = div(l - 1, ny * nz) + 1
    return (i, j, k)
end

"""
    triplet_to_singlet(gi, ny, nz)

Convert 3D grid indices (i,j,k) to linear index (1-based).
"""
function triplet_to_singlet(gi::NTuple{3, Int}, ny::Int, nz::Int)
    return (gi[1] - 1) * ny * nz + (gi[2] - 1) * nz + gi[3]
end

"""
    fft_third_grid_index(i1, i2, dims)

Returns the index on the grid that gives q3 = -q1 - q2 (momentum conservation).
For three-phonon processes: q1 + q2 + q3 = 0 (mod G).
"""
function fft_third_grid_index(i1::Int, i2::Int, dims::NTuple{3, Int})
    gi1 = singlet_to_triplet(i1, dims[2], dims[3])
    gi2 = singlet_to_triplet(i2, dims[2], dims[3])
    
    # Momentum conservation: q3 = -q1 - q2
    # In grid indices (1-based): 3 is the offset for "zero"
    gi3 = (3 .- gi1 .- gi2)
    
    # Periodic boundary conditions
    gi3 = ntuple(i -> begin
        val = gi3[i]
        while val < 1
            val += dims[i]
        end
        while val > dims[i]
            val -= dims[i]
        end
        val
    end, 3)
    
    return triplet_to_singlet(gi3, dims[2], dims[3])
end

"""
    fft_fourth_grid_index(i1, i2, i3, dims)

Returns the index on the grid that gives q4 = -q1 - q2 - q3 (momentum conservation).
For four-phonon processes: q1 + q2 + q3 + q4 = 0 (mod G).
"""
function fft_fourth_grid_index(i1::Int, i2::Int, i3::Int, dims::NTuple{3, Int})
    gi1 = singlet_to_triplet(i1, dims[2], dims[3])
    gi2 = singlet_to_triplet(i2, dims[2], dims[3])
    gi3 = singlet_to_triplet(i3, dims[2], dims[3])
    
    # Momentum conservation: q4 = -q1 - q2 - q3
    # In grid indices (1-based): 4 is the offset for "zero"
    gi4 = (4 .- gi1 .- gi2 .- gi3)
    
    # Periodic boundary conditions
    gi4 = ntuple(i -> begin
        val = gi4[i]
        while val < 1
            val += dims[i]
        end
        while val > dims[i]
            val -= dims[i]
        end
        val
    end, 3)
    
    return triplet_to_singlet(gi4, dims[2], dims[3])
end

######################

# Smaller version of dispersion data that
# does not have eigenvector or full q-grid info
struct DispersionDataSimple{N} # N = N branch
    ibz::IBZMesh
    freqs::Vector{SVector{N, Float64}}  # each entry sorted smallest --> largest
    vels::Vector{SMatrix{3, N, Float64}} # Length is Nq
end 

# Constructor in dispersion.jl
######################


"""
    PhononDispersionPoint{N}

Phonon data at a single q-point. N = number of modes (3 * n_atoms).
"""
struct PhononDispersionPoint{N}
    omega::SVector{N, Float64}  # frequencies (angular, atomic units)
    vel::SMatrix{3, N, Float64, 3N}  # group velocities
    egv::SMatrix{N, N, ComplexF64, N*N}  # eigenvectors (column b is mode b)
end

"""
    PhononDispersions{N}

Full phonon dispersion data on irreducible and full q-meshes.
N = number of modes (3 * n_atoms).
"""
struct PhononDispersions{N}
    n_mode::Int  # number of phonon modes
    iq::Vector{PhononDispersionPoint{N}}  # irreducible q-points
    aq::Vector{PhononDispersionPoint{N}}  # all/full q-points
    default_smearing::Vector{Float64}  # default smearing per band
end

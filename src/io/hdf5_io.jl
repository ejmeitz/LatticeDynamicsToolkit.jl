# HDF5 read/write for amorphous force constants

export write_ifc2, read_ifc2_matrix

"""
    write_ifc2(path::String, ifc2::AmorphousIFC2; compress::Int=5)

Save AmorphousIFC2 to an HDF5 file. Stores the dense 3N×3N force constant matrix
with deflate compression to efficiently handle sparsity.

# Arguments
- `path`: Output HDF5 file path
- `ifc2`: AmorphousIFC2 to save
- `compress`: Compression level 0-9 (0=none, 9=max, default=5)

# Stored datasets (in eV/Å units for consistency with TDEP)
- `ifc2`: 3N×3N dense matrix (eV/Å²)
- `na`: Number of atoms
- `r_cut`: Cutoff radius (Å)
"""
function write_ifc2(path::String, ifc2::AmorphousIFC2; compress::Int=5)
    # Convert from internal units (Hartree/bohr²) to eV/Å²
    Φ = Matrix(ifc2) .* forceconstant_2nd_HartreeBohr_to_eVA

    h5open(path, "w") do f
        # Chunked + compressed dataset for the force constant matrix
        # Chunk size chosen for good compression of sparse-ish matrices
        n = size(Φ, 1)
        chunk_size = min(n, 256)  # reasonable chunk for most systems

        d = create_dataset(f, "ifc2", Float64, (n, n);
                          chunk=(chunk_size, chunk_size),
                          compress=compress)
        d[:, :] = Φ

        # Metadata (convert r_cut from bohr to Å)
        f["na"] = ifc2.na
        f["r_cut"] = ifc2.r_cut * bohr_to_A

        # Store units as attributes for clarity
        attrs(f)["units_ifc2"] = "eV/Ang^2"
        attrs(f)["units_r_cut"] = "Ang"
    end

    return nothing
end

"""
    read_ifc2_matrix(path::String) -> DenseIFC2

Read a dense force constant matrix from an HDF5 file saved by `write_ifc2`.

File is stored in eV/Å² units; returned DenseIFC2 is in internal units (Hartree/bohr²).

# Returns
- `DenseIFC2` containing the 3N×3N force constant matrix (Hartree/bohr²), na, and r_cut (bohr)
"""
function read_ifc2_matrix(path::String)
    h5open(path, "r") do f
        # Read and convert from eV/Å² to Hartree/bohr²
        Φ = read(f["ifc2"]) .* forceconstant_2nd_eVA_to_HartreeBohr
        na = read(f["na"])
        # Convert r_cut from Å to bohr
        r_cut = read(f["r_cut"]) * A_to_bohr
        return DenseIFC2(Φ, na, r_cut)
    end
end

"""
    read_ifc2(path::String, ::Type{AmorphousIFC2}; tol=1e-14) -> AmorphousIFC2

Read an HDF5 force constant file and return as sparse AmorphousIFC2 format.

# Arguments
- `path`: HDF5 file path (saved by `write_ifc2`)
- `tol`: Threshold for non-zero blocks (default 1e-14)

# Example
```julia
ifc2 = read_ifc2("forceconstants.h5", AmorphousIFC2)
```
"""
function read_ifc2(path::String, ::Type{AmorphousIFC2}; tol::Float64=1e-13)
    dense = read_ifc2_matrix(path)
    return AmorphousIFC2(dense; tol=tol)
end

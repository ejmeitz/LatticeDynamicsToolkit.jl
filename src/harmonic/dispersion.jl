export DOS

function DispersionDataSimple(
        uc::CrystalStructure,
        ifc2::IFC2,
        mesh;
        n_threads::Integer = Threads.nthreads()
    )
    ibz = SimpleMesh(uc, mesh)
    return DispersionDataSimple(uc, ifc2, ibz; n_threads = n_threads)
end

function DispersionDataSimple(
        uc::CrystalStructure,
        ifc2::IFC2,
        ibz::SimpleMesh; 
        n_threads::Integer = Threads.nthreads()
    )

    na = length(uc)
    nb = 3*na
    Nk = length(ibz.k_ibz) 

    # just allocate it all, not that much RAM
    freq_buf = zeros(SVector{nb, Float64}, Nk)
    vel_buf = zeros(SMatrix{3, nb, Float64}, Nk)

    # Get freqs at every q-point in IBZ
    @tasks for i in 1:Nk
        @set ntasks=n_threads

        @local begin
            ix = zeros(Int, nb)
            Dq = zeros(ComplexF64, nb, nb)
            ∂D∂q = zeros(ComplexF64, nb, nb, 3)
        end

        fill!(Dq, 0.0)
        fill!(∂D∂q, zero(ComplexF64))

        q_frac  = ibz.k_ibz[i]
        dynmat_and_derivative_q!(Dq, ∂D∂q, ifc2, uc, q_frac)
        freqs_sq, phi = get_modes(Hermitian(Dq), Val{is_gamma(q_frac)}())
        freqs_sq_real = clean_eigenvalue.(freqs_sq)
        freqs = SVector{nb, Float64}(negsqrt.(freqs_sq_real))

        sortperm!(ix, freqs)
        freqs = freqs[ix]
        phi = phi[:, ix]

        freq_buf[i] = freqs
        vel_buf[i] = group_velocities(freqs, phi, ∂D∂q)
    end

    return DispersionDataSimple{3*na}(ibz, freq_buf, vel_buf)

end
"""
    PhononDispersions(uc, ifc2, fft_mesh; n_threads, smearing_prefactor)

Construct PhononDispersions with dispersions on both irreducible and full q-meshes.

This constructor computes phonon dispersions (frequencies, group velocities, and 
eigenvectors) at all q-points in both the irreducible wedge and the full Brillouin zone.
This is required for anharmonic free energy calculations.

Uses the same functions as your existing DispersionDataSimple constructor:
`dynmat_and_derivative_q!`, `get_modes`, `group_velocities`, `is_gamma`, 
`clean_eigenvalue`, `negsqrt` (assumed to exist elsewhere).

# Arguments
- `uc::CrystalStructure`: Crystal structure
- `ifc2::IFC2`: Second order force constants
- `fft_mesh::FFTMesh`: FFT mesh (contains both IBZ and full mesh points)
- `n_threads::Integer`: Number of threads for parallel computation (default: all threads)
- `smearing_prefactor::Float64`: Prefactor for default smearing (default: 0.1)


# Notes
- The full mesh computation can be expensive for large meshes
- Default smearing is computed as: `smearing_prefactor * mean(|velocity|)` per band
- For acoustic modes with small velocities, uses frequency-based smearing instead
"""
function PhononDispersions(
    uc::CrystalStructure,
    ifc2::IFC2,
    fft_mesh::FFTMesh;
    n_threads::Integer = Threads.nthreads()
)
    
    na = length(uc)
    nb = 3 * na
    n_irr = n_irr_point(fft_mesh)
    n_full = n_full_q_point(fft_mesh)
    
    # Allocate arrays for irreducible points
    iq = Vector{PhononDispersionPoint{nb}}(undef, n_irr)
    
    # Compute dispersions at irreducible q-points
    @tasks for i in 1:n_irr
        @set ntasks = n_threads
        
        @local begin
            ix = zeros(Int, nb)
            Dq = zeros(ComplexF64, nb, nb)
            ∂D∂q = zeros(ComplexF64, nb, nb, 3)
        end
        
        fill!(Dq, 0.0)
        fill!(∂D∂q, zero(ComplexF64))
        
        q_frac = fft_mesh.k_ibz[i]
        dynmat_and_derivative_q!(Dq, ∂D∂q, ifc2, uc, q_frac)
        freqs_sq, phi = get_modes(Hermitian(Dq), Val{is_gamma(q_frac)}())
        freqs_sq_real = clean_eigenvalue.(freqs_sq)
        freqs = negsqrt.(freqs_sq_real)
        
        # Sort by frequency
        sortperm!(ix, freqs)
        omega = SVector{nb, Float64}(freqs[ix])
        egv_sorted = phi[:, ix]
        vels = group_velocities(freqs[ix], egv_sorted, ∂D∂q)
        
        # Store as PhononDispersionPoint
        iq[i] = PhononDispersionPoint{nb}(
            omega,
            SMatrix{3, nb, Float64}(vels),
            SMatrix{nb, nb, ComplexF64}(egv_sorted)
        )
    end
    
    # Allocate arrays for full mesh
    aq = Vector{PhononDispersionPoint{nb}}(undef, n_full)
    
    # Compute dispersions at all full q-points
    @tasks for i in 1:n_full
        @set ntasks = n_threads
        
        @local begin
            ix = zeros(Int, nb)
            Dq = zeros(ComplexF64, nb, nb)
            ∂D∂q = zeros(ComplexF64, nb, nb, 3)
        end
        
        fill!(Dq, 0.0)
        fill!(∂D∂q, zero(ComplexF64))
        
        q_frac = fft_mesh.k_full[i].r
        dynmat_and_derivative_q!(Dq, ∂D∂q, ifc2, uc, q_frac)
        freqs_sq, phi = get_modes(Hermitian(Dq), Val{is_gamma(q_frac)}())
        freqs_sq_real = clean_eigenvalue.(freqs_sq)
        freqs = negsqrt.(freqs_sq_real)
        
        # Sort by frequency
        sortperm!(ix, freqs)
        omega = SVector{nb, Float64}(freqs[ix])
        egv_sorted = phi[:, ix]
        vels = group_velocities(freqs[ix], egv_sorted, ∂D∂q)
        
        # Store as PhononDispersionPoint
        aq[i] = PhononDispersionPoint{nb}(
            omega,
            SMatrix{3, nb, Float64}(vels),
            SMatrix{nb, nb, ComplexF64}(egv_sorted)
        )
    end
    
    # Compute default smearing per band using frequency differences
    # Same approach as _default_smearing for DispersionDataSimple
    default_smearing = zeros(Float64, nb)
    for b in 1:nb
        f = 0.0
        freqs_band = [iq[i].omega[b] for i in 1:n_irr]
        sort!(freqs_band)
        for i in 1:(n_irr - 1)
            f = max(f, abs(freqs_band[i+1] - freqs_band[i]))
        end
        default_smearing[b] = f
    end
    
    # Help avoid issues when there are flat bands
    f = maximum(default_smearing) / 5.0
    default_smearing = max.(default_smearing, Ref(f))
    
    return PhononDispersions{nb}(nb, iq, aq, default_smearing)
end

function DOS(
        uc::CrystalStructure,
        ifc2::IFC2,
        mesh;
        n_threads::Integer = Threads.nthreads(),
        n_bins::Integer = 500,
        sigma_scale::Float64 = 1.0
    )

    ibz = SimpleMesh(uc, mesh)
    return DOS(
                uc, ifc2, ibz; 
                n_threads = n_threads,
                n_bins = n_bins,
                sigma_scale = sigma_scale
            )
end

function DOS(
        uc::CrystalStructure,
        ifc2::IFC2,
        ibz::SimpleMesh;
        n_threads::Integer = Threads.nthreads(),
        n_bins::Integer = 500,
        sigma_scale::Float64 = 1.0
    )

    k_ibz = getfield(ibz, :k_ibz)      # Vector{SVector{3,Float64}}
    w_ibz = getfield(ibz, :weights)    # Vector{Float64}
    @assert length(k_ibz) == length(w_ibz) "IBZ k-points and weights size mismatch"
    @assert abs(sum(w_ibz) - 1.0) < 1e-8 "IBZ weights should sum to 1 over the full BZ"

    na = length(uc)
    Nk = length(k_ibz) 

    dd = DispersionDataSimple(uc, ifc2, ibz; n_threads = n_threads)

    # -------- Build ω grid and Gaussian width --------
    max_freq = maximum(maximum.(dd.freqs))
    # Small pad so the rightmost Gaussian isn't truncated visually
    max_freq_pad = max_freq * 1.02

    σs_default = _default_smearing(dd)

    ω_grid = range(0.0, max_freq_pad; length=n_bins)
    Δω = step(ω_grid)
    ω_grid = collect(ω_grid)

    # -------- Pass 2: smear into DOS (threaded, windowed around each ωj) --------
    g = zeros(Float64, n_bins, n_threads)

    Threads.@threads for chunk_i in 1:n_threads
        for j in chunk_i:n_threads:Nk
        
            for (k, ωj) in enumerate(dd.freqs[j])

                σ = @views _adaptive_sigma(ibz.radius, dd.vels[j][:,k], σs_default[k], sigma_scale)
                σ² = σ*σ

                # restrict to ~4σ window for speed
                lo = ωj - 4σ
                hi = ωj + 4σ

                # find bin range (clamped)
                bmin = max(1, Int(cld(lo - ω_grid[1], Δω)) + 1)   # first bin with center >= lo
                bmax = min(n_bins, Int(fld(hi - ω_grid[1], Δω)) + 1)

                @inbounds for b in bmin:bmax
                    Δ = ω_grid[b] - ωj
                    # the extra divide by sigma from gaussian prefactor
                    g[b, chunk_i] += w_ibz[j] * exp(-0.5 * (Δ*Δ) / σ²) / (σ)
                end
            end
        end
    end

    # Reduce across threads
    g = vec(sum(g, dims=2)) ./ sqrt(2π)

    # -------- Normalize so ∫ g(ω) dω = 3N --------
    total = sum(g) * Δω
    if total > 0
        g .*= (3.0 * na) / total
    end

    return ω_grid, g
end

"""
    group_velocities(ω², U, ∂D∂q)


Inputs
---------
- `ω::AbstractVector{<:Real}`  # frequencies (sorted ascending)
- `U::Matrix{<:Complex} `      # eigenvectors as columns (nb×nb), orthonormal
- `∂D∂q::Array{ComplexF64,3}`  # ∂D/∂q, shape nb×nb×3 (Hermitian slices)

Returns
---------
- V :: Matrix{Float64} of size (3, nb), with V[α, i] = ∂ω_i/∂q_α
"""
function group_velocities(
        ω::AbstractVector{<:Real},
        U::Matrix{<:Complex},
        ∂D∂q::Array{ComplexF64,3}
    )

    nb = length(ω)

    # Detect subspaces (Fortran-style)
    subspaces, has_deg = find_degenerate_subspaces(ω)

    V_dwsq = zeros(Float64, 3, nb)  # d(ω²)/dq
    V      = @MMatrix zeros(Float64, 3, nb)  # dω/dq
    tmp    = similar(view(U, :, 1)) # Complex workspace (nb)

    if !has_deg
        # Non-degenerate path
        for i in 1:nb
            ui = @view U[:, i]
            for α in 1:3
                @views mul!(tmp, ∂D∂q[:,:,α], ui)       # tmp = (∂D/∂q_α) * e_i
                V_dwsq[α, i] = real(dot(ui, tmp))    # conj inner product
            end
        end
    else
        # Mixed: singleton subspaces use non-degenerate; larger ones use projected subspace
        for S in subspaces
            if length(S) == 1
                i = S[1]
                ui = @view U[:, i]
                for α in 1:3
                    @views mul!(tmp, ∂D∂q[:,:,α], ui)
                    V_dwsq[α, i] = real(dot(ui, tmp))
                end
            else
                US = @views U[:, S]     # nb×mb
                mb = size(US, 2)
                for α in 1:3
                    Hα = Hermitian(US' * (@view ∂D∂q[:,:,α]) * US)  # mb×mb
                    λα = eigen(Hα).values                         # real
                    avg = sum(λα) / mb
                    for i in S
                        V_dwsq[α, i] = avg
                    end
                end
            end
        end
    end

    # Convert dω²/dq -> dω/dq, force group velocities at Γ to be 0
    # Also forces group velocities of negative frequencies to 0
    for i in 1:nb
        ωi = ω[i]
        if ωi > lo_freqtol
            @views V[:, i] .= V_dwsq[:, i] ./ (2ωi)
        else
            @views V[:, i] .= 0.0
        end
    end

    # zero anything less than 1 nm / s 
    chop!(V, 1E-9/groupvel_Hartreebohr_to_ms)

    return SMatrix(V)
end


"""
    find_degenerate_subspaces(w2; freq_tol = lo_freqtol)

Given sorted frequencies ω (ascending), partition the mode indices `1:nb`
into contiguous subspaces where adjacent frequencies differ by less than `freq_tol`.

Returns
- subspaces :: Vector{Vector{Int}}   # e.g., [[1,2],[3],[4,5,6],...]
- has_degeneracies :: Bool           # true if any subspace has length > 1

"""
function find_degenerate_subspaces(ω::AbstractVector{<:Real})

    nb = length(ω)

    subspaces = Vector{Vector{Int}}()
    current = Int[]

    for i in 1:nb
        if isempty(current)
            push!(current, i)
        else
            if (ω[i] - ω[current[end]]) < lo_freqtol
                push!(current, i)
            else
                push!(subspaces, current)
                current = [i]
            end
        end
    end
    !isempty(current) && push!(subspaces, current)

    has_deg = any(length(S) > 1 for S in subspaces)
    return subspaces, has_deg
end

function _default_smearing(dd::DispersionDataSimple{N_BRANCH}) where N_BRANCH
    
    n_q_point = length(dd.freqs)
    σs = zeros(Float64, N_BRANCH)

    for i in 1:N_BRANCH
        f = 0.0
        freqs_band = [dd.freqs[j][i] for j in 1:n_q_point]
        sort!(freqs_band)
        for j in 1:(n_q_point - 1)
            f = max(f, abs(freqs_band[j+1] - freqs_band[j]))
        end
        σs[i] = f
    end

    # Copied from TDEP, help avoid issues when there are flat bands
    f = maximum(σs) / 5.0
    σs = max.(σs, Ref(f))
    return σs
end


const adaptive_sigma_largefactor = 4.0
const adaptive_sigma_smallfactor = 0.25
const adaptive_sigma_prefactor = 2*pi/sqrt(2.0)
function _adaptive_sigma(
        radius::Float64,
        gradient::SVector{3, Float64},
        default_σ::Float64,
        scale::Float64
    )
    σ = scale * adaptive_sigma_prefactor * radius * norm(gradient)
    σ = max(default_σ*adaptive_sigma_smallfactor*scale, σ)
    σ = min(default_σ*adaptive_sigma_largefactor*scale, σ)
    return σ
end



# function all_freqs(
#         uc::CrystalStructure,
#         ifc2::IFC2,
#         mesh;
#         n_threads::Integer = Threads.nthreads()
#     )

#     Nk = prod(mesh)
#     freq_buf = zeros(Float64, 3*length(uc), Nk) # just allocate it all, probably not that much RAM


#     @inline function q_from_linidx(idx::Int)
#         m1, m2, m3 = mesh
#         # 0-based indices in each dim
#         k0 = (idx-1) % m3
#         j0 = ((idx-1) ÷ m3) % m2
#         i0 = ((idx-1) ÷ (m2*m3)) % m1
#         q1 = (i0 / m1) - 0.5
#         q2 = (j0 / m2) - 0.5
#         q3 = (k0 / m3) - 0.5
#         return SVector{3,Float64}(q1, q2, q3)
#     end

#     # Get freqs at every q-point in IBZ
#     @tasks for i in 1:Nk
#         @set ntasks=n_threads

#         q_frac = q_from_linidx(i)
#         Dq = dynmat_q(ifc2, uc, q_frac)
#         f2, _ = get_modes(Dq, Val{is_gamma(q_frac)}())
#         freq_buf[:, i] .= sort!(sqrt.(f2))
#     end

#     return freq_buf
# end
export DOS

function _default_smearing(dd::DispersionData{N_BRANCH}) where N_BRANCH
    
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


"""
    group_velocities(ω², U, ∂D∂q)

Faithful port of the TDEP group-velocity computation.

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

function DispersionData(
        uc::CrystalStructure,
        ifc2::IFC2,
        mesh;
        n_threads::Integer = Threads.nthreads()
    )
    ibz = IBZMesh(uc, mesh)
    return DispersionData(uc, ifc2, ibz; n_threads = n_threads)
end

function DispersionData(
        uc::CrystalStructure,
        ifc2::IFC2,
        ibz::IBZMesh; 
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
        freqs = SVector{nb, Float64}(negsqrt.(freqs_sq))

        sortperm!(ix, freqs)
        freqs = freqs[ix]
        phi = phi[:, ix]

        freq_buf[i] = freqs
        vel_buf[i] = group_velocities(freqs, phi, ∂D∂q)
    end

    return DispersionData{3*na}(ibz, freq_buf, vel_buf)

end

function DOS(
        uc::CrystalStructure,
        ifc2::IFC2,
        mesh;
        n_threads::Integer = Threads.nthreads(),
        n_bins::Integer = 500,
        sigma_scale::Float64 = 1.0
    )

    ibz = IBZMesh(uc, mesh)
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
        ibz::IBZMesh;
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

    dd = DispersionData(uc, ifc2, ibz; n_threads = n_threads)

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
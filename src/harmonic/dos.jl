
#! NEED TO FIGURE OUT SMOOTHING PARAMETER BETTER

# This seems to oversmooth
function _default_smearing(freqs_buf_sorted::AbstractMatrix)
    
    n_branch, n_q_point = size(freqs_buf_sorted)
    σs = zeros(Float64, n_branch)

    for i in 1:n_branch
        f = 0.0
        for j in 1:(n_q_point - 1)
            f = max(f, abs(freqs_buf_sorted[i,j+1] - freqs_buf_sorted[i,j]))
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
    σ = max(default_σ*adaptive_sigma_smallfactor, σ)
    σ = min(default_σ*adaptive_sigma_largefactor, σ)
    return σ
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

    freq_buf = zeros(Float64, 3*na, Nk) # just allocate it all, probably not that much RAM

    # Get freqs at every q-point in IBZ
    @tasks for i in 1:Nk
        @set ntasks=n_threads
        q_frac  = k_ibz[i]
        Dq = dynmat_q(ifc2, uc, q_frac)
        f2, _ = get_modes(Dq, Val{is_gamma(q_frac)}())
        freq_buf[:, i] .= sort!(sqrt.(f2))
    end

    # -------- Build ω grid and Gaussian width --------
    max_freq = maximum(freq_buf)
    # Small pad so the rightmost Gaussian isn't truncated visually
    max_freq_pad = max_freq * 1.02

    #! TDEP HAS SOME ALGO FOR PICKING THIS
    σ = 0.01 * max_freq
    σ2 = σ * σ

    # σs = _default_smearing(freq_buf)
    # σs_sq = σs .* σs

    ω_grid = range(0.0, max_freq_pad; length=n_bins)
    Δω = step(ω_grid)
    ω_grid = collect(ω_grid)

    # -------- Pass 2: smear into DOS (threaded, windowed around each ωj) --------
    g = zeros(Float64, n_bins, n_threads)

    Threads.@threads for chunk_i in 1:n_threads
        for j in chunk_i:n_threads:Nk
        
            for (k, ωj) in enumerate(view(freq_buf, :, j))
                # σ = σs[k] 
                # σ2 = σs_sq[k]

                # restrict to ~4σ window for speed
                lo = ωj - 4σ
                hi = ωj + 4σ

                # find bin range (clamped)
                bmin = max(1, Int(cld(lo - ω_grid[1], Δω)) + 1)   # first bin with center >= lo
                bmax = min(n_bins, Int(fld(hi - ω_grid[1], Δω)) + 1)

                @inbounds for b in bmin:bmax
                    Δ = ω_grid[b] - ωj
                    g[b, chunk_i] += w_ibz[j] * exp(-0.5 * (Δ*Δ) / σ2) / (σ)
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

    return ω_grid, g, freq_buf
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
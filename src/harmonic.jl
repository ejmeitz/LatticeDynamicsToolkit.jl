export 
    dynmat_gamma,
    get_modes, 
    harmonic_properties

q_cart_from_frac(cell::CrystalStructure, q_frac::SVector{3,Float64}) = 2pi .* (cell.L_inv' * q_frac)

function dynmat_gamma(fc_sc::IFC2, sc::CrystalStructure)

    na = length(sc)
    nb = 3*na

    @assert na == fc_sc.na "Failed building dynmat. IFCs build on $(fc_sc.na) cell, but supercell has $(na) atoms"

    D = zeros(nb, nb)

    @inbounds for a1 in 1:na
        r1 = 3*(a1-1)
        w1 = sc.invsqrtm[a1]
        for pair in get_interactions(fc_sc, a1)
            a2 = pair.idxs[2]
            r2 = 3*(a2-1)
            D[r1+1:r1+3, r2+1:r2+3] .= pair.ifcs .* (w1 * sc.invsqrtm[a2])  
        end
    end

    # enforce exact symmetry
    @inbounds for j in 1:nb, i in j+1:nb
        s = 0.5 * (D[i,j] + D[j,i])
        D[i,j] = s; D[j,i] = s
    end

    return Hermitian(D)

end


function dynmat_q(
        fc_uc::IFC2,
        uc::CrystalStructure,
        q_frac::SVector{3,Float64};
    )

    na = length(uc)
    nb = 3*na

    @assert na == fc_uc.na "Failed building dynmat at q-point. IFCs build on $(fc_uc.na) cell, but unitcell has $(na) atoms"

    D_q = zeros(ComplexF64, nb, nb)

    q_cart = q_cart_from_frac(uc, q_frac) 

    @inbounds for a1 in 1:na
        inv_m_a1 = uc.invsqrtm[a1]
        r1 = 3*(a1 - 1)
        for fc in get_interactions(fc_uc, a1)        
            a2 = fc.idxs[2]
            r2 = 3*(a2-1)

            phase = cis(dot(q_cart, fc.lvs[2]))  # exp(im*(q · r))
            D_q[r1+1:r1+3, r2+1:r2+3] .+= fc.ifcs .* (phase * (inv_m_a1 * uc.invsqrtm[a2]))
        end
    end
    
    # Check its actually Hermetian
    c0 = 0.0
    for i in 1:nb
         for j in (i+1):nb
            c0 += D_q[i,j] - conj(D_q[j,i])
         end
    end
    if abs(c0) > lo_sqtol*nb
        @warn "Dynamical matrix at q=$(q_frac) is not Hermitian within tolerance, |Σ D_ij - D_ji*| = $(abs(c0))"
    end

    return Hermitian(D_q)
end

function get_modes(D::Hermitian, gamma_point::Val{true})
    eig_stuff = eigen(D)
    freqs_sq = eig_stuff.values
    idx_rt = sortperm(abs.(freqs_sq))
    @views freqs_sq[idx_rt[1:3]] .= 0.0
    return freqs_sq, eig_stuff.vectors
end

function get_modes(D::Hermitian, gamma_point::Val{false})
    eig_stuff = eigen(D)
    return eig_stuff.values, eig_stuff.vectors
end

function harmonic_properties(
        T, 
        uc::CrystalStructure,
        ifc2::IFC2,
        mesh,
        ::Type{L};
        n_threads::Integer = Threads.nthreads()
    ) where {L <: Limit}

    ibz = IBZMesh(uc, mesh)
    
    F₀ = sum_over_freqs(
        (ω) -> F_harmonic_single(ω, kB_Hartree*T, L), 
        ibz,
        uc,
        ifc2;
        n_threads = n_threads
    )

    S₀ = sum_over_freqs(
        (ω) -> S_harmonic_single(ω, kB_Hartree*T, L), 
        ibz,
        uc,
        ifc2;
        n_threads = n_threads
    )

    U₀ = sum_over_freqs(
        (ω) -> U_harmonic_single(ω, kB_Hartree*T, L), 
        ibz,
        uc,
        ifc2;
        n_threads = n_threads
    )

    Cᵥ₀ = sum_over_freqs(
        (ω) -> Cv_harmonic_single(ω, kB_Hartree*T, L), 
        ibz,
        uc,
        ifc2;
        n_threads = n_threads
    )

    #! IF THE USER GIVES NON-PRIMITIVE CELL, THIS MIGHT CAUSE ISSUES?
    #! SPGLIB GENERATES K-MESH FOR PRIMITIVE CELL
    N = length(uc)

    return F₀ / N, S₀ / N, U₀ / N, Cᵥ₀ / N
end


function harmonic_properties(T,  ω::AbstractVector, ::Type{L}) where {L <: Limit}
    N_atoms = Int(length(ω) / 3)
    F₀ = F_harmonic(ω, T, L)
    S₀ = S_harmonic(ω, T, L)
    U₀ = U_harmonic(ω, T, L)
    Cᵥ₀ = Cv_harmonic(ω, T, L)
    return F₀ / N_atoms, S₀ / N_atoms, U₀ / N_atoms, Cᵥ₀ / N_atoms
end


#* FIX CONTRIBUTION FROM Zero-Point MOTION ON RIGID TRANSLATION MODES??
#* SEE HOW TDEP IMPLEMENTS THINGS LIKE FREE ENERGY

function sum_over_freqs(freqs, f::Function; kwargs...)
    res = 0.0
    for freq in freqs
        if freq > lo_freqtol
            res += f(freq)
        end
    end
    return res
end

function sum_over_freqs(
        f, 
        uc::CrystalStructure,
        ifc2::IFC2,
        mesh;
        n_threads::Integer = Threads.nthreads()
    )

    @assert length(uc) == ifc2.na "Failed summing over freqs. IFCs build on $(ifc2.na) cell, but unitcell has $(length(uc)) atoms"

    ibz = IBZMesh(uc, mesh)

    return sum_over_freqs(f, ibz, uc, ifc2, n_threads = n_threads)
end

function sum_over_freqs(
        f, 
        ibz::IBZMesh,
        uc::CrystalStructure,
        ifc2::IFC2;
        n_threads::Integer = Threads.nthreads()
    )

    @assert length(uc) == ifc2.na "Failed summing over freqs. IFCs build on $(ifc2.na) cell, but unitcell has $(length(uc)) atoms"

    is_gamma(q) = sqnorm(q) < lo_sqtol

    res = @tasks for i in eachindex(ibz.weights)

        w = ibz.weights[i]
        q = ibz.k_ibz[i]

        @set begin
            ntasks=n_threads
            reducer=+
        end
        D_q = dynmat_q(ifc2, uc, q)
        freqs_sq, _ = get_modes(D_q, Val{is_gamma(q)}())
        freqs = sqrt.(freqs_sq)

        res_local = 0.0
        for freq in freqs
            if freq > lo_freqtol
                res_local += f(freq)
            end
        end

        w * res_local
    end
    
    return res
end

U_harmonic_single = (ω, kBT, ::Type{Quantum}) -> ω * ((1 / (exp(ω/(kBT)) - 1)) + 0.5)
U_harmonic_single = (ω, kBT, ::Type{Classical}) -> kBT

F_harmonic_single = (ω, kBT, ::Type{Quantum}) -> (0.5*ω) + (kBT * log(1 - exp(-ω/kBT)))
F_harmonic_single = (ω, kBT, ::Type{Classical}) -> kBT * log(ω/kBT)

S_harmonic_single = (ω, kBT, ::Type{Quantum}) -> kB_Hartree * (((ω/kBT) / (exp(ω/kBT) - 1)) - log(1 - exp(-ω/kBT)))
S_harmonic_single = (ω, kBT, ::Type{Classical}) -> kB_Hartree * (1 - log(ω/kBT))

Cv_harmonic_single = (ω, kBT, ::Type{Quantum}) -> kB_Hartree * ((ω/(2*kBT))^2) * (csch(ω/(2*kBT))^2)
Cv_harmonic_single = (ω, kBT, ::Type{Classical}) -> kB_Hartree

U_harmonic(freqs, T, ::Type{L}) where {L <: Limit} = sum_over_freqs(freqs, ω -> U_harmonic_single(ω, kB_Hartree*T, L))
F_harmonic(freqs, T, ::Type{L}) where {L <: Limit} = sum_over_freqs(freqs, ω -> F_harmonic_single(ω, kB_Hartree*T, L))
S_harmonic(freqs, T, ::Type{L}) where {L <: Limit} = sum_over_freqs(freqs, ω -> S_harmonic_single(ω, kB_Hartree*T, L))
Cv_harmonic(freqs, T, ::Type{L}) where {L <: Limit} = sum_over_freqs(freqs, ω -> Cv_harmonic_single(ω, kB_Hartree*T, L))

# function V_harmonic(ifc2::AbstractMatrix, u::AbstractVector)
#     return 0.5 * ((transpose(u) * ifc2) * u)
# end

# function U_harmonic(ω, T, ::Type{Quantum})
#     f = (freq) -> freq * ((1 / (exp(freq/(kB_Hartree*T)) - 1)) + 0.5)
#     return sum_over_freqs(ω, f)
# end

# function U_harmonic(ω, T, ::Type{Classical})
#     n_nonzero = count(freq -> freq > lo_freqtol, ω)
#     return n_nonzero*kB_Hartree*T
# end

# function F_harmonic(ω, T, ::Type{Quantum})
#     kBT = kB_Hartree * T
#     f = (freq) -> (0.5*freq) + kBT * log(1 - exp(-freq/kBT))
#     return sum_over_freqs(ω, f)
# end

# function F_harmonic(ω, T, ::Type{Classical})
#     kBT = kB_Hartree * T
#     f = (freq) -> log(freq/kBT)
#     return kBT * sum_over_freqs(ω, f)
# end

# function S_harmonic(ω, T, ::Type{Quantum})
#     kBT = kB_Hartree * T
#     f = (freq) -> ((freq/kBT) / (exp(freq/kBT) - 1)) - log(1 - exp(-freq/kBT))
#     return kB * sum_over_freqs(ω, f)
# end

# function S_harmonic(ω, T, ::Type{Classical})
#     kBT = kB_Hartree * T
#     f = (freq) -> (1 - log(freq/kBT))
#     return kB_Hartree * sum_over_freqs(ω, f)
# end

# function Cᵥ_harmonic(ω, T, ::Type{Quantum})
#     tkBT =  2 * kB_Hartree * T
#     f = (freq) -> ((freq/tkBT)^2) * (csch(freq/tkBT)^2)
#     return kB_Hartree * sum_over_freqs(ω, f)
# end

# function Cᵥ_harmonic(ω, T, ::Type{Classical})
#     n_nonzero = count(freq -> freq > lo_freqtol, ω)
#     return n_nonzero*kB_Hartree
# end
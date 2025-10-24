export dynmat_gamma, get_modes, harmonic_properties

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

function get_modes(D::Hermitian)
    eig_stuff = eigen(D)
    freqs_sq = eig_stuff.values
    idx_rt = sortperm(abs.(freqs_sq))
    freqs_sq[idx_rt[1:3]] .= 0.0
    return freqs_sq, eig_stuff.vectors
end


function harmonic_properties(T,  ω::AbstractVector, ::Type{L}) where {L <: Limit}
    F₀ = F_harmonic(ω, T, L)
    S₀ = S_harmonic(ω, T, L)
    U₀ = U_harmonic(ω, T, L)
    Cᵥ₀ = Cᵥ_harmonic(ω, T, L)
    return F₀, S₀, U₀, Cᵥ₀
end


#* FIX CONTRIBUTION FROM Zero-Point MOTION ON RIGID TRANSLATION MODES??
#* SEE HOW TDEP IMPLEMENTS THINGS LIKE FREE ENERGY

function V_harmonic(ifc2::AbstractMatrix, u::AbstractVector)
    return 0.5 * ((transpose(u) * ifc2) * u)
end

function sum_over_freqs(freqs, f::Function)
    res = 0.0
    for freq in freqs
        if freq > lo_freqtol
            res += f(freq)
        end
    end
    return res
end

function U_harmonic(ω, T, ::Type{Quantum})
    f = (freq) -> freq * ((1 / (exp(freq/(kB_Hartree*T)) - 1)) + 0.5)
    return sum_over_freqs(ω, f)
end

function U_harmonic(ω, T, ::Type{Classical})
    n_nonzero = count(freq -> freq > lo_freqtol, ω)
    return n_nonzero*kB_Hartree*T
end

function F_harmonic(ω, T, ::Type{Quantum})
    kBT = kB_Hartree * T
    f = (freq) -> (0.5*freq) + kBT * log(1 - exp(-freq/kBT))
    return sum_over_freqs(ω, f)
end

function F_harmonic(ω, T, ::Type{Classical})
    kBT = kB_Hartree * T
    f = (freq) -> log(freq/kBT)
    return kBT * sum_over_freqs(ω, f)
end

function S_harmonic(ω, T, ::Type{Quantum})
    kBT = kB_Hartree * T
    f = (freq) -> ((freq/kBT) / (exp(freq/kBT) - 1)) - log(1 - exp(-freq/kBT))
    return kB * sum_over_freqs(ω, f)
end

function S_harmonic(ω, T, ::Type{Classical})
    kBT = kB_Hartree * T
    f = (freq) -> (1 - log(freq/kBT))
    return kB_Hartree * sum_over_freqs(ω, f)
end

function Cᵥ_harmonic(ω, T, ::Type{Quantum})
    tkBT =  2 * kB_Hartree * T
    f = (freq) -> ((freq/tkBT)^2) * (csch(freq/tkBT)^2)
    return kB_Hartree * sum_over_freqs(ω, f)
end

function Cᵥ_harmonic(ω, T, ::Type{Classical})
    n_nonzero = count(freq -> freq > lo_freqtol, ω)
    return n_nonzero*kB_Hartree
end
export harmonic_properties

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

    # This is the cell used by Spglib when
    # building the IBZ, so this we use this 
    # to get the values per atom.
    N = ibz.n_atoms_prim

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


function sum_over_freqs_dense(
    f::F,
    mesh::SVector{3,Int},
    uc::CrystalStructure,
    ifc2::IFC2;
    n_threads::Integer = Threads.nthreads(),
) where {F}

    @assert length(uc) == ifc2.na "IFCs built on $(ifc2.na) atom cell, but unitcell has $(length(uc)) atoms"

    Nk = prod(mesh)

    # map a linear index -> (i,j,k) -> q (SVector{3,Float64})
    @inline function q_from_linidx(idx::Int)
        m1, m2, m3 = mesh
        # 0-based indices in each dim
        k0 = (idx-1) % m3
        j0 = ((idx-1) ÷ m3) % m2
        i0 = ((idx-1) ÷ (m2*m3)) % m1
        q1 = (i0 / m1) - 0.5
        q2 = (j0 / m2) - 0.5
        q3 = (k0 / m3) - 0.5
        return SVector{3,Float64}(q1, q2, q3)
    end

    res = @tasks for lin in 1:Nk
        @set begin
            ntasks = n_threads
            reducer = +
        end

        q = q_from_linidx(lin)

        D_q = dynmat_q(ifc2, uc, q)
        freqs_sq, _ = get_modes(D_q, Val{is_gamma(q)}())
        # Guard against tiny negative roundoff before sqrt
        freqs = sqrt.(max.(freqs_sq, 0.0))

        acc = 0.0
        @inbounds for ω in freqs
            if ω > lo_freqtol
                acc += f(ω)
            end
        end

        acc
    end

    return res / Nk / length(uc)
end

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
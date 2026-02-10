export 
    dynmat_gamma,
    get_modes,
    dynmat_q,
    dynmat_q!,
    dynmat_and_derivative_q,
    dynmat_and_derivative_q!

#! probably some speed to be gained by using SMatrix for dynmats

q_cart_from_frac(cell::CrystalStructure, q_frac::SVector{3,Float64}) = 2pi .* (cell.L_inv' * q_frac)

function precompute_ewald_parameters(ifc2::IFC2, uc::CrystalStructure)
    if !ifc2.has_polar_data
        return nothing
    end

    ew = EwaldParameters()
    λ = ifc2.polar.lambda > 0 ? ifc2.polar.lambda : 0.5
    set_ewald_parameters!(ew, uc, ifc2.polar.eps; strategy=3, lambda_forced=λ)
    return ew
end

function check_hermetian(D, nb; name = "Dynamical matrix")
    c0 = 0.0
    for i in 1:nb
         for j in (i+1):nb
            c0 += D[i,j] - conj(D[j,i])
         end
    end
    if abs(c0) > lo_sqtol*nb
        @warn "$(name) at q=$(q_frac) is not Hermitian within tolerance, |Σ D_ij - D_ji*| = $(abs(c0))"
    end
end

"""
    add_nonanalytic_gamma!(D_q, ifc2, uc, q_cart)

Add the non-analytic long-range dipole–dipole correction at Γ (LO–TO splitting)
to an in-place dynamical matrix `D_q`, using the polar data in `ifc2`.

This mirrors TDEP's `nonanalytical_dynamical_matrix`:

    D_{ij,αβ} = (q ⋅ Z*_{α,i}) (q ⋅ Z*_{β,j}) * 4π / (V * (qᵀ ε q))

where `q = q_cart / ‖q_cart‖`, `Z*` are Born effective charges, and `ε` is the
high-frequency dielectric tensor. `D_q` is assumed to be *mass-weighted*
already, so the contribution is scaled by `1/√(m_i m_j)` internally.
"""
function add_nonanalytic_gamma!(
        D_q::AbstractMatrix{ComplexF64},
        ifc2::IFC2,
        uc::CrystalStructure,
        q_cart::SVector{3,Float64},
    )

    # Require polar data and non-zero direction
    if !ifc2.has_polar_data
        return D_q
    end

    norm_q = sqrt(sqnorm(q_cart))
    norm_q < lo_sqtol && return D_q

    q = q_cart / norm_q

    na = ifc2.na
    nb = 3 * na

    eps   = ifc2.polar.eps
    bornZ = ifc2.polar.born_Z  # shape (3,3,na), matching TDEP: bornZ[α, :, a]

    V = volume(uc)  # bohr^3

    # qᵀ ε q
    den = dot(q, eps * q)
    if abs(den) < lo_sqtol
        return D_q
    end

    f0 = 4pi / (V * den)

    @inbounds for a1 in 1:na
        invm1 = uc.invsqrtm[a1]
        for a2 in 1:na
            invm2 = uc.invsqrtm[a2]
            mass_scale = invm1 * invm2

            # Precompute q·Z* for this atom pair
            zi = ntuple(i -> q[1]*bornZ[i,1,a1] + q[2]*bornZ[i,2,a1] + q[3]*bornZ[i,3,a1], 3)
            zj = ntuple(j -> q[1]*bornZ[j,1,a2] + q[2]*bornZ[j,2,a2] + q[3]*bornZ[j,3,a2], 3)

            r1 = 3*(a1-1)
            r2 = 3*(a2-1)

            for i in 1:3, j in 1:3
                Δ = zi[i] * zj[j] * f0 * mass_scale
                D_q[r1+i, r2+j] += ComplexF64(Δ, 0.0)
            end
        end
    end

    return D_q
end

"""
    add_longrange_ewald!(D_q, fc_uc, uc, q_frac; ewald=nothing, dDdq=nothing)

Add the full Ewald long-range dipole-dipole contribution to the in-place dynamical matrix
`D_q` at arbitrary q. The matrix is assumed mass-weighted; the long-range block is scaled
by 1/√(m_i m_j) before adding.

Use this for q ≠ Γ when polar data are present. At Γ, `add_nonanalytic_gamma!` is
preferred (direction-dependent limit).

- `ewald`: optional precomputed `EwaldParameters` for efficiency when sampling many q.
- `dDdq`: optional (nb,nb,3) array to accumulate ∂D/∂q_cart for group velocities etc.
"""
function add_longrange_ewald!(
        D_q::AbstractMatrix{ComplexF64},
        fc_uc::IFC2,
        uc::CrystalStructure,
        q_frac::SVector{3,Float64};
        ewald::Union{Nothing,EwaldParameters}=nothing,
        dDdq::Union{Nothing,AbstractArray{ComplexF64,3}}=nothing,
    )
    fc_uc.has_polar_data || return D_q

    ew = if ewald === nothing
        precompute_ewald_parameters(fc_uc, uc)
    else
        ewald
    end
    ew === nothing && return D_q

    na = length(uc)
    compute_grad = dDdq !== nothing

    if compute_grad
        D_lr = zeros(ComplexF64, 3, 3, na, na)
        Dx_lr = zeros(ComplexF64, 3, 3, na, na)
        Dy_lr = zeros(ComplexF64, 3, 3, na, na)
        Dz_lr = zeros(ComplexF64, 3, 3, na, na)
        longrange_dynamical_matrix_with_gradient!(D_lr, Dx_lr, Dy_lr, Dz_lr,
            ew, uc, q_frac, fc_uc.polar.born_Z, fc_uc.polar.eps)
    else
        D_lr = longrange_dynamical_matrix(ew, uc, q_frac, fc_uc.polar.born_Z, fc_uc.polar.eps)
    end

    @inbounds for a2 in 1:na, a1 in 1:na
        w = uc.invsqrtm[a1] * uc.invsqrtm[a2]
        r1, r2 = 3*(a1-1)+1, 3*(a2-1)+1
        for j in 1:3, i in 1:3
            D_q[r1+i-1, r2+j-1] += D_lr[i, j, a1, a2] * w
            if compute_grad
                dDdq[r1+i-1, r2+j-1, 1] += Dx_lr[i, j, a1, a2] * w
                dDdq[r1+i-1, r2+j-1, 2] += Dy_lr[i, j, a1, a2] * w
                dDdq[r1+i-1, r2+j-1, 3] += Dz_lr[i, j, a1, a2] * w
            end
        end
    end

    return D_q
end

function dynmat_gamma(fc_sc::IFC2, sc::CrystalStructure; include_polar::Bool=fc_sc.has_polar_data)

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

    # Add long-range dipole contribution at Γ when polar data present
    if include_polar && fc_sc.has_polar_data
        ew = precompute_ewald_parameters(fc_sc, sc)
        Φ_lr = supercell_longrange_forceconstant(ew, fc_sc.polar.born_Z, fc_sc.polar.eps, sc)
        @inbounds for a2 in 1:na, a1 in 1:na
            w = sc.invsqrtm[a1] * sc.invsqrtm[a2]
            r1, r2 = 3*(a1-1)+1, 3*(a2-1)+1
            D[r1:r1+2, r2:r2+2] .+= Φ_lr[:, :, a1, a2] .* w
        end
    end

    # enforce exact symmetry
    @inbounds for j in 1:nb, i in j+1:nb
        s = 0.5 * (D[i,j] + D[j,i])
        D[i,j] = s; D[j,i] = s
    end

    return Hermitian(D)

end


"""
    dynmat_gamma(ifc::AmorphousIFC2, crystal::CrystalStructure) -> Hermitian{Float64}

Build the dynamical matrix at the gamma point (q=0) for an amorphous system.

The dynamical matrix is: D_ij,αβ = Φ_ij,αβ / sqrt(m_i * m_j)

For amorphous systems, only the gamma point is meaningful since there's no 
translational symmetry / Brillouin zone.
"""
function dynmat_gamma(ifc::AmorphousIFC2, sc::CrystalStructure)
    na = ifc.na
    
    @assert na == length(sc) "IFC has $(na) atoms but crystal has $(length(sc))"
    
    # Get dense force constant matrix
    Φ = Matrix(ifc)
    
    return _mass_weight_ifc2_matrix_gamma(Φ, sc)
end

function _mass_weight_ifc2_matrix_gamma(Φ::AbstractMatrix, sc::CrystalStructure)
    na = length(sc)
    D = zeros(eltype(Φ), size(Φ))

    @inbounds for i in 1:na
        w_i = sc.invsqrtm[i]
        ri = 3*(i-1)
        for j in 1:na
            w_ij = w_i * sc.invsqrtm[j]
            rj = 3*(j-1)
            for α in 1:3
                for β in 1:3
                    D[ri+α, rj+β] = Φ[ri+α, rj+β] * w_ij
                end
            end
        end
    end
    
    return Hermitian(D)
end

function dynmat_q(
        fc_uc::IFC2,
        uc::CrystalStructure,
        q_frac::SVector{3,Float64};
        include_polar::Bool = fc_uc.has_polar_data,
        qdir_gamma::SVector{3,Float64} = SVector{3,Float64}(1.0, 0.0, 0.0),
        ewald::Union{Nothing,EwaldParameters}=nothing,
    )

    na = length(uc)
    nb = 3*na
    D_q = zeros(ComplexF64, nb, nb)

    return dynmat_q!(D_q, fc_uc, uc, q_frac; include_polar=include_polar, qdir_gamma=qdir_gamma, ewald=ewald)
end

function dynmat_q!(
        D_q::Matrix{ComplexF64},
        fc_uc::IFC2,
        uc::CrystalStructure,
        q_frac::SVector{3,Float64};
        include_polar::Bool = fc_uc.has_polar_data,
        qdir_gamma::SVector{3,Float64} = SVector{3,Float64}(1.0, 0.0, 0.0),
        ewald::Union{Nothing,EwaldParameters}=nothing,
    )

    na = length(uc)
    nb = 3*na

    @assert na == fc_uc.na "Failed building dynmat at q-point. IFCs build on $(fc_uc.na) cell, but unitcell has $(na) atoms"

    q_cart = q_cart_from_frac(uc, q_frac) 

    @inbounds for a1 in 1:na
        inv_m_a1 = uc.invsqrtm[a1]
        r1 = 3*(a1 - 1) + 1
        for fc in get_interactions(fc_uc, a1)        
            a2 = fc.idxs[2]
            r2 = 3*(a2-1) + 1

            phase = cis(dot(q_cart, fc.lvs[2]))  # exp(im*(q · r))
            D_q[r1:r1+2, r2:r2+2] .+= fc.ifcs .* (phase * (inv_m_a1 * uc.invsqrtm[a2]))
        end
    end

    # Polar long-range correction
    if include_polar && fc_uc.has_polar_data
        if is_gamma(q_frac)
            add_nonanalytic_gamma!(D_q, fc_uc, uc, q_cart_from_frac(uc, qdir_gamma))
        else
            add_longrange_ewald!(D_q, fc_uc, uc, q_frac; ewald=ewald)
        end
    end

    check_hermetian(D_q, nb)

    return D_q
end

function dynmat_and_derivative_q(
        fc_uc::IFC2,
        uc::CrystalStructure,
        q_frac::SVector{3,Float64};
        include_polar::Bool = fc_uc.has_polar_data,
        qdir_gamma::SVector{3,Float64} = SVector{3,Float64}(1.0, 0.0, 0.0),
        ewald::Union{Nothing,EwaldParameters}=nothing,
    )

    na = length(uc)
    nb = 3*na

    D_q = zeros(ComplexF64, nb, nb)
    ∂D∂q = zeros(ComplexF64, nb, nb, 3)

    return dynmat_and_derivative_q!(D_q, ∂D∂q, fc_uc, uc, q_frac; include_polar=include_polar, qdir_gamma=qdir_gamma, ewald=ewald)
end

function dynmat_and_derivative_q!(
        D_q::Matrix{ComplexF64},
        ∂D∂q::Array{ComplexF64, 3},
        fc_uc::IFC2,
        uc::CrystalStructure,
        q_frac::SVector{3,Float64};
        include_polar::Bool = fc_uc.has_polar_data,
        qdir_gamma::SVector{3,Float64} = SVector{3,Float64}(1.0, 0.0, 0.0),
        ewald::Union{Nothing,EwaldParameters}=nothing,
    )

    na = length(uc)
    nb = 3*na

    @assert na == fc_uc.na "Failed building dynmat at q-point. IFCs build on $(fc_uc.na) cell, but unitcell has $(na) atoms"

    q_cart = q_cart_from_frac(uc, q_frac) 

    @inbounds for a1 in 1:na
        inv_m_a1 = uc.invsqrtm[a1]
        r1 = 3*(a1 - 1) + 1
        for fc in get_interactions(fc_uc, a1)        
            a2 = fc.idxs[2]
            r2 = 3*(a2-1) + 1

            phase = cis(dot(q_cart, fc.lvs[2]))  # exp(im*(q · r))
            block = fc.ifcs .* (phase * (inv_m_a1 * uc.invsqrtm[a2]))
            D_q[r1:r1+2, r2:r2+2] .+= block
            ∂D∂q[r1:r1+2, r2:r2+2, 1] .+= (im * fc.lvs[2][1]) .* block
            ∂D∂q[r1:r1+2, r2:r2+2, 2] .+= (im * fc.lvs[2][2]) .* block
            ∂D∂q[r1:r1+2, r2:r2+2, 3] .+= (im * fc.lvs[2][3]) .* block
        end
    end

    # Polar long-range correction
    if include_polar && fc_uc.has_polar_data
        if is_gamma(q_frac)
            add_nonanalytic_gamma!(D_q, fc_uc, uc, q_cart_from_frac(uc, qdir_gamma))
        else
            add_longrange_ewald!(D_q, fc_uc, uc, q_frac; ewald=ewald, dDdq=∂D∂q)
        end
    end

    check_hermetian(D_q, nb)
    @views check_hermetian(∂D∂q[:,:,1], nb; name = "∂D∂q_x")
    @views check_hermetian(∂D∂q[:,:,2], nb; name = "∂D∂q_y")
    @views check_hermetian(∂D∂q[:,:,3], nb; name = "∂D∂q_z")

    return D_q, ∂D∂q
end

function get_modes(D, gamma_point::Val{true})
    eig_stuff = eigen(D)
    freqs_sq = eig_stuff.values
    idx_rt = sortperm(abs.(freqs_sq))
    @views freqs_sq[idx_rt[1:3]] .= 0.0
    return freqs_sq, eig_stuff.vectors
end

function get_modes(D, gamma_point::Val{false})
    eig_stuff = eigen(D)
    return eig_stuff.values, eig_stuff.vectors
end

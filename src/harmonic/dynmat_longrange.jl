# -----------------------------------------------------------------------------
# Long-range dynamical matrix (TDEP: longrange_dynamical_matrix)
# -----------------------------------------------------------------------------

"""
    add_nonanalytic_gamma!(D_q, ifc2, uc, q_cart)

Add the non-analytic long-range dipole-dipole correction at Γ (LO-TO splitting)
to an in-place dynamical matrix `D_q`, using the polar data in `ifc2`.

This mirrors TDEP's `nonanalytical_dynamical_matrix`:

    D_{ij,ab} = (q ⋅ Z*_{a,i}) (q ⋅ Z*_{b,j}) * 4π / (V * (qᵀ ε q))

where `q = q_cart / ||q_cart||`, `Z*` are Born effective charges, and `ε` is the
high-frequency dielectric tensor. `D_q` is assumed to be mass-weighted
already, so the contribution is scaled by `1/sqrt(m_i m_j)` internally.
"""
function add_nonanalytic_gamma!(
    D_q::AbstractMatrix{ComplexF64},
    ifc2::IFC2,
    uc::CrystalStructure,
    q_cart::SVector{3,Float64},
)
    if !ifc2.has_polar_data
        return D_q
    end

    norm_q = sqrt(sqnorm(q_cart))
    norm_q < lo_sqtol && return D_q
    q = q_cart / norm_q

    na = ifc2.na
    eps = ifc2.polar.eps
    bornZ = ifc2.polar.born_Z

    V = volume(uc)
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

            zi = ntuple(i -> q[1]*bornZ[i,1,a1] + q[2]*bornZ[i,2,a1] + q[3]*bornZ[i,3,a1], 3)
            zj = ntuple(j -> q[1]*bornZ[j,1,a2] + q[2]*bornZ[j,2,a2] + q[3]*bornZ[j,3,a2], 3)

            r1 = 3 * (a1 - 1)
            r2 = 3 * (a2 - 1)
            for i in 1:3, j in 1:3
                Δ = zi[i] * zj[j] * f0 * mass_scale
                D_q[r1 + i, r2 + j] += ComplexF64(Δ, 0.0)
            end
        end
    end

    return D_q
end

"""
    add_longrange_ewald!(D_q, fc_uc, uc, q_frac; ewald=nothing, dDdq=nothing)

Add the full Ewald long-range dipole-dipole contribution to the in-place dynamical matrix
`D_q` at arbitrary q. The matrix is assumed mass-weighted; the long-range block is scaled
by `1/sqrt(m_i m_j)` before adding.

Use this for q != Γ when polar data are present. At Γ, `add_nonanalytic_gamma!` is
preferred (direction-dependent limit).

- `ewald`: optional precomputed `EwaldParameters` for efficiency when sampling many q.
- `dDdq`: optional `(nb,nb,3)` array to accumulate `∂D/∂q_cart`.
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
    q_gamma = SVector{3,Float64}(0.0, 0.0, 0.0)

    # TDEP parity: include born_onsite_correction in long-range dynamical matrix.
    D_gamma = longrange_dynamical_matrix(
        ew, uc, q_gamma, fc_uc.polar.born_Z, fc_uc.polar.eps;
        reconly=true,
    )
    born_onsite = zeros(Float64, 3, 3, na)
    @inbounds for a1 in 1:na
        m = zeros(Float64, 3, 3)
        for a2 in 1:na
            m .+= real.(D_gamma[:, :, a1, a2])
        end
        born_onsite[:, :, a1] .= -m
    end

    if compute_grad
        D_lr = zeros(ComplexF64, 3, 3, na, na)
        Dx_lr = zeros(ComplexF64, 3, 3, na, na)
        Dy_lr = zeros(ComplexF64, 3, 3, na, na)
        Dz_lr = zeros(ComplexF64, 3, 3, na, na)
        longrange_dynamical_matrix_with_gradient!(
            D_lr, Dx_lr, Dy_lr, Dz_lr,
            ew, uc, q_frac, fc_uc.polar.born_Z, fc_uc.polar.eps;
            reconly=true, born_onsite=born_onsite,
        )
    else
        D_lr = longrange_dynamical_matrix(
            ew, uc, q_frac, fc_uc.polar.born_Z, fc_uc.polar.eps;
            reconly=true, born_onsite=born_onsite,
        )
    end

    @inbounds for a2 in 1:na, a1 in 1:na
        w = uc.invsqrtm[a1] * uc.invsqrtm[a2]
        r1, r2 = 3 * (a1 - 1) + 1, 3 * (a2 - 1) + 1
        for j in 1:3, i in 1:3
            D_q[r1 + i - 1, r2 + j - 1] += D_lr[i, j, a1, a2] * w
            if compute_grad
                dDdq[r1 + i - 1, r2 + j - 1, 1] += Dx_lr[i, j, a1, a2] * w
                dDdq[r1 + i - 1, r2 + j - 1, 2] += Dy_lr[i, j, a1, a2] * w
                dDdq[r1 + i - 1, r2 + j - 1, 3] += Dz_lr[i, j, a1, a2] * w
            end
        end
    end

    return D_q
end

"""
    longrange_dynamical_matrix!(D, ew, uc, q_frac, born_Z, eps; reconly=false, born_onsite=nothing)

Compute dipole-dipole long-range dynamical matrix D(q) in-place (no derivatives).
For D and ∂D/∂q_cart, use `longrange_dynamical_matrix_with_gradient!`.

- `reconly`: if `true`, only compute reciprocal-space part.
"""
function longrange_dynamical_matrix!(
    D::AbstractArray{ComplexF64,4},
    ew::EwaldParameters,
    uc::CrystalStructure,
    q_frac::SVector{3,Float64},
    born_Z::AbstractArray{<:Real,3},
    eps::AbstractMatrix;
    reconly::Bool=false,
    born_onsite::Union{Nothing,AbstractArray{<:Real,3}}=nothing,
)
    _longrange_dynamical_matrix_impl!(
        D, ew, uc, q_frac, born_Z, eps;
        born_onsite=born_onsite,
        reconly=reconly,
        chgmult=true,
        compute_grad=false,
        Dx=nothing,
        Dy=nothing,
        Dz=nothing,
    )
end

"""
    longrange_dynamical_matrix_with_gradient!(D, Dx, Dy, Dz, ew, uc, q_frac, born_Z, eps; reconly=false, born_onsite=nothing)

Compute dipole-dipole long-range dynamical matrix D(q) and its Cartesian derivatives
∂D/∂q_cart (returned in Dx, Dy, Dz) in-place.

- `reconly`: if `true`, only compute reciprocal-space part (TDEP uses `true`).
"""
function longrange_dynamical_matrix_with_gradient!(
    D::AbstractArray{ComplexF64,4},
    Dx::AbstractArray{ComplexF64,4},
    Dy::AbstractArray{ComplexF64,4},
    Dz::AbstractArray{ComplexF64,4},
    ew::EwaldParameters,
    uc::CrystalStructure,
    q_frac::SVector{3,Float64},
    born_Z::AbstractArray{<:Real,3},
    eps::AbstractMatrix;
    reconly::Bool=false,
    born_onsite::Union{Nothing,AbstractArray{<:Real,3}}=nothing,
)
    _longrange_dynamical_matrix_impl!(
        D, ew, uc, q_frac, born_Z, eps;
        born_onsite=born_onsite,
        reconly=reconly,
        chgmult=true,
        compute_grad=true,
        Dx=Dx,
        Dy=Dy,
        Dz=Dz,
    )
end

function _longrange_dynamical_matrix_impl!(
    D::AbstractArray{ComplexF64,4},
    ew::EwaldParameters,
    uc::CrystalStructure,
    q_frac::SVector{3,Float64},
    born_Z::AbstractArray{<:Real,3},
    eps::AbstractMatrix;
    born_onsite::Union{Nothing,AbstractArray{<:Real,3}},
    reconly::Bool,
    chgmult::Bool,
    compute_grad::Bool,
    Dx::Union{Nothing,AbstractArray{ComplexF64,4}},
    Dy::Union{Nothing,AbstractArray{ComplexF64,4}},
    Dz::Union{Nothing,AbstractArray{ComplexF64,4}},
)
    na = length(uc)
    @assert size(D) == (3, 3, na, na)
    @assert size(born_Z) == (3, 3, na)

    if compute_grad
        @assert Dx !== nothing && Dy !== nothing && Dz !== nothing
        @assert size(Dx) == size(Dy) == size(Dz) == size(D)
    end

    # Use shared conversion helper, then map to TDEP's internal no-2π basis.
    q_cart_twopi = q_cart_from_frac(uc, q_frac)       # includes 2π
    q_cart = q_cart_twopi / TWOPI                     # no-2π basis
    x_cart = uc.x_cart
    inveps = SMatrix{3,3,Float64}(inv(eps))
    inveps = SMatrix{3,3,Float64}(chop.(inveps, Ref(lo_sqtol)))
    dete = 1.0 / sqrt(det(eps))
    mult = chgmult
    if born_onsite !== nothing
        @assert size(born_onsite) == (3, 3, na)
    end
    onsite = something(born_onsite, zeros(3, 3, na))

    npair = na * (na + 1) ÷ 2
    ucvl = [chop3(x_cart[a2] - x_cart[a1], lo_sqtol) for a1 in 1:na for a2 in a1:na]

    inv4lam2 = 1.0 / (4.0 * ew.lambda^2)
    lambdacub = ew.lambda^3

    DL1 = [@SMatrix zeros(ComplexF64, 3, 3) for _ in 1:npair]
    DL2 = [@SMatrix zeros(ComplexF64, 3, 3) for _ in 1:npair]
    DL3 = [@SMatrix zeros(ComplexF64, 3, 3) for _ in 1:na]

    DQLx1 = compute_grad ? [@SMatrix zeros(ComplexF64, 3, 3) for _ in 1:npair] : nothing
    DQLy1 = compute_grad ? [@SMatrix zeros(ComplexF64, 3, 3) for _ in 1:npair] : nothing
    DQLz1 = compute_grad ? [@SMatrix zeros(ComplexF64, 3, 3) for _ in 1:npair] : nothing
    DQLx2 = compute_grad ? [@SMatrix zeros(ComplexF64, 3, 3) for _ in 1:npair] : nothing
    DQLy2 = compute_grad ? [@SMatrix zeros(ComplexF64, 3, 3) for _ in 1:npair] : nothing
    DQLz2 = compute_grad ? [@SMatrix zeros(ComplexF64, 3, 3) for _ in 1:npair] : nothing

    # Reciprocal space
    for ig in 1:size(ew.Gvec, 2)
        G = SVector{3,Float64}(ew.Gvec[1,ig], ew.Gvec[2,ig], ew.Gvec[3,ig])
        K = chop3((G + q_cart) * TWOPI, lo_sqtol)
        if sqnorm(K) < lo_sqtol^2
            continue
        end
        knorm = dot(K, eps * K)
        invknorm = 1.0 / knorm
        partialChi = exp(-knorm * inv4lam2) * invknorm

        Kx, Ky, Kz = K[1], K[2], K[3]
        kk = K * K'

        if compute_grad
            kkx = @SMatrix [2*Kx Ky Kz; Ky 0.0 0.0; Kz 0.0 0.0]
            kky = @SMatrix [0.0 Kx 0.0; Kx 2*Ky Kz; 0.0 Kz 0.0]
            kkz = @SMatrix [0.0 0.0 Kx; 0.0 0.0 Ky; Kx Ky 2*Kz]
            Keps_vec = SVector{3,Float64}(2 * (K' * eps)) .* (invknorm + inv4lam2)
        end

        for ip in 1:npair
            tau = -ucvl[ip]
            ikr = dot(K, tau)
            expikr = cis(ikr)
            Chi = partialChi * expikr
            DL1[ip] += kk * Chi

            if compute_grad
                v0 = im .* tau .- Keps_vec
                DQLx1[ip] += kk * v0[1] * Chi + kkx * Chi
                DQLy1[ip] += kk * v0[2] * Chi + kky * Chi
                DQLz1[ip] += kk * v0[3] * Chi + kkz * Chi
            end
        end
    end
    f0 = 4π / volume(uc)
    for ip in 1:npair
        DL1[ip] *= f0
        if compute_grad
            DQLx1[ip] *= f0
            DQLy1[ip] *= f0
            DQLz1[ip] *= f0
        end
    end

    # Real space
    if !reconly
        for ir in 1:size(ew.Rvec, 2)
            R = SVector{3,Float64}(ew.Rvec[1,ir], ew.Rvec[2,ir], ew.Rvec[3,ir])
            for ip in 1:npair
                Dvec = chop3(R + ucvl[ip], lo_sqtol)
                ikr = dot(q_cart_twopi, R)
                expikr = cis(ikr)
                delta = inveps * Dvec
                bigD = sqrt(dot(delta, Dvec))
                H = ewald_H_thingy(ew.lambda .* delta, ew.lambda * bigD, inveps)
                Hp = H * (lambdacub * expikr * dete)
                DL2[ip] += Hp
                if compute_grad
                    DQLx2[ip] += im * R[1] * Hp
                    DQLy2[ip] += im * R[2] * Hp
                    DQLz2[ip] += im * R[3] * Hp
                end
            end
        end

        for a1 in 1:na
            DL3[a1] = (ew.lambda^3) * inveps * dete * FOUR_OVER_THREE_SQRTPI
        end
    end

    # Assemble D0 and Dx0, Dy0, Dz0
    D0 = zeros(ComplexF64, 3, 3, na, na)
    Dx0 = compute_grad ? zeros(ComplexF64, 3, 3, na, na) : nothing
    Dy0 = compute_grad ? zeros(ComplexF64, 3, 3, na, na) : nothing
    Dz0 = compute_grad ? zeros(ComplexF64, 3, 3, na, na) : nothing

    ip = 0
    for a1 in 1:na, a2 in a1:na
        ip += 1
        if a1 == a2
            D0[:, :, a1, a2] = DL1[ip] - DL2[ip] - DL3[a1]
            if compute_grad
                Dx0[:, :, a1, a2] = DQLx1[ip] - DQLx2[ip]
                Dy0[:, :, a1, a2] = DQLy1[ip] - DQLy2[ip]
                Dz0[:, :, a1, a2] = DQLz1[ip] - DQLz2[ip]
            end
        else
            D0[:, :, a1, a2] = DL1[ip] - DL2[ip]
            if compute_grad
                Dx0[:, :, a1, a2] = DQLx1[ip] - DQLx2[ip]
                Dy0[:, :, a1, a2] = DQLy1[ip] - DQLy2[ip]
                Dz0[:, :, a1, a2] = DQLz1[ip] - DQLz2[ip]
            end
            for i in 1:3, j in 1:3
                D0[j, i, a2, a1] = conj(D0[i, j, a1, a2])
                if compute_grad
                    Dx0[j, i, a2, a1] = conj(Dx0[i, j, a1, a2])
                    Dy0[j, i, a2, a1] = conj(Dy0[i, j, a1, a2])
                    Dz0[j, i, a2, a1] = conj(Dz0[i, j, a1, a2])
                end
            end
        end
    end

    if mult
        D .= 0
        if compute_grad
            Dx .= 0
            Dy .= 0
            Dz .= 0
        end
        for a2 in 1:na, a1 in 1:na
            B1 = @view born_Z[:, :, a1]
            B2 = @view born_Z[:, :, a2]
            D[:, :, a1, a2] .= B1 * (@view D0[:, :, a1, a2]) * transpose(B2)
            if compute_grad
                Dx[:, :, a1, a2] .= B1 * (@view Dx0[:, :, a1, a2]) * transpose(B2)
                Dy[:, :, a1, a2] .= B1 * (@view Dy0[:, :, a1, a2]) * transpose(B2)
                Dz[:, :, a1, a2] .= B1 * (@view Dz0[:, :, a1, a2]) * transpose(B2)
            end
        end
        for a1 in 1:na
            D[:, :, a1, a1] += onsite[:, :, a1]
        end
    else
        D .= D0
        if compute_grad
            Dx .= Dx0
            Dy .= Dy0
            Dz .= Dz0
        end
    end

    chop!(D, lo_sqtol)
    return D
end

"""
    longrange_dynamical_matrix(ew, uc, q_frac, born_Z, eps; reconly=false, born_onsite=nothing) -> Array{ComplexF64,4}

Allocating version returning (3,3,na,na) complex dynamical matrix.
"""
function longrange_dynamical_matrix(
    ew::EwaldParameters,
    uc::CrystalStructure,
    q_frac::SVector{3,Float64},
    born_Z::AbstractArray,
    eps::AbstractMatrix;
    reconly::Bool=false,
    born_onsite::Union{Nothing,AbstractArray{<:Real,3}}=nothing,
)
    na = length(uc)
    D = zeros(ComplexF64, 3, 3, na, na)
    longrange_dynamical_matrix!(D, ew, uc, q_frac, born_Z, eps; reconly=reconly, born_onsite=born_onsite)
    return D
end

"""
    longrange_dynamical_matrix(ifc2::IFC2, uc::CrystalStructure, q_frac; reconly=false)

Convenience: use polar data from IFC2. Returns nothing if `has_polar_data` is false.
"""
function longrange_dynamical_matrix(
    ifc2::IFC2,
    uc::CrystalStructure,
    q_frac::SVector{3,Float64};
    reconly::Bool=false,
)
    ifc2.has_polar_data || return nothing
    ew = EwaldParameters()
    λ = ifc2.polar.lambda
    λ > 0 || throw(ArgumentError("For TDEP parity, IFC2 polar lambda must be > 0; got λ=$(λ)."))
    set_ewald_parameters!(ew, uc, ifc2.polar.eps; strategy=3, lambda_forced=λ)
    longrange_dynamical_matrix(ew, uc, q_frac, ifc2.polar.born_Z, ifc2.polar.eps; reconly=reconly)
end

# -----------------------------------------------------------------------------
# Supercell long-range force constants at Gamma (TDEP: supercell_longrange_forceconstant)
# -----------------------------------------------------------------------------

"""
    supercell_longrange_forceconstant(ew, fc2, sc::CrystalStructure, uc::CrystalStructure; thres=1e-8) -> Array{Float64,4}

Checked convenience overload using polar data from `IFC2`.
Requires remapped supercell IFC2 (`fc2.na == length(sc)`) and supercell-indexed
Born charges (`size(fc2.polar.born_Z,3) == length(sc)`).
"""
function supercell_longrange_forceconstant(
    ew::EwaldParameters,
    fc2::IFC2,
    sc::CrystalStructure,
    uc::CrystalStructure;
    thres::Float64=1e-8,
)
    fc2.has_polar_data || throw(ArgumentError("supercell_longrange_forceconstant(ew, fc2, sc, uc) requires polar IFC2 data."))
    fc2.na == length(sc) || throw(ArgumentError(
        "Expected IFC2 remapped to supercell: fc2.na=$(fc2.na), length(sc)=$(length(sc))."
    ))
    size(fc2.polar.born_Z, 3) == length(sc) || throw(ArgumentError(
        "Expected supercell-indexed Born charges in remapped IFC2: size(born_Z,3)=$(size(fc2.polar.born_Z,3)), length(sc)=$(length(sc))."
    ))
    return _sc_lr_ifc_impl(ew, fc2.polar.born_Z, fc2.polar.eps, sc, uc; thres=thres)
end

"""
    supercell_longrange_forceconstant(ew, born_Z, eps, sc::CrystalStructure, uc::CrystalStructure; thres=1e-8) -> Array{Float64,4}

Long-range contribution to force constants at Γ for supercell.
For TDEP-style supercell workflows, `born_Z` is expected to be supercell-indexed
with size `(3,3,na_sc)`.
Returns real (3,3,na,na) - add to short-range IFCs for full force constants.
"""
function _sc_lr_ifc_impl(
    ew::EwaldParameters,
    born_Z::AbstractArray,
    eps::AbstractMatrix,
    sc::CrystalStructure,
    uc::CrystalStructure;
    thres::Float64=1e-8,
)
    na = length(sc)
    # Keep reciprocal basis convention consistent with TDEP (no 2π here).
    L_rec = inv(sc.L)'
    inv4lam2 = 1.0 / (4.0 * ew.lambda^2)
    size(born_Z, 3) == na || error("For TDEP-style supercell flow, born_Z must be supercell-indexed with size(born_Z,3)=length(sc)=$na; got $(size(born_Z, 3)).")
    s2u = map_super_to_unitcell(uc, sc)

    # G-vectors in reciprocal space (Gamma: q=0)
    pts = _points_on_sphere(20)
    krad = inscribed_sphere_in_box(L_rec) * 0.5
    krad_converged = false
    for _ in 1:100000
        dampsum = 0.0
        for i in 1:size(pts,2)
            v = SVector{3,Float64}(pts[1,i], pts[2,i], pts[3,i]) * krad * TWOPI
            knorm = dot(v, eps * v)
            dampsum += exp(-knorm * inv4lam2) * knorm
        end
        if dampsum <= thres
            krad_converged = true
            break
        end
        krad += krad * 0.25
    end
    krad_converged || error("Could not converge krad in supercell_longrange_forceconstant")

    ndim = 0
    for i in 1:100
        ndim += 1
        m0 = [L_rec[:,j] * (2*ndim + 1) for j in 1:3]
        m0 = hcat(m0...)
        if inscribed_sphere_in_box(m0) > krad + lo_tol
            break
        end
    end

    k2 = krad^2
    qvecs = SVector{3,Float64}[]
    nvec_total = (2 * ndim + 1)^3
    nvec_stop = (nvec_total + 1) ÷ 2
    ctrtot = 0
    done = false
    for i in -ndim:ndim
        for j in -ndim:ndim
            for k in -ndim:ndim
                v = L_rec * SVector{3,Float64}(Float64(i), Float64(j), Float64(k))
                if sqnorm(v) < k2 && sqnorm(v) > lo_sqtol
                    push!(qvecs, v)
                end
                ctrtot += 1
                if ctrtot == nvec_stop
                    done = true
                    break
                end
            end
            done && break
        end
        done && break
    end
    nvec = length(qvecs)

    # Build delta vectors like TDEP: wrapped fractional deltas with unit-cell indices in key
    pair_keys = SVector{5,Float64}[]
    pair_delta_cart = SVector{3,Float64}[]
    for a1 in 1:na, a2 in a1:na
        dfrac = sc.x_frac[a2] - sc.x_frac[a1]
        dfrac = clean_fractional_coordinates.(dfrac .+ 0.5; tol=lo_tol) .- 0.5
        uca1 = s2u[a1]
        uca2 = s2u[a2]
        key = SVector{5,Float64}(dfrac[1], dfrac[2], dfrac[3], Float64(uca1), Float64(uca2))
        push!(pair_keys, key)
        push!(pair_delta_cart, sc.L * SVector{3,Float64}(dfrac[1], dfrac[2], dfrac[3]))
    end

    deltavecind = Int[]
    unique_keys = SVector{5,Float64}[]
    deltavec = SVector{3,Float64}[]
    for ip in eachindex(pair_keys)
        idx = findfirst(k -> sum(abs.(k - pair_keys[ip])) <= lo_tol, unique_keys)
        if idx === nothing
            push!(unique_keys, pair_keys[ip])
            push!(deltavec, pair_delta_cart[ip])
            push!(deltavecind, length(unique_keys))
        else
            push!(deltavecind, idx)
        end
    end

    # Reciprocal sum only (Gamma): K = q (Cartesian), dKK = (K⊗K) * exp(-K²/4λ²) / K²
    dKK = [let K = q  # q already Cartesian reciprocal
        knorm = dot(K, eps * K)
        f0 = exp(-knorm * inv4lam2) / knorm
        SMatrix{3,3,Float64}(K * K' * f0)
    end for q in qvecs]

    # exp(i K·r) with r in Cartesian; K·r is dimensionless
    rDL = [sum(dKK[qi] * cos(-dot(qvecs[qi], d)) for qi in 1:nvec) for d in deltavec]

    f0 = 2.0 * (4π / volume(sc))
    rDL = [m * f0 for m in rDL]

    # Map back to full (a1,a2)
    dm = zeros(Float64, 3, 3, na, na)
    ip = 0
    for a1 in 1:na, a2 in a1:na
        ip += 1
        pp = deltavecind[ip]
        dm[:, :, a1, a2] = rDL[pp]
        if a1 != a2
            for i in 1:3, j in 1:3
                dm[j, i, a2, a1] = dm[i, j, a1, a2]
            end
        end
    end

    chop!(dm, 1e-13)

    # Multiply in Born charges
    forceconstant = zeros(Float64, 3, 3, na, na)
    for a2 in 1:na, a1 in 1:na
        uca1 = s2u[a1]
        uca2 = s2u[a2]
        B1 = @view born_Z[:, :, uca1]
        B2 = @view born_Z[:, :, uca2]
        Dblock = @view dm[:, :, a1, a2]
        forceconstant[:, :, a1, a2] .= B1 * Dblock * transpose(B2)
    end

    # On-site correction (same formula/order as TDEP)
    dm0 = copy(forceconstant)
    for a1 in 1:na
        m0 = zeros(3, 3)
        for a2 in 1:na
            m0 .+= dm0[:, :, a1, a2]
        end
        forceconstant[:, :, a1, a1] .-= m0
    end

    return forceconstant
end



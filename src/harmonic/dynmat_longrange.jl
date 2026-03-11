# Long-range dipole-dipole dynamical matrix via Ewald summation
#
# Faithful port of TDEP's lo_longrange_electrostatics and lo_longrange_electrostatics_dynmat
# (libolle) for computing the non-analytic part of the dynamical matrix at arbitrary q.


# -----------------------------------------------------------------------------
# Geometry helpers (TDEP: geometryfunctions)
# -----------------------------------------------------------------------------

"""Radius of largest sphere inscribed in parallelepiped with columns of L as edges."""
function inscribed_sphere_in_box(L::AbstractMatrix)
    v1, v2, v3 = eachcol(L)
    A1 = norm(cross(Vector(v2), Vector(v3)))
    A2 = norm(cross(Vector(v3), Vector(v1)))
    A3 = norm(cross(Vector(v1), Vector(v2)))
    V = abs(det(L))
    min(V / (2 * A1), V / (2 * A2), V / (2 * A3))
end

"""Radius of smallest sphere containing the parallelepiped (half max vertex distance)."""
function bounding_sphere_of_box(L::AbstractMatrix)
    v1, v2, v3 = eachcol(L)
    rmax = 0.0
    for s1 in (-1, 1), s2 in (-1, 1), s3 in (-1, 1)
        v = SVector{3,Float64}(s1 * v1[1] + s2 * v2[1] + s3 * v3[1],
                              s1 * v1[2] + s2 * v2[2] + s3 * v3[2],
                              s1 * v1[3] + s2 * v2[3] + s3 * v3[3])
        rmax = max(rmax, sqrt(sqnorm(v)))
    end
    0.5 * rmax
end

# -----------------------------------------------------------------------------
# Ewald real-space dipole kernel (TDEP: ewald_H_thingy)
# -----------------------------------------------------------------------------

"""
    ewald_H_thingy(x, y, inveps) -> SMatrix{3,3,Float64}

Real-space dipole-dipole interaction kernel in Ewald summation.
`x` = λ * δ (scaled displacement in dielectric metric),
`y` = λ * D (scaled distance),
`inveps` = inverse dielectric tensor.
Returns 3×3 matrix H.
"""
function ewald_H_thingy(x::SVector{3,T}, y::T, inveps::AbstractMatrix) where {T<:AbstractFloat}
    H = @SMatrix zeros(T, 3, 3)
    if abs(y) < 1e-10
        return H
    end
    erfcy = erfc(y)
    invy2 = 1.0 / (y * y)
    invy3 = invy2 / y
    expminvy2 = exp(-y * y)
    f0 = invy2 * (3 * erfcy * invy3 + TWO_OVER_SQRTPI * expminvy2 * (3 * invy2 + 2))
    f1 = erfcy * invy3 + TWO_OVER_SQRTPI * expminvy2 * invy2
    H = SMatrix{3,3,T}([x[i] * x[j] * f0 - inveps[i, j] * f1 for i in 1:3, j in 1:3])
    return H
end

# -----------------------------------------------------------------------------
# Ewald parameters (TDEP: lo_ewald_parameters)
# -----------------------------------------------------------------------------

"""
    EwaldParameters

Container for dipole Ewald summation: λ, dielectric tensor, R and G lattice vectors.
"""
mutable struct EwaldParameters
    lambda::Float64
    eps::SMatrix{3,3,Float64,9}
    inveps::SMatrix{3,3,Float64,9}
    Rvec::Matrix{Float64}  # (3, n_R)
    Gvec::Matrix{Float64}  # (3, n_G)
end

function EwaldParameters()
    EwaldParameters(-1.0, @SMatrix(zeros(3, 3)), @SMatrix(zeros(3, 3)),
                    zeros(3, 0), zeros(3, 0))
end

"""
    _brent_root(f, a, b; atol, maxiters=1000, context="")

Bracketed root solve using Brent with explicit residual check. Keeps failure
messages contextual and consistent across Ewald parameter solves.
"""
function _brent_root(f, a::Real, b::Real; atol::Real, maxiters::Int=1000, context::AbstractString="")
    root = Roots.find_zero(
        f, (Float64(a), Float64(b)), Roots.Brent();
        atol=Float64(atol), rtol=0.0,
        maxiters=maxiters,
    )
    abs(f(root)) < atol || error("Could not converge $(context)")
    return root
end

# -----------------------------------------------------------------------------
# Ewald dipole K and R radii (TDEP: ewald_dipole_k_r)
# -----------------------------------------------------------------------------

function _ewald_dipole_k_r(lambda::Float64, pts::AbstractMatrix, eps::AbstractMatrix,
                           inveps::AbstractMatrix, L::AbstractMatrix, L_rec::AbstractMatrix,
                           tol::Float64)::Tuple{Float64,Float64}
    inv4lam2 = 1.0 / (4.0 * lambda * lambda)
    invdete = 1.0 / sqrt(det(eps))
    npts = size(pts, 2)

    # Find krad: radius where sum exp(-K²/4λ²)*K² over sphere = tol
    k_inscribed = inscribed_sphere_in_box(L_rec)
    f_k(r) = sum(i -> begin
        v = SVector{3,Float64}(pts[1,i]*r*TWOPI, pts[2,i]*r*TWOPI, pts[3,i]*r*TWOPI)
        knorm = dot(v, eps * v)
        exp(-knorm * inv4lam2) * knorm
    end, 1:npts) - tol

    k0 = k_inscribed * 1e-9
    k1 = k_inscribed * 1e9
    fk0 = f_k(k0)
    fk1 = f_k(k1)
    fk0 * fk1 > 0 && error("Could not bracket K-space radius in _ewald_dipole_k_r")
    kroot = _brent_root(
        f_k, k0, k1;
        atol=tol * 1e-3,
        maxiters=1000,
        context="K-space radius in _ewald_dipole_k_r",
    )
    krad = kroot + bounding_sphere_of_box(L_rec)

    # Find rrad: real-space decay radius
    f_r(r) = sum(i -> begin
        v = inveps * SVector{3,Float64}(pts[1,i]*r, pts[2,i]*r, pts[3,i]*r)
        rr = SVector{3,Float64}(pts[1,i]*r, pts[2,i]*r, pts[3,i]*r)
        D = sqrt(dot(v, rr))
        (lambda^3) * sum(abs, ewald_H_thingy(SVector{3,Float64}(lambda .* v), lambda * D, inveps)) * invdete
    end, 1:npts) - tol

    r_inscribed = inscribed_sphere_in_box(L)
    r0 = r_inscribed * 1e-4
    r1 = r_inscribed * 1e4
    fr0 = f_r(r0)
    fr1 = f_r(r1)
    fr0 * fr1 > 0 && error("Could not bracket R-space radius in _ewald_dipole_k_r")
    rroot = _brent_root(
        f_r, r0, r1;
        atol=tol * 1e-3,
        maxiters=1000,
        context="R-space radius in _ewald_dipole_k_r",
    )
    rrad = rroot + bounding_sphere_of_box(L)
    return krad, rrad
end

# Points on unit sphere (simple Fibonacci / random; TDEP uses lo_points_on_sphere)
function _points_on_sphere(n::Int=40)
    pts = zeros(3, n)
    for i in 1:n
        θ = acos(1 - 2 * (i - 0.5) / n)
        φ = π * (1 + sqrt(5)) * (i - 1)
        pts[1,i] = sin(θ) * cos(φ)
        pts[2,i] = sin(θ) * sin(φ)
        pts[3,i] = cos(θ)
    end
    pts
end

"""Increment cell dimensions to enlarge inscribed sphere (TDEP: increment_dimensions)."""
function _increment_dimensions(bd::AbstractVector{Int}, L::AbstractMatrix)
    m0 = [L[:,j] * (2*bd[j] + 1) for j in 1:3]
    m0 = hcat(m0...)
    r0 = inscribed_sphere_in_box(m0)
    best = 0
    rbest = 0.0
    for i in 1:3
        bnew = copy(bd)
        bnew[i] += 1
        m1 = [L[:,j] * (2*bnew[j] + 1) for j in 1:3]
        m1 = hcat(m1...)
        r1 = inscribed_sphere_in_box(m1)
        if r1 > rbest && abs(r1 - r0) > lo_tol
            rbest = r1
            best = i
        end
    end
    if best > 0
        bd2 = copy(bd)
        bd2[best] += 1
        return bd2
    end
    k = argmin(bd)
    bd2 = copy(bd)
    bd2[k] += 1
    return bd2
end

"""
    set_ewald_parameters!(ew::EwaldParameters, uc::CrystalStructure, eps::AbstractMatrix;
                         strategy::Int=1, tol::Float64=1e-20, lambda_forced::Union{Nothing,Float64}=nothing)

Fill Ewald parameters for unit cell and dielectric tensor.
- strategy 1: optimize λ for speed (balance K vs R sums)
- strategy 2: range-separation at ~2× nn distance
- strategy 3: use lambda_forced
"""
function set_ewald_parameters!(ew::EwaldParameters, uc::CrystalStructure, eps::AbstractMatrix;
                              strategy::Int=1, tol::Float64=1e-20,
                              lambda_forced::Union{Nothing,Float64}=nothing)
    L = uc.L
    # TDEP's reciprocal basis here is WITHOUT 2π; factors are applied inside kernels.
    L_rec = inv(L)'  # reciprocal lattice vectors as columns
    ew.eps = SMatrix{3,3}(eps)
    ew.inveps = SMatrix{3,3}(inv(eps))

    pts = _points_on_sphere(40)
    lambda = 0.0
    if strategy == 3
        isnothing(lambda_forced) && error("strategy=3 requires lambda_forced")
        lambda = lambda_forced
    elseif strategy == 1
        # Balance n_R and n_G: find lambda where (nR - nG)/(nR + nG) ≈ 0
        y1(x) = begin
            kr, rr = _ewald_dipole_k_r(x, pts, eps, ew.inveps, L, L_rec, tol)
            nR = 4 * (4π/3) * rr^3 / abs(det(L))
            nG = (4π/3) * kr^3 * abs(det(L))
            (nR - nG) / (nR + nG + 1e-30)
        end
        x0, x1 = 0.3, 3.0
        has_bracket = false
        for _ in 1:30
            ylo, yhi = y1(x0), y1(x1)
            if ylo * yhi < 0
                has_bracket = true
                break
            end
            ylo > 0 ? (x0 *= 0.5) : (x1 *= 2)
        end
        has_bracket || error("Could not bracket lambda for strategy=1")
        lambda = _brent_root(
            y1, x0, x1;
            atol=1e-3,
            maxiters=1000,
            context="lambda for strategy=1",
        )
    else
        # strategy 2: range separation
        nn_dist = minimum(norm(uc.x_cart[i] - uc.x_cart[j]) for i in 1:length(uc) for j in (i+1):length(uc))
        r_target = 2.0 * nn_dist
        invdete = 1.0 / sqrt(det(eps))
        pts = _points_on_sphere(40)
        l0, l1 = 10.0, 1e-8
        y2(lmid) = begin
            acc = 0.0
            for i in 1:size(pts,2)
                v = ew.inveps * (pts[:,i] * r_target)
                D = sqrt(dot(v, pts[:,i] * r_target))
                acc += (lmid^3) * sum(abs, ewald_H_thingy(SVector{3}(lmid .* v), lmid*D, ew.inveps)) * invdete
            end
            acc / size(pts,2) - 1e-5
        end
        y0, y1 = y2(l0), y2(l1)
        y0 * y1 < 0 || error("Could not bracket lambda for strategy=2")
        lambda = _brent_root(
            y2, l0, l1;
            atol=1e-11,
            maxiters=1000,
            context="lambda for strategy=2",
        )
    end

    ew.lambda = lambda

    # Get krad, rrad
    krad, rrad = _ewald_dipole_k_r(lambda, pts, eps, ew.inveps, L, L_rec, tol)

    # Build R vectors
    bd = [1, 1, 1]
    while true
        m0 = hcat([L[:,j] * (2*bd[j] + 1) for j in 1:3]...)
        if inscribed_sphere_in_box(m0) > rrad
            break
        end
        bd = _increment_dimensions(bd, L)
    end
    r2 = rrad^2
    Rlist = Vector{SVector{3,Float64}}()
    for i in -bd[1]:bd[1], j in -bd[2]:bd[2], k in -bd[3]:bd[3]
        v_frac = SVector{3,Float64}(Float64(i), Float64(j), Float64(k))
        v_cart = L * v_frac
        if sqnorm(v_cart) < r2
            push!(Rlist, v_cart)
        end
    end
    sort!(Rlist, by=v->-sqnorm(v))  # largest first (TDEP convention)
    # Avoid `hcat()` on empty (=> Matrix{Any}); ensure at least R=0 is present.
    if isempty(Rlist)
        ew.Rvec = zeros(3, 0)
    else
        ew.Rvec = hcat(Rlist...)
    end

    # Build G vectors
    bd = [1, 1, 1]
    while true
        m0 = hcat([L_rec[:,j] * (2*bd[j] + 1) for j in 1:3]...)
        if inscribed_sphere_in_box(m0) > krad
            break
        end
        bd = _increment_dimensions(bd, L_rec)
    end
    k2 = krad^2
    Glist = Vector{SVector{3,Float64}}()
    for i in -bd[1]:bd[1], j in -bd[2]:bd[2], k in -bd[3]:bd[3]
        v_frac = SVector{3,Float64}(Float64(i), Float64(j), Float64(k))
        v_cart = L_rec * v_frac
        if sqnorm(v_cart) < k2
            push!(Glist, v_cart)
        end
    end
    sort!(Glist, by=v->-sqnorm(v))
    # Avoid `hcat()` on empty (=> Matrix{Any}); return empty Matrix{Float64} instead
    if isempty(Glist)
        ew.Gvec = zeros(3, 0)
    else
        ew.Gvec = hcat(Glist...)
    end

    return ew
end

# -----------------------------------------------------------------------------
# Long-range dynamical matrix (TDEP: longrange_dynamical_matrix)
# -----------------------------------------------------------------------------

"""
    longrange_dynamical_matrix!(D, ew, uc, q_frac, born_Z, eps; reconly=false, born_onsite=nothing)

Compute dipole-dipole long-range dynamical matrix D(q) in-place (no derivatives).
For D and ∂D/∂q_cart, use `longrange_dynamical_matrix_with_gradient!`.

- `reconly`: if `true`, only compute reciprocal-space part (TDEP uses `true`).
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
    born_onsite,
    reconly::Bool,
    chgmult::Bool,
    compute_grad::Bool,
    Dx,
    Dy,
    Dz,
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
            D0[:, :, a2, a1] = adjoint(D0[:, :, a1, a2])
            if compute_grad
                Dx0[:, :, a2, a1] = adjoint(Dx0[:, :, a1, a2])
                Dy0[:, :, a2, a1] = adjoint(Dy0[:, :, a1, a2])
                Dz0[:, :, a2, a1] = adjoint(Dz0[:, :, a1, a2])
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
function longrange_dynamical_matrix(ew::EwaldParameters, uc::CrystalStructure,
    q_frac::SVector{3,Float64}, born_Z::AbstractArray, eps::AbstractMatrix;
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
function longrange_dynamical_matrix(ifc2::IFC2, uc::CrystalStructure, q_frac::SVector{3,Float64}; reconly::Bool=false)
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
    supercell_longrange_forceconstant(ew, born_Z, eps, sc::CrystalStructure, uc::CrystalStructure; thres=1e-8) -> Array{Float64,4}

Long-range contribution to force constants at Γ for supercell.
Returns real (3,3,na,na) - add to short-range IFCs for full force constants.
"""
function supercell_longrange_forceconstant(ew::EwaldParameters, born_Z::AbstractArray,
    eps::AbstractMatrix, sc::CrystalStructure, uc::CrystalStructure;
    thres::Float64=1e-8,
)
    na = length(sc)
    # Keep reciprocal basis convention consistent with TDEP (no 2π here).
    L_rec = inv(sc.L)'
    inv4lam2 = 1.0 / (4.0 * ew.lambda^2)
    na_uc = length(uc)
    size(born_Z, 3) == na_uc || error("For TDEP parity, born_Z must be unit-cell indexed: size(born_Z,3)=$(size(born_Z,3)) but length(uc)=$na_uc.")
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
            dm[:, :, a2, a1] = transpose(dm[:, :, a1, a2])
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



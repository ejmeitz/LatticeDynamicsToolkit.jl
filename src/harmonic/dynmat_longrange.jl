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

# -----------------------------------------------------------------------------
# Ewald dipole K and R radii (TDEP: ewald_dipole_k_r)
# -----------------------------------------------------------------------------

function _ewald_dipole_k_r(lambda::Float64, pts::AbstractMatrix, eps::AbstractMatrix,
                           inveps::AbstractMatrix, L::AbstractMatrix, L_rec::AbstractMatrix,
                           tol::Float64)::Tuple{Float64,Float64}
    inv4lam2 = 1.0 / (4.0 * lambda * lambda)
    invdete = 1.0 / sqrt(abs(det(eps)))
    npts = size(pts, 2)

    # Find krad: radius where sum exp(-K²/4λ²)*K² over sphere = tol (f increases with r)
    k_inscribed = inscribed_sphere_in_box(L_rec)
    f_k(r) = sum(i -> begin
        v = SVector{3,Float64}(pts[1,i]*r*TWOPI, pts[2,i]*r*TWOPI, pts[3,i]*r*TWOPI)
        knorm = dot(v, eps * v)
        exp(-knorm * inv4lam2) * knorm
    end, 1:npts) - tol

    k0 = k_inscribed * 1e-9
    k1 = k_inscribed * 0.5
    for _ in 1:20
        f_k(k1) >= tol && break
        k1 = min(k1 * 2, k_inscribed * 10)
    end
    for _ in 1:100
        kmid = 0.5 * (k0 + k1)
        abs(f_k(kmid)) < tol * 1e-2 && (k0 = kmid; break)
        f_k(k0) * f_k(kmid) < 0 ? (k1 = kmid) : (k0 = kmid)
    end
    krad = 0.5 * (k0 + k1) + bounding_sphere_of_box(L_rec)

    # Find rrad: real-space decay radius
    f_r(r) = (sum(i -> begin
        v = inveps * SVector{3,Float64}(pts[1,i]*r, pts[2,i]*r, pts[3,i]*r)
        rr = SVector{3,Float64}(pts[1,i]*r, pts[2,i]*r, pts[3,i]*r)
        D = sqrt(max(1e-20, dot(v, rr)))
        (lambda^3) * sum(abs, ewald_H_thingy(SVector{3,Float64}(lambda .* v), lambda * D, inveps)) * invdete
    end, 1:npts) / npts) - tol

    r_inscribed = inscribed_sphere_in_box(L)
    r0 = r_inscribed * 0.01
    r1 = r_inscribed * 4
    for _ in 1:15
        f_r(r1) <= 0 && break
        r1 = min(r1 * 2, r_inscribed * 100)
    end
    for _ in 1:80
        rmid = 0.5 * (r0 + r1)
        abs(f_r(rmid)) < tol * 0.1 && (r0 = rmid; break)
        f_r(r0) * f_r(rmid) < 0 ? (r1 = rmid) : (r0 = rmid)
    end
    rrad = max(0.5 * (r0 + r1), r_inscribed * 0.5) + bounding_sphere_of_box(L)
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
                         strategy::Int=1, tol::Float64=1e-8, lambda_forced::Union{Nothing,Float64}=nothing)

Fill Ewald parameters for unit cell and dielectric tensor.
- strategy 1: optimize λ for speed (balance K vs R sums)
- strategy 2: range-separation at ~2× nn distance
- strategy 3: use lambda_forced
"""
function set_ewald_parameters!(ew::EwaldParameters, uc::CrystalStructure, eps::AbstractMatrix;
                              strategy::Int=1, tol::Float64=1e-8,
                              lambda_forced::Union{Nothing,Float64}=nothing)
    L = uc.L
    L_rec = TWOPI * inv(L)'  # reciprocal lattice vectors as columns
    # Cap radii for very small cells to avoid huge loops
    k_insc = inscribed_sphere_in_box(L_rec)
    r_insc = inscribed_sphere_in_box(L)
    ew.eps = SMatrix{3,3}(eps)
    ew.inveps = SMatrix{3,3}(inv(eps))

    pts = _points_on_sphere(40)
    lambda = 0.0
    if strategy == 3
        isnothing(lambda_forced) && error("strategy=3 requires lambda_forced")
        lambda = lambda_forced
    elseif strategy == 1
        # Balance n_R and n_G: find lambda where (nR - nG)/(nR + nG) ≈ 0
        y(x) = begin
            kr, rr = _ewald_dipole_k_r(x, pts, eps, ew.inveps, L, L_rec, tol)
            nR = 4 * (4π/3) * rr^3 / abs(det(L))
            nG = (4π/3) * kr^3 * abs(det(L))
            (nR - nG) / (nR + nG + 1e-30)
        end
        x0, x1 = 0.3, 3.0
        for _ in 1:30
            y0, y1 = y(x0), y(x1)
            if y0 * y1 < 0
                break
            end
            y0 > 0 ? (x0 *= 0.5) : (x1 *= 2)
        end
        for _ in 1:80
            xmid = 0.5 * (x0 + x1)
            ymid = y(xmid)
            abs(ymid) < 1e-3 && (lambda = xmid; break)
            y(x0) * ymid < 0 ? (x1 = xmid) : (x0 = xmid)
            lambda = xmid
        end
        lambda <= 0 && (lambda = 0.5 * (x0 + x1))
    else
        # strategy 2: range separation
        nn_dist = minimum(norm(uc.x_cart[i] - uc.x_cart[j]) for i in 1:length(uc) for j in (i+1):length(uc))
        r_target = 2.0 * nn_dist
        invdete = 1.0 / sqrt(det(eps))
        pts = _points_on_sphere(40)
        l0, l1 = 10.0, 1e-8
        for _ in 1:100
            lmid = sqrt(l0 * l1)
            y = 0.0
            for i in 1:size(pts,2)
                v = ew.inveps * (pts[:,i] * r_target)
                D = sqrt(dot(v, pts[:,i] * r_target))
                y += (lmid^3) * sum(abs, ewald_H_thingy(SVector{3}(lmid .* v), lmid*D, ew.inveps)) * invdete
            end
            y = y / size(pts,2) - 1e-5
            if abs(y) < 1e-8
                lambda = lmid
                break
            end
            y > 0 ? (l0 = lmid) : (l1 = lmid)
            lambda = lmid
        end
    end

    ew.lambda = lambda

    # Get krad, rrad
    if strategy == 3
        # Heuristic radii when lambda is provided (skip expensive root-finding)
        # Cap to avoid huge R/G sets for small cells
        krad = min(2.0 / max(lambda, 0.2), k_insc * 4) + bounding_sphere_of_box(L_rec)
        rrad = min(r_insc * 2.5, r_insc * 6) + bounding_sphere_of_box(L)
    else
        krad, rrad = _ewald_dipole_k_r(lambda, pts, eps, ew.inveps, L, L_rec, tol)
    end

    # Build R vectors
    bd = [1, 1, 1]
    for _ in 1:50
        m0 = hcat([L[:,j] * (2*bd[j] + 1) for j in 1:3]...)
        if inscribed_sphere_in_box(m0) > rrad || maximum(bd) > 15
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
    for _ in 1:50
        m0 = hcat([L_rec[:,j] * (2*bd[j] + 1) for j in 1:3]...)
        if inscribed_sphere_in_box(m0) > krad || maximum(bd) > 15
            break
        end
        bd = _increment_dimensions(bd, L_rec)
    end
    k2 = krad^2
    Glist = Vector{SVector{3,Float64}}()
    for i in -bd[1]:bd[1], j in -bd[2]:bd[2], k in -bd[3]:bd[3]
        v_frac = SVector{3,Float64}(Float64(i), Float64(j), Float64(k))
        v_cart = L_rec * v_frac
        if sqnorm(v_cart) < k2 && sqnorm(v_cart) > lo_sqtol
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
    longrange_dynamical_matrix!(D, ew, uc, q_frac, born_Z, eps; reconly=false)

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
)
    _longrange_dynamical_matrix_impl!(
        D, ew, uc, q_frac, born_Z, eps;
        born_onsite=nothing,
        reconly=reconly,
        chgmult=true,
        compute_grad=false,
        Dx=nothing,
        Dy=nothing,
        Dz=nothing,
    )
end

"""
    longrange_dynamical_matrix_with_gradient!(D, Dx, Dy, Dz, ew, uc, q_frac, born_Z, eps; reconly=false)

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
)
    _longrange_dynamical_matrix_impl!(
        D, ew, uc, q_frac, born_Z, eps;
        born_onsite=nothing,
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

    q_cart = TWOPI * (uc.L_inv' * q_frac)
    x_cart = uc.x_cart
    inveps = ew.inveps
    dete = 1.0 / sqrt(det(eps))
    mult = chgmult
    onsite = something(born_onsite, zeros(3, 3, na))

    npair = na * (na + 1) ÷ 2
    ucvl = [chop3(x_cart[a2] - x_cart[a1], sqrt(lo_sqtol)) for a1 in 1:na for a2 in a1:na]

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
        K = chop3(G + q_cart, sqrt(lo_sqtol))
        if sqnorm(K) < lo_sqtol
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
                Dvec = chop3(R + ucvl[ip], sqrt(lo_sqtol))
                ikr = dot(q_cart, R)
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
        for a2 in 1:na, a1 in 1:na, j in 1:3, i in 1:3
            for jj in 1:3, ii in 1:3
                D[i, j, a1, a2] += born_Z[i, ii, a1] * born_Z[j, jj, a2] * D0[ii, jj, a1, a2]
                if compute_grad
                    Dx[i, j, a1, a2] += born_Z[i, ii, a1] * born_Z[j, jj, a2] * Dx0[ii, jj, a1, a2]
                    Dy[i, j, a1, a2] += born_Z[i, ii, a1] * born_Z[j, jj, a2] * Dy0[ii, jj, a1, a2]
                    Dz[i, j, a1, a2] += born_Z[i, ii, a1] * born_Z[j, jj, a2] * Dz0[ii, jj, a1, a2]
                end
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

    chop!(D, sqrt(lo_sqtol))
    if compute_grad
        chop!(Dx, sqrt(lo_sqtol))
        chop!(Dy, sqrt(lo_sqtol))
        chop!(Dz, sqrt(lo_sqtol))
    end
    return D
end

"""
    longrange_dynamical_matrix(ew, uc, q_frac, born_Z, eps; reconly=false) -> Array{ComplexF64,4}

Allocating version returning (3,3,na,na) complex dynamical matrix.
"""
function longrange_dynamical_matrix(ew::EwaldParameters, uc::CrystalStructure,
    q_frac::SVector{3,Float64}, born_Z::AbstractArray, eps::AbstractMatrix; reconly::Bool=false)
    na = length(uc)
    D = zeros(ComplexF64, 3, 3, na, na)
    longrange_dynamical_matrix!(D, ew, uc, q_frac, born_Z, eps; reconly=reconly)
    return D
end

"""
    longrange_dynamical_matrix(ifc2::IFC2, uc::CrystalStructure, q_frac; reconly=false)

Convenience: use polar data from IFC2. Returns nothing if `has_polar_data` is false.
"""
function longrange_dynamical_matrix(ifc2::IFC2, uc::CrystalStructure, q_frac::SVector{3,Float64}; reconly::Bool=false)
    ifc2.has_polar_data || return nothing
    ew = EwaldParameters()
    λ = ifc2.polar.lambda > 0 ? ifc2.polar.lambda : 0.5  # fallback if unset
    set_ewald_parameters!(ew, uc, ifc2.polar.eps; strategy=3, lambda_forced=λ)
    longrange_dynamical_matrix(ew, uc, q_frac, ifc2.polar.born_Z, ifc2.polar.eps; reconly=reconly)
end

# -----------------------------------------------------------------------------
# Supercell long-range force constants at Gamma (TDEP: supercell_longrange_forceconstant)
# -----------------------------------------------------------------------------

"""
    supercell_longrange_forceconstant(ew, born_Z, eps, sc::CrystalStructure; thres=1e-8) -> Array{Float64,4}

Long-range contribution to force constants at Γ for supercell.
Returns real (3,3,na,na) - add to short-range IFCs for full force constants.
"""
function supercell_longrange_forceconstant(ew::EwaldParameters, born_Z::AbstractArray,
    eps::AbstractMatrix, sc::CrystalStructure; thres::Float64=1e-8)
    na = length(sc)
    L_rec = TWOPI * inv(sc.L)'
    inv4lam2 = 1.0 / (4.0 * ew.lambda^2)

    # G-vectors in reciprocal space (Gamma: q=0)
    krad = inscribed_sphere_in_box(L_rec) * 0.5
    for iter in 1:1000
        pts = _points_on_sphere(20)
        dampsum = 0.0
        for i in 1:size(pts,2)
            v = SVector{3,Float64}(pts[1,i], pts[2,i], pts[3,i]) * krad * TWOPI
            knorm = dot(v, eps * v)
            dampsum += exp(-knorm * inv4lam2) * knorm
        end
        if dampsum <= thres
            break
        end
        krad += krad * 0.25
    end

    ndim = 0
    for i in 1:100
        ndim = i
        m0 = [L_rec[:,j] * (2*ndim + 1) for j in 1:3]
        m0 = hcat(m0...)
        if inscribed_sphere_in_box(m0) > krad + lo_tol
            break
        end
    end

    k2 = krad^2
    qvecs = SVector{3,Float64}[]
    for i in -ndim:ndim, j in -ndim:ndim, k in -ndim:ndim
        v = L_rec * SVector{3,Float64}(Float64(i), Float64(j), Float64(k))
        if sqnorm(v) < k2 && sqnorm(v) > lo_sqtol
            push!(qvecs, v)
        end
    end
    nvec = length(qvecs)

    # All pair deltas (Cartesian); build mapping for unique deltas
    pair_deltas = [chop3(sc.x_cart[a2] - sc.x_cart[a1], sqrt(lo_sqtol))
                   for a1 in 1:na for a2 in a1:na]
    # Unique deltas with tolerance
    deltavec = SVector{3,Float64}[]
    deltavecind = Int[]
    for (ip, d) in enumerate(pair_deltas)
        idx = findfirst(dv -> sqnorm(dv - d) < lo_sqtol, deltavec)
        if idx === nothing
            push!(deltavec, d)
            push!(deltavecind, length(deltavec))
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

    # Multiply in Born charges (use unit-cell indices for polar data)
    # For supercell we need index_in_unitcell; if sc is just a repetition, uca = ((a-1) % na_uc) + 1
    na_uc = size(born_Z, 3)
    forceconstant = zeros(Float64, 3, 3, na, na)
    for a2 in 1:na, a1 in 1:na
        uca1 = ((a1 - 1) % na_uc) + 1
        uca2 = ((a2 - 1) % na_uc) + 1
        m0 = zeros(3, 3)
        for j in 1:3, i in 1:3, jj in 1:3, ii in 1:3
            m0[i, j] += born_Z[i, ii, uca1] * born_Z[j, jj, uca2] * dm[ii, jj, a1, a2]
        end
        forceconstant[:, :, a1, a2] = m0
    end

    # On-site correction: Φ_ii = -Σ_{j≠i} Φ_ij (ASR-like for long-range)
    for a1 in 1:na
        s = sum(forceconstant[:, :, a1, a2] for a2 in 1:na if a2 != a1)
        forceconstant[:, :, a1, a1] -= s
    end

    return forceconstant
end



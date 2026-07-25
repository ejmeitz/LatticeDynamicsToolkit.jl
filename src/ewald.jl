# -----------------------------------------------------------------------------
# Geometry helpers (TDEP: geometryfunctions)
# -----------------------------------------------------------------------------

"""Radius of largest sphere inscribed in parallelepiped with columns of L as edges."""
function inscribed_sphere_in_box(L::AbstractMatrix)
    v1, v2, v3 = eachcol(L)
    A1 = norm(cross(v2, v3))
    A2 = norm(cross(v3, v1))
    A3 = norm(cross(v1, v2))
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
function ewald_H_thingy(
    x::SVector{3,T},
    y::T,
    inveps::AbstractMatrix,
) where {T<:AbstractFloat}
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
function _brent_root(
    f::Function,
    a::Real,
    b::Real;
    atol::Real,
    maxiters::Int=1000,
    context::AbstractString="",
)
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

function _ewald_dipole_k_r(
    lambda::Float64,
    pts::AbstractMatrix,
    eps::AbstractMatrix,
    inveps::AbstractMatrix,
    L::AbstractMatrix,
    L_rec::AbstractMatrix,
    tol::Float64,
)

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
function _increment_dimensions(
    bd::MVector{3,Int},
    L::AbstractMatrix,
)
    m0 = [L[:,j] * (2*bd[j] + 1) for j in 1:3]
    m0 = hcat(m0...)
    r0 = inscribed_sphere_in_box(m0)
    best = 0
    rbest = 0.0
    for i in 1:3
        bnew = MVector{3,Int}(bd)
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
        bd[best] += 1
        return bd
    end
    k = argmin(bd)
    bd[k] += 1
    return bd
end

"""
    set_ewald_parameters!(ew::EwaldParameters, uc::CrystalStructure, eps::AbstractMatrix;
                         strategy::Int=1, tol::Float64=1e-20, lambda_forced::Union{Nothing,Float64}=nothing)

Fill Ewald parameters for unit cell and dielectric tensor.
- strategy 1: optimize λ for speed (balance K vs R sums)
- strategy 2: range-separation at ~2× nn distance
- strategy 3: use lambda_forced
"""
function set_ewald_parameters!(
    ew::EwaldParameters,
    uc::CrystalStructure,
    eps::AbstractMatrix;
    strategy::Int=1,
    tol::Float64=1e-20,
    lambda_forced::Union{Nothing,Float64}=nothing,
)
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
    bd = MVector{3,Int}(1, 1, 1)
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
    if isempty(Rlist)
        ew.Rvec = zeros(3, 0)
    else
        ew.Rvec = hcat(Rlist...)
    end

    # Build G vectors
    bd = MVector{3,Int}(1, 1, 1)
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
    if isempty(Glist)
        ew.Gvec = zeros(3, 0)
    else
        ew.Gvec = hcat(Glist...)
    end

    return ew
end

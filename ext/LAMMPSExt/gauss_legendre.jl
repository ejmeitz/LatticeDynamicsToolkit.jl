# THIS CODE IS FROM FastGaussQuadrature.jl
# One day if package extensions can have dependencies just do that
# https://github.com/JuliaApproximation/FastGaussQuadrature.jl/blob/master/src/gausslegendre.jl
# https://github.com/JuliaApproximation/FastGaussQuadrature.jl/blob/master/src/besselroots.jl
# https://github.com/JuliaApproximation/FastGaussQuadrature.jl/blob/master/src/constants.jl



const BESSELJ0_ROOTS = @SVector [
    # runic: off
    2.4048255576957728,
    5.5200781102863106,
    8.6537279129110122,
    11.791534439014281,
    14.930917708487785,
    18.071063967910922,
    21.211636629879258,
    24.352471530749302,
    27.493479132040254,
    30.634606468431975,
    33.775820213573568,
    36.917098353664044,
    40.058425764628239,
    43.199791713176730,
    46.341188371661814,
    49.482609897397817,
    52.624051841114996,
    55.765510755019979,
    58.906983926080942,
    62.048469190227170,
    # runic: on
]


const PIESSENS_C = @SMatrix [
    # runic: off
     2.883975316228  8.263194332307 11.493871452173 14.689036505931 17.866882871378 21.034784308088
     0.767665211539  4.209200330779  4.317988625384  4.387437455306  4.435717974422  4.471319438161
    -0.086538804759 -0.164644722483 -0.130667664397 -0.109469595763 -0.094492317231 -0.083234240394
     0.020433979038  0.039764618826  0.023009510531  0.015359574754  0.011070071951  0.008388073020
    -0.006103761347 -0.011799527177 -0.004987164201 -0.002655024938 -0.001598668225 -0.001042443435
     0.002046841322  0.003893555229  0.001204453026  0.000511852711  0.000257620149  0.000144611721
    -0.000734476579 -0.001369989689 -0.000310786051 -0.000105522473 -0.000044416219 -0.000021469973
     0.000275336751  0.000503054700  0.000083834770  0.000022761626  0.000008016197  0.000003337753
    -0.000106375704 -0.000190381770 -0.000023343325 -0.000005071979 -0.000001495224 -0.000000536428
     0.000042003336  0.000073681222  0.000006655551  0.000001158094  0.000000285903  0.000000088402
    -0.000016858623 -0.000029010830 -0.000001932603 -0.000000269480 -0.000000055734 -0.000000014856
     0.000006852440  0.000011579131  0.000000569367  0.000000063657  0.000000011033  0.000000002536
    -0.000002813300 -0.000004672877 -0.000000169722 -0.000000015222 -0.000000002212 -0.000000000438
     0.000001164419  0.000001903082  0.000000051084  0.000000003677  0.000000000448  0.000000000077
    -0.000000485189 -0.000000781030 -0.000000015501 -0.000000000896 -0.000000000092 -0.000000000014
     0.000000203309  0.000000322648  0.000000004736  0.000000000220  0.000000000019  0.000000000002
    -0.000000085602 -0.000000134047 -0.000000001456 -0.000000000054 -0.000000000004               0
     0.000000036192  0.000000055969  0.000000000450  0.000000000013               0               0
    -0.000000015357 -0.000000023472 -0.000000000140 -0.000000000003               0               0
     0.000000006537  0.000000009882  0.000000000043  0.000000000001               0               0
    -0.000000002791 -0.000000004175 -0.000000000014               0               0               0
     0.000000001194  0.000000001770  0.000000000004               0               0               0
    -0.000000000512 -0.000000000752               0               0               0               0
     0.000000000220  0.000000000321               0               0               0               0
    -0.000000000095 -0.000000000137               0               0               0               0
     0.000000000041  0.000000000059               0               0               0               0
    -0.000000000018 -0.000000000025               0               0               0               0
     0.000000000008  0.000000000011               0               0               0               0
    -0.000000000003 -0.000000000005               0               0               0               0
     0.000000000001  0.000000000002               0               0               0               0
    # runic: on
]


const BESSELJ1_ON_BESSELJ0_ROOTS = @SVector [
    # runic: off
    0.2695141239419169,
    0.1157801385822037,
    0.07368635113640822,
    0.05403757319811628,
    0.04266142901724309,
    0.03524210349099610,
    0.03002107010305467,
    0.02614739149530809,
    0.02315912182469139,
    0.02078382912226786,
    # runic: on
]


const AIRY_ROOTS = @SVector [
    # runic: off
    -2.338107410459767,
    -4.08794944413097,
    -5.520559828095551,
    -6.786708090071759,
    -7.944133587120853,
    -9.022650853340981,
    -10.04017434155809,
    -11.00852430373326,
    -11.93601556323626,
    -12.828776752865757,
    -13.69148903521072,
    # runic: on
]

const CUMSUMMAT_10 = [
    # runic: off
    0                    0                   0                   0                   0                    0                   0                    0                    0                    0
    0.019080722834519    0.0496969890549313 -0.0150585059796021  0.0126377679164575 -0.0118760811432484   0.0115424841953298 -0.0113725236133433   0.0112812076497144  -0.011235316890839    0.00561063519017238
    0.000812345683614654 0.14586999854807    0.0976007154946748 -0.0146972757610091  0.00680984376276729 -0.00401953146146086 0.00271970678005437 -0.00205195604894289  0.00172405556686793 -0.000812345683614662
    0.017554012345679    0.103818185816131   0.249384588781868   0.149559082892416  -0.0321899366961563   0.0210262631520163 -0.0171075837742504   0.0153341224604243  -0.0145160806571407   0.00713734567901234
    0.00286927716087872  0.136593368810421   0.201074970443365   0.339479954969535   0.164397864607267   -0.0260484364615523  0.0127399306249393  -0.00815620454308202  0.00627037388217603 -0.00286927716087872
    0.0152149561732244   0.110297082689861   0.233440527881186   0.289200104648429   0.369910942265696    0.179464641196877  -0.0375399196961666   0.0242093528947391  -0.0200259122383839   0.00947640185146695
    0.00520833333333334  0.131083537229178   0.20995020087768    0.319047619047619   0.322836242652128    0.376052442500301   0.152380952380952   -0.024100265443764    0.0127492707559062  -0.00520833333333333
    0.0131580246959603   0.114843401005169   0.227336279387047   0.299220328493314   0.347882037265605    0.337052662041377   0.316637311034378    0.12768360784343    -0.0293025419760333   0.011533333328731
    0.00673504382217329  0.127802773462876   0.21400311568839    0.313312558886712   0.332320021608814    0.355738586947393   0.289302267356911    0.240342829317707    0.0668704675171058  -0.00673504382217329
    0.0123456790123457   0.116567456572037   0.225284323338104   0.301940035273369   0.343862505804144    0.343862505804144   0.301940035273369    0.225284323338104    0.116567456572037    0.0123456790123457
    # runic: off
]

const DIFFMAT_10 = [
    # runic: off
    -27.1666666666667    33.1634374775264   -8.54863217041303   4.0               -2.42027662546121    1.70408819104185   -1.33333333333333   1.13247433143179   -1.03109120412576   0.5;
     -8.29085936938159    4.01654328417507   5.75877048314363  -2.27431608520652   1.30540728933228   -0.898197570222574   0.694592710667722 -0.586256827714545   0.532088886237956 -0.257772801031441;
      2.13715804260326   -5.75877048314363   0.927019729872654  3.75877048314364  -1.68805925749197    1.06417777247591   -0.789861687269397  0.652703644666139  -0.586256827714545  0.283118582857949;
     -1.0                 2.27431608520652  -3.75877048314364   0.333333333333335  3.06417777247591   -1.48445439793712    1.0               -0.789861687269397   0.694592710667722 -0.333333333333333;
      0.605069156365302  -1.30540728933228   1.68805925749197  -3.06417777247591   0.0895235543024196  2.87938524157182   -1.48445439793712   1.06417777247591   -0.898197570222574  0.426022047760462;
     -0.426022047760462   0.898197570222574 -1.06417777247591   1.48445439793712  -2.87938524157182   -0.0895235543024196  3.06417777247591  -1.68805925749197    1.30540728933228  -0.605069156365302;
      0.333333333333333  -0.694592710667722  0.789861687269397 -1.0                1.48445439793712   -3.06417777247591   -0.333333333333335  3.75877048314364   -2.27431608520652   1.0;
     -0.283118582857949   0.586256827714545 -0.652703644666139  0.789861687269397 -1.06417777247591    1.68805925749197   -3.75877048314364  -0.927019729872654   5.75877048314363  -2.13715804260326;
      0.257772801031441  -0.532088886237956  0.586256827714545 -0.694592710667722  0.898197570222574  -1.30540728933228    2.27431608520652  -5.75877048314363   -4.01654328417507   8.29085936938159;
     -0.5                 1.03109120412576  -1.13247433143179   1.33333333333333  -1.70408819104185    2.42027662546121   -4.0                8.54863217041303  -33.1634374775264    27.1666666666667
    # runic: on
]

@inline function gausslegendre(n::Integer)
    # GAUSSLEGENDRE(n) COMPUTE THE GAUSS-LEGENDRE NODES AND WEIGHTS IN O(n) time.

    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    elseif n == 0
        return Float64[], Float64[]
    elseif n == 1
        return [0.0], [2.0]
    elseif n == 2
        return [-1 / sqrt(3), 1 / sqrt(3)], [1.0, 1.0]
    elseif n == 3
        return [-sqrt(3 / 5), 0.0, sqrt(3 / 5)], [5 / 9, 8 / 9, 5 / 9]
    elseif n == 4
        a = 2 / 7 * sqrt(6 / 5)
        return (
            [-sqrt(3 / 7 + a), -sqrt(3 / 7 - a), sqrt(3 / 7 - a), sqrt(3 / 7 + a)],
            [
                (18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36,
                (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36,
            ],
        )
    elseif n == 5
        b = 2 * sqrt(10 / 7)
        return (
            [
                -sqrt(5 + b) / 3, -sqrt(5 - b) / 3, 0.0,
                sqrt(5 - b) / 3, sqrt(5 + b) / 3,
            ],
            [
                (322 - 13 * sqrt(70)) / 900, (322 + 13 * sqrt(70)) / 900, 128 / 225,
                (322 + 13 * sqrt(70)) / 900, (322 - 13 * sqrt(70)) / 900,
            ],
        )
    elseif n ≤ 60
        # NEWTON'S METHOD WITH THREE-TERM RECURRENCE:
        return rec(n)
    else
        # USE ASYMPTOTIC EXPANSIONS:
        return asy(n)
    end
end

function asy(n)
    # COMPUTE GAUSS-LEGENDRE NODES AND WEIGHTS USING ASYMPTOTIC EXPANSIONS.
    # COMPLEXITY O(n).

    # Nodes and weights:
    m = (n + 1) >> 1
    a = besselZeroRoots(m)
    rmul!(a, 1 / (n + 0.5))
    x = legpts_nodes(n, a)
    w = legpts_weights(n, m, a)
    # Use symmetry to get the others:
    resize!(x, n)
    resize!(w, n)
    @inbounds for i in 1:m
        x[n + 1 - i] = x[i]
        w[n + 1 - i] = w[i]
    end
    @inbounds for i in 1:m
        x[i] = -x[i]
    end
    @inbounds isodd(n) && (x[m] = 0.0)

    return x, w
end

function legpts_nodes(n, a)
    # ASYMPTOTIC EXPANSION FOR THE GAUSS-LEGENDRE NODES.
    vn = 1 / (n + 0.5)
    m = length(a)
    nodes = cot.(a)
    vn² = vn * vn
    vn⁴ = vn² * vn²
    @inbounds if n ≤ 255
        vn⁶ = vn⁴ * vn²
        for i in 1:m
            u = nodes[i]
            u² = u^2
            ai = a[i]
            ai² = ai * ai
            ai³ = ai² * ai
            ai⁵ = ai² * ai³
            node = ai + (u - 1 / ai) / 8 * vn²
            v1 = (6 * (1 + u²) / ai + 25 / ai³ - u * muladd(31, u², 33)) / 384
            v2 = u * evalpoly(u², (2595 / 15360, 6350 / 15360, 3779 / 15360))
            v3 = (1 + u²) * (
                -muladd(31 / 1024, u², 11 / 1024) / ai +
                    u / 512 / ai² + -25 / 3072 / ai³
            )
            v4 = (v2 - 1073 / 5120 / ai⁵ + v3)
            node = muladd(v1, vn⁴, node)
            node = muladd(v4, vn⁶, node)
            nodes[i] = node
        end
    elseif n ≤ 3950
        for i in 1:m
            u = nodes[i]
            u² = u^2
            ai = a[i]
            ai² = ai * ai
            ai³ = ai² * ai
            node = ai + (u - 1 / ai) / 8 * vn²
            v1 = (6 * (1 + u²) / ai + 25 / ai³ - u * muladd(31, u², 33)) / 384
            node = muladd(v1, vn⁴, node)
            nodes[i] = node
        end
    else
        for i in 1:m
            u = nodes[i]
            ai = a[i]
            node = ai + (u - 1 / ai) / 8 * vn²
            nodes[i] = node
        end
    end
    @inbounds for jj in 1:m
        nodes[jj] = cos(nodes[jj])
    end

    return nodes
end

function legpts_weights(n, m, a)
    # ASYMPTOTIC EXPANSION FOR THE GAUSS-LEGENDRE WEIGHTS.
    vn = 1 / (n + 0.5)
    vn² = vn^2
    weights = Array{Float64}(undef, m)
    if n ≤ 850000
        @inbounds for i in 1:m
            weights[i] = cot(a[i])
        end
    end
    # Split out the part that can be vectorized by llvm
    @inbounds if n ≤ 170
        for i in 1:m
            u = weights[i]
            u² = u^2
            ai = a[i]
            ai⁻¹ = 1 / ai
            ai² = ai^2
            ai⁻² = 1 / ai²
            ua = u * ai
            W1 = muladd(ua - 1, ai⁻², 1.0) / 8
            W2 = evalpoly(
                ai⁻², (
                    evalpoly(u², (-27.0, -84.0, -56.0)),
                    muladd(-3.0, muladd(u², -2.0, 1.0), 6 * ua),
                    muladd(ua, -31.0, 81.0),
                )
            ) / 384
            W3 = evalpoly(
                ai⁻¹, (
                    evalpoly(u², (153 / 1024, 295 / 256, 187 / 96, 151 / 160)),
                    evalpoly(u², (-65 / 1024, -119 / 768, -35 / 384)) * u,
                    evalpoly(u², (5 / 512, 15 / 512, 7 / 384)),
                    muladd(u², 1 / 512, -13 / 1536) * u,
                    muladd(u², -7 / 384, + 53 / 3072),
                    3749 / 15360 * u, -1125 / 1024,
                )
            )
            weights[i] = evalpoly(vn², (1 / vn² + W1, W2, W3))
        end
    elseif n ≤ 1500
        for i in 1:m
            u = weights[i]
            u² = u^2
            ai = a[i]
            ai² = ai^2
            ai⁻² = 1 / ai²
            ua = u * ai
            W1 = muladd(ua - 1, ai⁻², 1.0) / 8
            W2 = evalpoly(
                ai⁻², (
                    evalpoly(u², (-27.0, -84.0, -56.0)),
                    muladd(-3.0, muladd(u², -2.0, 1.0), 6 * ua),
                    muladd(ua, -31.0, 81.0),
                )
            ) / 384
            weights[i] = muladd(vn², W2, 1 / vn² + W1)
        end
    elseif n ≤ 850000
        for i in 1:m
            u = weights[i]
            u² = u^2
            ai = a[i]
            ai² = ai^2
            ai⁻² = 1 / ai²
            ua = u * ai
            W1 = muladd(ua - 1, ai⁻², 1.0) / 8
            weights[i] = 1 / vn² + W1
        end
    else
        for i in 1:m
            weights[i] = 1 / vn²
        end
    end
    bJ1 = besselJ1(m)
    @inbounds for i in 1:m
        weights[i] = 2 / (bJ1[i] * (a[i] / sin(a[i])) * weights[i])
    end

    return weights
end

function rec(n)
    # COMPUTE GAUSS-LEGENDRE NODES AND WEIGHTS USING NEWTON'S METHOD.
    # THREE-TERM RECURENCE IS USED FOR EVALUATION. COMPLEXITY O(n^2).

    # Initial guess:
    x = leg_initial_guess(n)

    # Perform Newton to find zeros of Legendre polynomial:
    PP1, PP2 = similar(x), similar(x)

    # Two iterations of Newton's method
    for _ in 1:2
        innerRec!(PP1, PP2, n, x)
        newt_step!(x, PP1, PP2)
    end

    # Use symmetry to get the other Legendre nodes and weights:
    m = length(x)
    resize!(x, n)
    resize!(PP2, n)
    @inbounds for i in 1:(m - 1)
        x[n + 1 - i] = -x[i]
        PP2[n + 1 - i] = -PP2[i]
    end
    @inbounds for i in 1:n
        PP2[i] = 2 / ((1 - x[i]^2) * PP2[i]^2)
    end

    return x, PP2
end

@inline function innerRec!(myPm1, myPPm1, n, x)
    # EVALUATE LEGENDRE AND ITS DERIVATIVE USING THREE-TERM RECURRENCE RELATION.
    N = size(x, 1)
    @inbounds for j in 1:N
        xj = x[j]
        Pm2 = 1.0
        Pm1 = xj
        PPm1 = 1.0
        PPm2 = 0.0
        for k in 1:(n - 1)
            Pm2, Pm1 = Pm1, muladd((2 * k + 1) * Pm1, xj, - k * Pm2) / (k + 1)
            PPm2, PPm1 = PPm1, (
                    (2 * k + 1) * muladd(xj, PPm1, Pm2) -
                    k * PPm2
                ) / (k + 1)
        end
        myPm1[j] = Pm1
        myPPm1[j] = PPm1
    end
    return myPm1, myPPm1
end

@inline function newt_step!(x, PP1, PP2)
    # In-place iteration of Newton's method.
    return @inbounds @simd for i in eachindex(x)
        x[i] -= PP1[i] / PP2[i]
    end
end

@inline function leg_initial_guess(n)
    # Returns an approximation of the first n÷2+1 roots of the Legendre polynomial.
    #  The following is equivalent to "x0 = asy(n)[1]; x = x0[1:n ÷ 2 + 1]" but it avoids unnecessary calculations.

    m = (n ÷ 2) + 1
    a = besselZeroRoots(m)
    rmul!(a, 1 / (n + 0.5))
    x = legpts_nodes(n, a)
    rmul!(x, -1.0)
    return x
end

function approx_besselroots(ν::Real, n::Integer)
    # FIXME (related issue #22 and #80)
    return approx_besselroots(Float64(ν), n)
end

function approx_besselroots(ν::Float64, n::Integer)
    # DEVELOPERS NOTES:
    #   ν = 0 --> Full Float64 precision for n ≤ 20 (Wolfram Alpha), and very
    #     accurate approximations for n > 20 (McMahon's expansion)
    #   -1 ≤ ν ≤ 5 : ν ~= 0 --> 12 decimal figures for the 6 first zeros
    #     (Piessens's Chebyshev series approximations), and very accurate
    #     approximations for the others (McMahon's expansion)
    #   ν > 5 --> moderately accurate for the 6 first zeros and good
    #     approximations for the others (McMahon's expansion)

    # This code was originally written by L. L. Peixoto in MATLAB.
    # Later modified by A. Townsend to work in Julia

    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end

    x = zeros(n)
    if ν == 0
        for k in 1:min(n, 20)
            x[k] = BESSELJ0_ROOTS[k]
        end
        for k in (min(n, 20) + 1):n
            x[k] = McMahon(ν, k)
        end
    elseif -1 ≤ ν ≤ 5
        correctFirstFew = piessens(ν)
        for k in 1:min(n, 6)
            x[k] = correctFirstFew[k]
        end
        for k in (min(n, 6) + 1):n
            x[k] = McMahon(ν, k)
        end
    elseif 5 < ν
        for k in 1:n
            x[k] = McMahon(ν, k)
        end
    end
    return x
end

function McMahon(ν::Real, k::Integer)
    # FIXME (related issue #22 and #80)
    return McMahon(Float64(ν), k)
end

function McMahon(ν::Float64, k::Integer)
    # McMahon's expansion. This expansion gives very accurate approximation
    # for the sth zero (s ≥ 7) in the whole region ν ≥ -1, and moderate
    # approximation in other cases.
    μ = 4ν^2
    a1 = 1 / 8
    a3 = (7μ - 31) / 384
    a5 = 4 * (3779 + μ * (-982 + 83μ)) / 61440 # Evaluate via Horner's method.
    a7 = 6 * (-6277237 + μ * (1585743 + μ * (-153855 + 6949μ))) / 20643840
    a9 = 144 * (2092163573 + μ * (-512062548 + μ * (48010494 + μ * (-2479316 + 70197μ)))) / 11890851840
    a11 = 720 * (
        -8249725736393 + μ * (
            1982611456181 + μ * (
                -179289628602 + μ * (
                    8903961290 + μ * (
                        -287149133 + μ * (
                            5592657
                        )
                    )
                )
            )
        )
    ) / 10463949619200
    a13 = 576 * (
        423748443625564327 + μ * (
            -100847472093088506 + μ * (
                8929489333108377 + μ * (
                    -426353946885548 + μ * (
                        13172003634537 + μ * (
                            -291245357370 + μ * (
                                4148944183
                            )
                        )
                    )
                )
            )
        )
    ) / 13059009124761600
    b = 0.25 * (2ν + 4k - 1) * π
    # Evaluate using Horner's scheme:
    x = b - (μ - 1) * (((((((a13 / b^2 + a11) / b^2 + a9) / b^2 + a7) / b^2 + a5) / b^2 + a3) / b^2 + a1) / b)
    return x
end

function _piessens_chebyshev30(ν::Float64)
    # Piessens's Chebyshev series approximations (1984). Calculates the 6 first
    # zeros to at least 12 decimal figures in region -1 ≤ ν ≤ 5:
    pt = (ν - 2) / 3

    T1 = 1.0
    T2 = pt
    T3 = 2pt * T2 - T1
    T4 = 2pt * T3 - T2
    T5 = 2pt * T4 - T3
    T6 = 2pt * T5 - T4
    T7 = 2pt * T6 - T5
    T8 = 2pt * T7 - T6
    T9 = 2pt * T8 - T7
    T10 = 2pt * T9 - T8
    T11 = 2pt * T10 - T9
    T12 = 2pt * T11 - T10
    T13 = 2pt * T12 - T11
    T14 = 2pt * T13 - T12
    T15 = 2pt * T14 - T13
    T16 = 2pt * T15 - T14
    T17 = 2pt * T16 - T15
    T18 = 2pt * T17 - T16
    T19 = 2pt * T18 - T17
    T20 = 2pt * T19 - T18
    T21 = 2pt * T20 - T19
    T22 = 2pt * T21 - T20
    T23 = 2pt * T22 - T21
    T24 = 2pt * T23 - T22
    T25 = 2pt * T24 - T23
    T26 = 2pt * T25 - T24
    T27 = 2pt * T26 - T25
    T28 = 2pt * T27 - T26
    T29 = 2pt * T28 - T27
    T30 = 2pt * T29 - T28

    T = SVector(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30)
    return T
end

function piessens(ν::Float64)
    # Piessens's Chebyshev series approximations (1984). Calculates the 6 first
    # zeros to at least 12 decimal figures in region -1 ≤ ν ≤ 5:
    C = PIESSENS_C
    T = _piessens_chebyshev30(ν)
    y = C' * T
    return Base.setindex(y, y[1] * sqrt(ν + 1), 1)
end


function besselZeroRoots(m)
    # BESSEL0ROOTS ROOTS OF BESSELJ(0,x). USE ASYMPTOTICS.
    # Use McMahon's expansion for the remainder (NIST, 10.21.19):
    jk = Array{Float64}(undef, m)
    p = (
        1071187749376 / 315, 0.0, -401743168 / 105, 0.0, 120928 / 15,
        0.0, -124 / 3, 0.0, 1.0, 0.0,
    )
    # First 20 are precomputed:
    @inbounds for jj in 1:min(m, 20)
        jk[jj] = BESSELJ0_ROOTS[jj]
    end
    @inbounds for jj in 21:min(m, 47)
        ak = π * (jj - 0.25)
        ak82 = (0.125 / ak)^2
        jk[jj] = ak + 0.125 / ak * evalpoly(ak82, (1.0, p[7], p[5], p[3]))
    end
    @inbounds for jj in 48:min(m, 344)
        ak = π * (jj - 0.25)
        ak82 = (0.125 / ak)^2
        jk[jj] = ak + 0.125 / ak * evalpoly(ak82, (1.0, p[7], p[5]))
    end
    @inbounds for jj in 345:min(m, 13191)
        ak = π * (jj - 0.25)
        ak82 = (0.125 / ak)^2
        jk[jj] = ak + 0.125 / ak * muladd(ak82, p[7], 1.0)
    end
    @inbounds for jj in 13192:m
        ak = π * (jj - 0.25)
        jk[jj] = ak + 0.125 / ak
    end
    return jk
end

function besselJ1(m)
    # BESSELJ1 EVALUATE BESSELJ(1,x)^2 AT ROOTS OF BESSELJ(0,x).
    # USE ASYMPTOTICS. Use Taylor series of (NIST, 10.17.3) and McMahon's
    # expansion (NIST, 10.21.19):
    Jk2 = Array{Float64}(undef, m)
    c = (
        -171497088497 / 15206400, 461797 / 1152, -172913 / 8064,
        151 / 80, -7 / 24, 0.0, 2.0,
    )
    # First 10 are precomputed:
    @inbounds for jj in 1:min(m, 10)
        Jk2[jj] = BESSELJ1_ON_BESSELJ0_ROOTS[jj]
    end
    @inbounds for jj in 11:min(m, 15)
        ak = π * (jj - 0.25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(evalpoly(ak2, (c[5], c[4], c[3], c[2], c[1])), ak2^2, c[7])
    end
    @inbounds for jj in 16:min(m, 21)
        ak = π * (jj - 0.25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(evalpoly(ak2, (c[5], c[4], c[3], c[2])), ak2^2, c[7])
    end
    @inbounds for jj in 22:min(m, 55)
        ak = π * (jj - 0.25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(evalpoly(ak2, (c[5], c[4], c[3])), ak2^2, c[7])
    end
    @inbounds for jj in 56:min(m, 279)
        ak = π * (jj - 0.25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(muladd(ak2, c[4], c[5]), ak2^2, c[7])
    end
    @inbounds for jj in 280:min(m, 2279)
        ak = π * (jj - 0.25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(ak2^2, c[5], c[7])
    end
    @inbounds for jj in 2280:m
        ak = π * (jj - 0.25)
        Jk2[jj] = 1 / (π * ak) * c[7]
    end
    return Jk2
end

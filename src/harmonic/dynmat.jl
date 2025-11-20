export 
    dynmat_gamma,
    get_modes

#! probably some speed to be gained by using SMatrix for dynmats

q_cart_from_frac(cell::CrystalStructure, q_frac::SVector{3,Float64}) = 2pi .* (cell.L_inv' * q_frac)


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
    D_q = zeros(ComplexF64, nb, nb)

    return dynmat_q!(D_q, fc_uc, uc, q_frac)
end

function dynmat_q!(
        D_q::Matrix{ComplexF64},
        fc_uc::IFC2,
        uc::CrystalStructure,
        q_frac::SVector{3,Float64};
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
    
    check_hermetian(D_q, nb)

    return D_q
end

function dynmat_and_derivative_q(
        fc_uc::IFC2,
        uc::CrystalStructure,
        q_frac::SVector{3,Float64};
    )

    na = length(uc)
    nb = 3*na

    D_q = zeros(ComplexF64, nb, nb)
    ∂D∂q = zeros(ComplexF64, nb, nb, 3)

    return dynmat_and_derivative_q!(D_q, ∂D∂q, fc_uc, uc, q_frac)
end

function dynmat_and_derivative_q!(
        D_q::Matrix{ComplexF64},
        ∂D∂q::Array{ComplexF64, 3},
        fc_uc::IFC2,
        uc::CrystalStructure,
        q_frac::SVector{3,Float64};
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

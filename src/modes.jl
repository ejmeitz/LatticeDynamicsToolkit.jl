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
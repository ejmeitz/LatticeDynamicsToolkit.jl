export free_energy_corrections

function free_energy_corrections(
        T::Float64,
        uc::CrystalStructure,
        ifc2::IFC2,
        ifc3::IFC3,
        ifc4::IFC4;
        mesh = [25, 25, 25],
        quantum::Bool = false,
        n_threads::Integer = Threads.nthreads()
    )

    ifc_nas = [ifc2.na, ifc3.na, ifc4.na]
    if !all(length(uc) .== ifc_nas)
        error(ArgumentError("Not all IFCs are built from the unitcell you passed which has $(length(uc)) atoms."))
    end

    k_mesh_frac = FFTMesh(uc, mesh)

    pd = PhononDispersions(uc, ifc2, k_mesh_frac; n_threads = n_threads)

    # Precompute planck derivatives on IBZ and full mesh
    p_ibz, p_all = precompute_planck_data(k_mesh_frac, pd, T, n_threads)

    # First order correction is <V4>
    F4, S4, Cv4 = free_energy_fourthorder(
                    uc,
                    ifc4,
                    k_mesh_frac,
                    pd,
                    p_ibz,
                    p_all,
                    T,
                    quantum;
                    n_threads = n_threads
                )

    # # Second order correction is <V3*V3>
    F3, S3, Cv3 = free_energy_thirdorder(
                    uc,
                    ifc3,
                    k_mesh_frac,
                    pd,
                    p_ibz,
                    p_all,
                    T,
                    quantum;
                    n_threads = n_threads
                )

    # F = U - TS
    U3 = (F3 + (T*S3)) * Hartree_to_eV
    U4 = (F4 + (T*S4)) * Hartree_to_eV

    F3 *= Hartree_to_eV
    F4 *= Hartree_to_eV

    S3 /= kB_Hartree
    S4 /= kB_Hartree

    Cv3 /= kB_Hartree
    Cv4 /= kB_Hartree

    return (F3 = F3, F4 = F4, S3 = S3, S4 = S4, U3 = U3, U4 = U4, Cv3 = Cv3, Cv4 = Cv4)

end

struct PlanckData
    n::Array{Float64}
    dn::Array{Float64}
    ddn::Array{Float64}
end

function precompute_planck_data(qp::FFTMesh, pd::PhononDispersions, T, n_threads)

    n_mode = pd.n_mode
    n_irr = n_irr_point(qp)
    n_full = n_full_q_point(qp)

    n_iq = zeros(n_irr, n_mode)
    dn_iq = zeros(n_irr, n_mode)
    ddn_iq = zeros(n_irr, n_mode)
    @tasks for q in 1:n_irr
        @set ntasks=n_threads
        for b in 1:n_mode
            ω = pd.iq[q].omega[b]
            n_iq[q, b] = planck(T, ω)
            dn_iq[q, b] = planck_deriv(T, ω)
            ddn_iq[q, b] = planck_secondderiv(T, ω)
        end
    end

    n_aq = zeros(n_full, n_mode)
    dn_aq = zeros(n_full, n_mode)
    ddn_aq = zeros(n_full, n_mode)
    @tasks for q in 1:n_full
        @set ntasks=n_threads
        for b in 1:n_mode
            ω = pd.aq[q].omega[b]
            n_aq[q, b] = planck(T, ω)
            dn_aq[q, b] = planck_deriv(T, ω)
            ddn_aq[q, b] = planck_secondderiv(T, ω)
        end
    end

    return PlanckData(n_iq, dn_iq, ddn_iq), PlanckData(n_aq, dn_aq, ddn_aq)
end

"""
    free_energy_thirdorder(uc, fct, qp, pd, temperature, quantum)

Calculate the second order contribution to free energy, entropy, and heat capacity.

# Arguments
- `uc::CrystalStructure`: crystal structure
- `fct::IFC3`: third order force constants
- `qp::FFTMesh`: q-point mesh (must be FFT mesh)
- `pd::PhononDispersions`: phonon dispersions
- `temperature::Float64`: temperature in Kelvin
- `quantum::Bool`: use quantum statistics (true) or classical (false)

# Returns
- `fe3`: third order free energy
- `s3`: third order entropy  
- `cv3`: third order heat capacity

"""
function free_energy_thirdorder(
    uc::CrystalStructure,
    fct::IFC3,
    qp::FractionalFFTMesh,
    pd::PhononDispersions,
    p_i::PlanckData,
    p_a::PlanckData,
    temperature::Float64,
    quantum::Bool;
    n_threads::Integer = Threads.nthreads()
)
    
    n_mode = pd.n_mode
    n_irr = n_irr_point(qp)
    n_full = n_full_q_point(qp)
    dims = Tuple(qp.mesh)
       
    # Compute broadening parameter for each mode
    sigsq = zeros(Float64, n_irr, n_mode)
    for q1 in 1:n_irr
        for b1 in 1:n_mode
            sigsq[q1, b1] = @views _adaptive_sigma(qp.radius, pd.iq[q1].vel[:, b1], pd.default_smearing[b1], 1.0)^2
        end
    end
    
    # Ω = (2pi)^3 / volume(uc)
    one_im = ComplexF64(1.0, 0.0)
    
    # Main loop over irreducible and full q-points
    p = Progress(n_irr, "Second Order Correction")
    (f3, s3, cv3) = @tasks for q1 in 1:n_irr

        @set ntasks = n_threads
        @set reducer = +

        @local begin
            ptf = zeros(ComplexF64, n_mode^3)
            evp1 = zeros(ComplexF64, n_mode^2)
            evp2 = zeros(ComplexF64, n_mode^3)
            egv1 = zeros(ComplexF64, n_mode)
            egv2 = zeros(ComplexF64, n_mode)
            egv3 = zeros(ComplexF64, n_mode)
            q2_full_cart = @MVector zeros(Float64, 3)
            q3_full_cart = @MVector zeros(Float64, 3)
        end

        res = [0.0, 0.0, 0.0] # f3, s3, cv3

        for q2 in 1:n_full
            q2_full_cart .= q_cart_from_frac(uc, qp.k_full[q2].r)

            # Find q3 such that q1 + q2 + q3 = 0
            q3 = fft_third_grid_index(qp.full_index_ibz[q1], q2, dims)
            
            # Skip if q3 < q2 (permutation symmetry)
            q3 < q2 && continue

            # Convert q-vector to cartesian
            q3_full_cart .= q_cart_from_frac(uc, qp.k_full[q3].r)


            # Multiplicity factor
            mult = (q2 == q3) ? 1.0 : 2.0
            
            # Prefactor including integration weights
            prefactor = qp.weights_ibz[q1] * qp.k_full[q2].weight * mult / length(uc)
            
            # Pre-transform the matrix element
            pretransform_phi3!(
                ptf,
                fct, 
                q2_full_cart,
                q3_full_cart,
                uc
            )
            
            for b1 in 1:n_mode
                # Get first phonon
                ω1 = pd.iq[q1].omega[b1]
                ω1 < lo_freqtol && continue

                n1 = p_i.n[q1, b1]
                dn1 = p_i.dn[q1, b1]
                ddn1 = p_i.ddn[q1, b1]
                egv1 .= view(pd.iq[q1].egv, :, b1) ./ sqrt(ω1)
                
                for b2 in 1:n_mode
                    # Get second phonon
                    ω2 = pd.aq[q2].omega[b2]
                    ω2 < lo_freqtol && continue

                    n2 = p_a.n[q2, b2]
                    dn2 = p_a.dn[q2, b2]
                    ddn2 = p_a.ddn[q2, b2]
                    egv2 .= view(pd.aq[q2].egv, :, b2) ./ sqrt(ω2)
                    
                    # Outer product: evp1 = egv2 * egv1^T
                    # zgeru(n_mode, n_mode, 1.0, egv2, 1, egv1, 1, evp1, n_mode)
                    fill!(evp1, 0.0)
                    evp1_mat = reshape(evp1, n_mode, n_mode)
                    BLAS.geru!(one_im, egv2, egv1, evp1_mat)
                    
                    for b3 in 1:n_mode
                        # Get third phonon
                        ω3 = pd.aq[q3].omega[b3]
                        ω3 < lo_freqtol && continue

                        n3 = p_a.n[q3, b3]
                        dn3 = p_a.dn[q3, b3]
                        ddn3 = p_a.ddn[q3, b3]
                        egv3 .= view(pd.aq[q3].egv, :, b3) ./ sqrt(ω3)
                        
                        # Project on third phonon: evp2 = egv3 * evp1^T
                        # zgeru(n_mode, n_mode^2, 1.0, egv3, 1, evp1, 1, evp2, n_mode)
                        fill!(evp2, 0.0)
                        evp2_mat = reshape(evp2, n_mode, n_mode^2)
                        BLAS.geru!(one_im, egv3, evp1, evp2_mat)
                        evp2 .= conj.(evp2)
                        
                        # Compute the scattering matrix element
                        c0 = dot(evp2, ptf)
                        psisq = abs2(c0) * prefactor
                        
                        if quantum
                            # Get smearing parameter
                            sig1 = sigsq[q1, b1]
                            sig2 = sigsq[qp.k_full[q2].irreducible_index, b2]
                            sig3 = sigsq[qp.k_full[q3].irreducible_index, b3]
                            sigma = sqrt(sig1 + sig2 + sig3)
                            
                            # Free energy formulas
                            f1 = (n1 + 1.0) * (n2 + n3 + 1.0) + n2 * n3
                            f2 = n1 * n2 + n1 * n3 - n2 * n3 + 1.0
                            
                            # Entropy formulas
                            df1 = dn1 * (n2 + n3 + 1.0) + (n1 + 1.0) * (dn2 + dn3) + dn2 * n3 + n2 * dn3
                            df2 = dn1 * n2 + n1 * dn2 + dn1 * n3 + n1 * dn3 - dn2 * n3 - n2 * dn3
                            
                            # Heat capacity formulas
                            ddf1 = ddn1 * (n2 + n3 + 1.0) + 2.0 * dn1 * (dn2 + dn3) + 
                                   (n1 + 1.0) * (ddn2 + ddn3) + ddn2 * n3 + n2 * ddn3 + 2.0 * dn2 * dn3
                            ddf2 = ddn1 * n2 + n1 * ddn2 + 2.0 * dn1 * dn2 + ddn1 * n3 + n1 * ddn3 + 
                                   2.0 * dn1 * dn3 - ddn2 * n3 - n2 * ddn3 - 2.0 * dn2 * dn3
                            
                            # Compute principal values with smearing
                            f1 = f1 * real(1.0 / (ω1 + ω2 + ω3 + im * sigma))
                            df1 = df1 * real(1.0 / (ω1 + ω2 + ω3 + im * sigma))
                            ddf1 = ddf1 * real(1.0 / (ω1 + ω2 + ω3 + im * sigma))
                            f2 = 3.0 * f2 * real(1.0 / (ω1 + ω2 - ω3 + im * sigma))
                            df2 = 3.0 * df2 * real(1.0 / (ω1 + ω2 - ω3 + im * sigma))
                            ddf2 = 3.0 * ddf2 * real(1.0 / (ω1 + ω2 - ω3 + im * sigma))
                            
                            # Combine
                            f0 = (f1 + f2) / 48.0
                            df0 = (df1 + df2) / 48.0
                            ddf0 = (ddf1 + ddf2) * temperature / 48.0
                        else
                            # Classical case
                            f0 = (kB_Hartree * temperature)^2 / (ω1 * ω2 * ω3) / 12.0
                            df0 = kB_Hartree^2 * temperature / (ω1 * ω2 * ω3) / 6.0
                            ddf0 = kB_Hartree^2 * temperature / (ω1 * ω2 * ω3) / 6.0
                        end
                        
                        # Accumulate
                        res[1] -= f0 * psisq
                        res[2] += df0 * psisq
                        res[3] += ddf0 * psisq
                    end
                end
            end
        end
        next!(p)
        res # the thing reduced
    end
    finish!(p)

    return f3, s3, cv3
end

"""
    free_energy_fourthorder(uc, fcf, qp, pd, temperature, quantum)

Calculate the first order correction to free energy, entropy, and heat capacity.
This is the first type (diagonal in q-space).

# Arguments
- `uc::CrystalStructure`: crystal structure
- `fcf::IFC4`: fourth order force constants
- `qp::FFTMesh`: q-point mesh (must be FFT mesh)
- `pd::PhononDispersions`: phonon dispersions
- `temperature::Float64`: temperature in Kelvin
- `quantum::Bool`: use quantum statistics (true) or classical (false)


# Returns
- `df4`: fourth order free energy
- `s4`: fourth order entropy
- `cv4`: fourth order heat capacity
"""
function free_energy_fourthorder(
    uc::CrystalStructure,
    fcf::IFC4,
    qp::FractionalFFTMesh,
    pd::PhononDispersions,
    p_i::PlanckData,
    p_a::PlanckData,
    temperature::Float64,
    quantum::Bool;
    n_threads::Integer = Threads.nthreads()
)
    
    n_mode = pd.n_mode
    n_irr = n_irr_point(qp)
    n_full = n_full_q_point(qp)   

    # Ω = (2pi)^3 / volume(uc)
    one_im = ComplexF64(1.0, 0.0)
    
    # Main loop
    p = Progress(n_irr, desc = "First Order Correction")
    (f4, s4, cv4) = @tasks for q1 in 1:n_irr
        @set ntasks = n_threads
        @set reducer = .+

        @local begin
            ptf = zeros(ComplexF64, n_mode^4)
            evp1 = zeros(ComplexF64, n_mode^2)
            evp2 = zeros(ComplexF64, n_mode^2)
            evp3 = zeros(ComplexF64, n_mode^4)
            egv1 = zeros(ComplexF64, n_mode)
            egv2 = zeros(ComplexF64, n_mode)
            q1_ibz_cart = @MVector zeros(Float64, 3)
            q2_full_cart = @MVector zeros(Float64, 3)
        end

        evp1_mat = reshape(evp1, n_mode, n_mode)
        evp2_mat = reshape(evp2, n_mode, n_mode)
        evp3_mat = reshape(evp3, n_mode^2, n_mode^2)

        # res = [0.0, 0.0, 0.0] # f4, s4, cv4
        f4_ = 0.0; s4_ = 0.0; cv4_ = 0.0

        q1_ibz_cart .= q_cart_from_frac(uc, qp.k_ibz[q1])

        for q2 in 1:n_full

            q2_full_cart .= q_cart_from_frac(uc, qp.k_full[q2].r)

            # Prefactor
            prefactor = qp.weights_ibz[q1] * qp.k_full[q2].weight / length(uc)
            
            pretransform_phi4!(
                ptf,
                fcf, 
                q1_ibz_cart,
                q2_full_cart,
                uc
            )
        
            for b1 in 1:n_mode
                ω1 = pd.iq[q1].omega[b1]
                ω1 < lo_freqtol && continue
                
                egv1 .= view(pd.iq[q1].egv, :, b1) ./ sqrt(ω1)

                # Outer products matching BLAS calls
                # zgerc(n_mode, n_mode, 1.0, egv1, 1, egv1, 1, evp1, n_mode)
                # evp1 = egv1 * conj(egv1)^T
                fill!(evp1, 0.0)
                BLAS.ger!(one_im, egv1, egv1, evp1_mat)
                
                for b2 in 1:n_mode
                    ω2 = pd.aq[q2].omega[b2]
                    ω2 < lo_freqtol && continue
                    
                    egv2 .= view(pd.aq[q2].egv, :, b2) ./ sqrt(ω2)
                    
                    # zgerc(n_mode, n_mode, 1.0, egv2, 1, egv2, 1, evp2, n_mode)
                    # evp2 = egv2 * conj(egv2)^T
                    fill!(evp2, 0.0)
                    BLAS.ger!(one_im, egv2, egv2, evp2_mat)
                    
                    # zgeru(n_mode^2, n_mode^2, 1.0, evp2, 1, evp1, 1, evp3, n_mode^2)
                    fill!(evp3, 0.0)
                    BLAS.geru!(one_im, evp2, evp1, evp3_mat)
                    evp3 .= conj.(evp3)
                    
                    psisq = real(dot(evp3, ptf))
                    
                    if quantum
                        n1 = p_i.n[q1, b1]
                        dn1 = p_i.dn[q1, b1]
                        ddn1 = p_i.ddn[q1, b1]
                        n2 = p_a.n[q2, b2]
                        dn2 = p_a.dn[q2, b2]
                        ddn2 = p_a.ddn[q2, b2]

                        # Free energy
                        #! SHOULD THIS PSISQ BE HERE??? ALREADY MULTIPLIED LATER
                        f0 = (2.0 * n1 + 1.0) * (2.0 * n2 + 1.0) * psisq * prefactor / 32.0
                        # Entropy
                        df0 = 2.0 * dn1 * (2.0 * n2 + 1.0) + (2.0 * n1 + 1.0) * 2.0 * dn2
                        df0 = df0 / 32.0
                        # Heat capacity
                        ddf0 = 2.0 * ddn1 * (2.0 * n2 + 1.0) + 2.0 * ddn2 * (2.0 * n1 + 1.0) + 8.0 * dn1 * dn2
                        ddf0 = ddf0 * temperature / 32.0
                    else
                        # Classical case
                        f0 = (kB_Hartree * temperature)^2 / (ω1 * ω2) / 8.0
                        df0 = kB_Hartree^2 * temperature / (ω1 * ω2) / 4.0
                        ddf0 = kB_Hartree^2 * temperature / (ω1 * ω2) / 4.0
                    end
                    
                    # Accumulate
                    f4_ += f0 * psisq * prefactor
                    s4_ -= df0 * psisq * prefactor
                    cv4_ -= ddf0 * psisq * prefactor
                end
            end
        end
        next!(p)
        (f4_, s4_, cv4_) # the thing to be reduced 
    end
    finish!(p)
    
    return f4, s4, cv4
end

"""
    pretransform_phi3(fct::IFC3, q2, q3)

Get the Fourier transform of the third order matrix element.
Returns flattened, pretransformed matrix element.
"""
function pretransform_phi3!(
        ptf::Vector{ComplexF64},
        fct::IFC3,
        q2_cart::S,
        q3_cart::S,
        uc::CrystalStructure
    ) where {S <: StaticVector{3, Float64}}

    nb = fct.na * 3
    fill!(ptf, zero(ComplexF64))
    
    for a1 in 1:fct.na
        triplets = get_interactions(fct, a1)
        m1 = uc.invsqrtm[a1]
        for trip in triplets
            a2 = trip.idxs[2]
            a3 = trip.idxs[3]

            m_factor = m1 * uc.invsqrtm[a2] * uc.invsqrtm[a3]
            
            rv2 = trip.rv2
            rv3 = trip.rv3
            
            iqr = -dot(q2_cart, rv2) - dot(q3_cart, rv3)
            expiqr = cis(iqr)  # exp(i*iqr)
            
            for i in 1:3, j in 1:3, k in 1:3
                ia = (a1 - 1) * 3 + i
                ib = (a2 - 1) * 3 + j
                ic = (a3 - 1) * 3 + k
                # Flattening scheme consistent with zgeru operations
                l = (ia - 1) * nb * nb + (ib - 1) * nb + ic
                ptf[l] += trip.ifcs[i, j, k] * expiqr * m_factor
            end
        end
    end
    
    return ptf
end

"""
    pretransform_phi4(fcf::IFC4, q1, q2)

Get the Fourier transform of the fourth order matrix element (first type).
Returns flattened, pretransformed matrix element.

This is a faithful port of the Fortran subroutine from fourthorder.f90.
"""
function pretransform_phi4!(
        ptf::Vector{ComplexF64},
        fcf::IFC4,
        q1_cart::S,
        q2_cart::S,
        uc::CrystalStructure
    ) where {S <: StaticVector{3, Float64}}
    
    nb = fcf.na * 3
    fill!(ptf, zero(ComplexF64))
    
    for a1 in 1:fcf.na
        m1 = uc.invsqrtm[a1]
        for quartet in get_interactions(fcf, a1)
            a2 = quartet.idxs[2]
            a3 = quartet.idxs[3]
            a4 = quartet.idxs[4]

            #! SHOULD PROBABLY PRE-MASS WEIGHT THE IFCs
            m_factor = m1 * uc.invsqrtm[a2] * uc.invsqrtm[a3] * uc.invsqrtm[a4]
            
            rv2 = quartet.rv2
            rv3 = quartet.rv3
            rv4 = quartet.rv4
            
            iqr = -dot(q1_cart, rv2) + dot(q2_cart, rv3) - dot(q2_cart, rv4)
            expiqr = cis(iqr)
            
            for l in 1:3, k in 1:3, j in 1:3, i in 1:3
                ia = (a1 - 1) * 3 + i
                ib = (a2 - 1) * 3 + j
                ic = (a3 - 1) * 3 + k
                id = (a4 - 1) * 3 + l
                # Flattening scheme consistent with zgeru operations
                m = (ia - 1) * nb^3 + (ib - 1) * nb^2 + (ic - 1) * nb + id
                ptf[m] += quartet.ifcs[i, j, k, l] * expiqr * (m_factor)
            end
        end
    end
    
    return ptf
end

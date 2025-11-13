

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

    k_mesh = FFTMesh(uc, mesh)
    N_prim = k_mesh.n_atoms_prim

    pd = PhononDispersions(uc, ifc2, k_mesh; n_threads = n_threads)

    F3, S3, Cv3 = free_energy_thirdorder(uc, ifc3, k_mesh, pd, T, quantum)
    F4, S4, Cv4 = free_energy_fourthorder(uc, ifc4, k_mesh, pd, T, quantum)

    # F = U - TS
    U3 = (F3 + (T*S3)) * Hartree_to_eV #! HOW TO NORMALIZE PER ATOM?
    U3 = (F4 + (T*S4)) * Hartree_to_eV

    F3 *= Hartree_to_eV
    F4 *= Hartree_to_eV

    S3 /= kB_Hartree
    S4 /= kB_Hartree

    Cv3 /= kB_Hartree
    Cv4 /= kB_Hartree

    #! TODO convert to eV / atom or per kB
    return (F3 = F3, F4 = F4, S3 = S3, S4 = S4, U3 = U3, U4 = U4, Cv3 = Cv3, Cv4 = Cv4)

end

"""
    free_energy_thirdorder(uc, fct, qp, pd, temperature, quantum)

Calculate the third order free energy, entropy, and heat capacity.

Faithful port of free_energy_thirdorder from thirdorder.f90.

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
    qp::FFTMesh,
    pd::PhononDispersions,
    temperature::Float64,
    quantum::Bool
)
    
    n_mode = pd.n_mode
    n_irr = n_irr_point(qp)
    n_full = n_full_q_point(qp)
    dims = Tuple(qp.mesh)
    
    # Allocate work arrays
    ptf = zeros(ComplexF64, n_mode^3)
    evp1 = zeros(ComplexF64, n_mode^2)
    evp2 = zeros(ComplexF64, n_mode^3)
    egv1 = zeros(ComplexF64, n_mode)
    egv2 = zeros(ComplexF64, n_mode)
    egv3 = zeros(ComplexF64, n_mode)
    
    # Compute broadening parameter for each mode
    sigsq = zeros(Float64, n_irr, n_mode)
    for q1 in 1:n_irr
        for b1 in 1:n_mode
            sigsq[q1, b1] = @views _adaptive_sigma(qp.radius, pd.iq[q1].vel[:, b1], pd.default_smearing[b1], 1.0)^2
        end
    end
    
    fe3 = 0.0
    s3 = 0.0
    cv3 = 0.0
    
    # Main loop over irreducible and full q-points
    for q1 in 1:n_irr
        for q2 in 1:n_full
            # Find q3 such that q1 + q2 + q3 = 0
            q3 = fft_third_grid_index(qp.full_index_ibz[q1], q2, dims)
            
            # Skip if q3 < q2 (permutation symmetry)
            if q3 < q2
                continue
            end
            
            # Multiplicity factor
            mult = (q2 == q3) ? 1.0 : 2.0
            
            # Prefactor including integration weights
            prefactor = qp.weights_ibz[q1] * qp.k_full[q2].weight * mult / qp.n_atoms_prim
            
            # Pre-transform the matrix element
            ptf .= pretransform_phi3(fct, qp.k_full[q2].r, qp.k_full[q3].r)
            
            for b1 in 1:n_mode
                # Get first phonon
                ω1 = pd.iq[q1].omega[b1]
                if ω1 < lo_freqtol
                    continue
                end
                
                n1 = planck(temperature, ω1)
                dn1 = planck_deriv(temperature, ω1)
                ddn1 = planck_secondderiv(temperature, ω1)
                egv1 .= view(pd.iq[q1].egv, :, b1) ./ sqrt(ω1)
                
                for b2 in 1:n_mode
                    # Get second phonon
                    ω2 = pd.aq[q2].omega[b2]
                    if ω2 < lo_freqtol
                        continue
                    end
                    
                    n2 = planck(temperature, ω2)
                    dn2 = planck_deriv(temperature, ω2)
                    ddn2 = planck_secondderiv(temperature, ω2)
                    egv2 .= view(pd.aq[q2].egv, :, b2) ./ sqrt(ω2)
                    
                    # Outer product: evp1 = egv2 * egv1^T
                    # zgeru(n_mode, n_mode, 1.0, egv2, 1, egv1, 1, evp1, n_mode)
                    # Flattened in column-major order
                    evp1 .= vec(egv2 * transpose(egv1))
                    
                    for b3 in 1:n_mode
                        # Get third phonon
                        ω3 = pd.aq[q3].omega[b3]
                        if ω3 < lo_freqtol
                            continue
                        end
                        
                        n3 = planck(temperature, ω3)
                        dn3 = planck_deriv(temperature, ω3)
                        ddn3 = planck_secondderiv(temperature, ω3)
                        egv3 .= view(pd.aq[q3].egv, :, b3) ./ sqrt(ω3)
                        
                        # Project on third phonon: evp2 = egv3 * evp1^T
                        # zgeru(n_mode, n_mode^2, 1.0, egv3, 1, evp1, 1, evp2, n_mode)
                        # Treat evp1 as a row vector of length n_mode^2
                        evp2 .= vec(egv3 * transpose(evp1))
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
                        fe3 -= f0 * psisq
                        s3 += df0 * psisq
                        cv3 += ddf0 * psisq
                    end
                end
            end
        end
    end
    
    return fe3, s3, cv3
end

"""
    free_energy_fourthorder(uc, fcf, qp, pd, temperature, quantum)

Calculate the fourth order free energy, entropy, and heat capacity.
This is the first type (diagonal in q-space).

Faithful port of free_energy_fourthorder from fourthorder.f90.

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
    qp::FFTMesh,
    pd::PhononDispersions,
    temperature::Float64,
    quantum::Bool
)
    
    n_mode = pd.n_mode
    n_irr = n_irr_point(qp)
    n_full = n_full_q_point(qp)
    
    # Allocate work arrays
    ptf = zeros(ComplexF64, n_mode^4)
    evp1 = zeros(ComplexF64, n_mode^2)
    evp2 = zeros(ComplexF64, n_mode^3)
    evp3 = zeros(ComplexF64, n_mode^4)
    egv1 = zeros(ComplexF64, n_mode)
    egv2 = zeros(ComplexF64, n_mode)
    
    df4 = 0.0
    s4 = 0.0
    cv4 = 0.0
    
    # Main loop
    for q1 in 1:n_irr
        for q2 in 1:n_full
            # Prefactor
            prefactor = qp.weights_ibz[q1] * qp.k_full[q2].weight / qp.n_atoms_prim
            
            # Pre-transform the matrix element
            ptf .= pretransform_phi4(fcf, qp.k_ibz[q1], qp.k_full[q2].r)
            
            for b1 in 1:n_mode
                ω1 = pd.iq[q1].omega[b1]
                if ω1 < lo_freqtol
                    continue
                end
                egv1 .= view(pd.iq[q1].egv, :, b1) ./ sqrt(ω1)
                
                for b2 in 1:n_mode
                    ω2 = pd.aq[q2].omega[b2]
                    if ω2 < lo_freqtol
                        continue
                    end
                    egv2 .= view(pd.aq[q2].egv, :, b2) ./ sqrt(ω2)
                    
                    # Outer products matching BLAS calls
                    # zgerc(n_mode, n_mode, 1.0, egv1, 1, egv1, 1, evp1, n_mode)
                    # evp1 = egv1 * conj(egv1)^T
                    evp1 .= vec(egv1 * adjoint(egv1))
                    
                    # zgerc(n_mode, n_mode, 1.0, egv2, 1, egv2, 1, evp2_mat, n_mode)
                    # evp2 = egv2 * conj(egv2)^T
                    evp2_mat = vec(egv2 * adjoint(egv2))
                    
                    # zgeru(n_mode^2, n_mode^2, 1.0, evp2, 1, evp1, 1, evp3, n_mode^2)
                    # Treat evp2 and evp1 as vectors and do outer product
                    evp3 .= vec(evp2_mat * transpose(evp1))
                    evp3 .= conj.(evp3)
                    
                    psisq = real(dot(evp3, ptf))
                    
                    if quantum
                        # Phonon occupation
                        n1 = planck(temperature, ω1)
                        n2 = planck(temperature, ω2)
                        # First derivative with respect to temperature
                        dn1 = planck_deriv(temperature, ω1)
                        dn2 = planck_deriv(temperature, ω2)
                        # Second derivative with respect to temperature
                        ddn1 = planck_secondderiv(temperature, ω1)
                        ddn2 = planck_secondderiv(temperature, ω2)
                        
                        # Free energy
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
                    df4 += f0 * psisq * prefactor
                    s4 -= df0 * psisq * prefactor
                    cv4 -= ddf0 * psisq * prefactor
                end
            end
        end
    end
    
    return df4, s4, cv4
end

"""
    pretransform_phi3(fct::IFC3, q2, q3)

Get the Fourier transform of the third order matrix element.
Returns flattened, pretransformed matrix element.

This is a faithful port of the Fortran subroutine from thirdorder.f90.
"""
function pretransform_phi3(fct::IFC3, q2::SVector{3, Float64}, q3::SVector{3, Float64})
    nb = fct.na * 3
    ptf = zeros(ComplexF64, nb^3)
    
    for a1 in 1:fct.na
        triplets = get_interactions(fct, a1)
        for trip in triplets
            a2 = trip.idxs[2]
            a3 = trip.idxs[3]
            
            rv2 = trip.rv2
            rv3 = trip.rv3
            
            iqr = dot(q2, rv2) + dot(q3, rv3)
            iqr = -iqr * lo_twopi
            expiqr = cis(iqr)  # exp(i*iqr)
            
            for i in 1:3, j in 1:3, k in 1:3
                ia = (a1 - 1) * 3 + i
                ib = (a2 - 1) * 3 + j
                ic = (a3 - 1) * 3 + k
                # Flattening scheme consistent with zgeru operations
                l = (ia - 1) * nb * nb + (ib - 1) * nb + ic
                ptf[l] += trip.ifcs[i, j, k] * expiqr
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
function pretransform_phi4(fcf::IFC4, q1::SVector{3, Float64}, q2::SVector{3, Float64})
    nb = fcf.na * 3
    ptf = zeros(ComplexF64, nb^4)
    
    for a1 in 1:fcf.na
        quartets = get_interactions(fcf, a1)
        for quart in quartets
            a2 = quart.idxs[2]
            a3 = quart.idxs[3]
            a4 = quart.idxs[4]
            
            rv2 = quart.rv2
            rv3 = quart.rv3
            rv4 = quart.rv4
            
            iqr = -dot(q1, rv2) + dot(q2, rv3) - dot(q2, rv4)
            iqr = iqr * lo_twopi
            expiqr = cis(iqr)
            
            for l in 1:3, k in 1:3, j in 1:3, i in 1:3
                ia = (a1 - 1) * 3 + i
                ib = (a2 - 1) * 3 + j
                ic = (a3 - 1) * 3 + k
                id = (a4 - 1) * 3 + l
                # Flattening scheme consistent with zgeru operations
                m = (ia - 1) * nb^3 + (ib - 1) * nb^2 + (ic - 1) * nb + id
                ptf[m] += quart.ifcs[i, j, k, l] * expiqr
            end
        end
    end
    
    return ptf
end

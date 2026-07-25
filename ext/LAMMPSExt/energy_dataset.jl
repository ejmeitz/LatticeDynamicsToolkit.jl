# make_calc is a function that takes in `sc` and
# returns a LAMMPSCalculator for that structure
function LatticeDynamicsToolkit.make_energy_dataset(
        cc_settings::ConfigSettings,
        uc::CrystalStructure,
        sc::CrystalStructure,
        make_calc::Function;
        ifc2::Union{IFC2, AmorphousIFC2}, # required, but pass as kwarg
        ifc3::Union{Nothing, IFC3} = nothing,
        ifc4::Union{Nothing, IFC4} = nothing,
        n_threads::Integer = Threads.nthreads(),
        verbose::Bool = true
    )

    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))

    if isa(ifc2, AmorphousIFC2) && (ifc3 !== nothing || ifc4 !== nothing)
        error(ArgumentError("Does not make sense to use AmorphousIFC2 with other higher order IFCs to build energy dataset"))
    end
    
    verbose && @info "Remapping IFCs to Supercell"
    valid_ifcs_remapped = remap(sc, uc, valid_ifcs...)
    valid_ifcs_remapped_kwargs = LatticeDynamicsToolkit.build_kwargs(valid_ifcs_remapped...)

    # Build dense long-range polar data once (TDEP fcp equivalent)
    fcp = nothing
    if isa(ifc2, IFC2) && valid_ifcs_remapped_kwargs.ifc2.has_polar_data
        fcp = DensePolarIFCs(valid_ifcs_remapped_kwargs.ifc2, uc, sc)
        verbose && @info "Built dense polar IFCs (LAMMPS dataset)" na=fcp.na lambda=valid_ifcs_remapped_kwargs.ifc2.polar.lambda
    end
    
    return _make_energy_dataset(
            cc_settings,
            sc, make_calc;
            valid_ifcs_remapped_kwargs...,
            fcp = fcp,
            n_threads = n_threads,
            verbose = verbose
        )
end

# Quantum weight: w = sech²(ħω/2kT) = 4n(n+1)/(2n+1)²
function _mode_weight(::QuantumConfigSettings, n_i::Float64)
    denom = 2*n_i + 1
    return 4 * n_i * (n_i + 1) / (denom * denom)
end
_mode_weight(::ClassicalConfigSettings, ::Float64) = 1.0

# Explicit ∂g/∂T factor: (x/T) tanh(x/2), with g = sech²(x/2).
# Classical g ≡ 1 ⇒ derivative vanishes.
function _mode_weight_dT(::ClassicalConfigSettings, ::Float64, ::Float64)
    return 0.0
end
function _mode_weight_dT(::QuantumConfigSettings, T::Float64, ω::Float64)
    if T < LatticeDynamicsToolkit.lo_temperaturetol
        return 0.0
    end
    x = ω / (LatticeDynamicsToolkit.kB_Hartree * T)
    if x > 1e2
        return 0.0
    end
    return (x / T) * tanh(x / 2)
end

# Precompute c_i = w_i * (1/2) * m_i * ω_i² * σ_i² and ∂c_i/∂T (explicit g(T) only).
# Mass factor needed for proper energy units (Hartree)
function _v2_tilde_coefficients(CS::ConfigSettings, freqs_view, sigma_sq, masses)
    T = CS.temperature
    coeffs = zeros(length(freqs_view))
    dcoeffs = zeros(length(freqs_view))
    for (i, (ω, σ², m)) in enumerate(zip(freqs_view, sigma_sq, masses))
        w = _mode_weight(CS, LatticeDynamicsToolkit.planck(T, ω))
        # the `mode_amplitude` function is mass normalized already
        # so we need to multiply by mass here to get correct units
        c = w * 0.5 * m * ω^2 * σ²
        coeffs[i] = c
        dcoeffs[i] = _mode_weight_dT(CS, T, ω) * c
    end
    return coeffs, dcoeffs
end

# Assumes IFCs are supercell already
# Compute true energy given `calc` via AtomsCalculators
function _make_energy_dataset(
    cc_settings::ConfigSettings,
    sc::CrystalStructure,
    make_calc::Function;
    ifc2::Union{IFC2, AmorphousIFC2},
    ifc3::Union{Nothing, IFC3} = nothing,
    ifc4::Union{Nothing, IFC4} = nothing,
    fcp::Union{Nothing,DensePolarIFCs} = nothing,
    n_threads::Integer = Threads.nthreads(),
    D::Int = 3,
    verbose::Bool = true     
)
    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))

    LatticeDynamicsToolkit.remap_checks(sc, valid_ifcs...)

    dynmat = dynmat_gamma(ifc2, sc)
    freqs_sq, phi = get_modes(dynmat, Val{true}())
    freqs = sqrt.(freqs_sq)  # Will error for negative frequencies which I am ok with

    N_atoms = Int(length(freqs) / D)

    # Prepare frequencies and amplitudes
    freqs_view, phi_view, atom_masses = LatticeDynamicsToolkit.prepare(freqs, phi, D, sc.m)
    
    phi_view_T = transpose(phi_view)
    atom_masses_T = transpose(atom_masses)
    mean_amplitude_matrix = LatticeDynamicsToolkit.mean_amplitude.(Ref(cc_settings), freqs_view, atom_masses_T)

    # Pre-scale modes by their amplitudes
    phi_A = phi_view_T .* mean_amplitude_matrix

    # Compute V2_tilde and ∂V2_tilde/∂T coefficients (include mass for proper energy units)
    sigma_sq = mean_amplitude_matrix .^ 2
    coeffs, dcoeffs = _v2_tilde_coefficients(cc_settings, freqs_view, sigma_sq, atom_masses)

    # f(config, z) returns (energies, v2_tilde, dv2_tilde_dT) - closure captures coeffs
    # energies returns (e2, e3, e4, ep); TEP uses (e2+ep, e3, e4)
    verbose && @info "Preparing LAMMPS-backed energy evaluations" has_polar_term=(fcp !== nothing) n_configs=cc_settings.n_configs
    f = (config, z) -> begin
        e2, e3, e4, ep = energies(config, ifc2; fc3=ifc3, fc4=ifc4, fcp=fcp, n_threads=1)
        tep = SVector{4,Float64}(e2, e3, e4, ep)
        v2t = sum(coeffs[i] * z[i]^2 for i in eachindex(z))
        dv2t = sum(dcoeffs[i] * z[i]^2 for i in eachindex(z))
        return (tep, v2t, dv2t)
    end

    # Storage arrays
    tep_energies = zeros(SVector{4, Float64}, cc_settings.n_configs)
    V = zeros(Float64, cc_settings.n_configs)
    V2_tilde = zeros(Float64, cc_settings.n_configs)
    dV2_tilde_dT = zeros(Float64, cc_settings.n_configs)

    # LAMMPSCalculator only supports metal units
    x_cart_eq_ang = copy(sc.x_cart) .* bohr_to_A

    # Make LAMMPSCalculator for each thread
    chnl = Channel{LAMMPSCalculator}(n_threads)
    foreach(1:n_threads) do _
        put!(chnl, make_calc(sc))
    end

    verbose && @info "Building Energy Dataset"
    p = Progress(cc_settings.n_configs; desc="Calculating Energies", dt = 0.25, color = :magenta, enabled = verbose)
    @tasks for n in 1:cc_settings.n_configs
        @set begin
            ntasks = n_threads
            scheduler = :static
        end
        @local begin
            tmp = zeros(size(phi_A))
            coord_storage = zeros(D*N_atoms)
            randn_storage = zeros(D*N_atoms - D)
            f_zero_buf = zeros(D, N_atoms)
        end

        calc = take!(chnl)

        randn!(randn_storage)
        copy!(tmp, phi_A)

        tmp .*= randn_storage

        # Evaluate user function with config AND z values
        coord_storage .= vec(sum(tmp, dims=1))
        cs = reinterpret(SVector{D, Float64}, coord_storage)
        tep_energies[n], V2_tilde[n], dV2_tilde_dT[n] = f(cs, randn_storage)

        # Calculate energy with provided calculator
        # all LAMMPSCalculators use metal units
        cs .*= bohr_to_A
        cs .+= x_cart_eq_ang
        V[n] = single_point_potential_energy(f_zero_buf, cs, calc)

        put!(chnl, calc)
        next!(p)
    end
    verbose && finish!(p)

    return Hartree_to_eV .* tep_energies, V, Hartree_to_eV .* V2_tilde, Hartree_to_eV .* dV2_tilde_dT
end

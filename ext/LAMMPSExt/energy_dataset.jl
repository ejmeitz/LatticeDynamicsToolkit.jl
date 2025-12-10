# make_calc is a function that takes in `sc` and
# returns a LAMMPSCalculator for that structure
function LatticeDynamicsToolkit.make_energy_dataset(
        cc_settings::ConfigSettings,
        uc::CrystalStructure,
        sc::CrystalStructure,
        make_calc::Function;
        ifc2::IFC2, # required, but pass as kwarg
        ifc3::Union{Nothing, IFC3} = nothing,
        ifc4::Union{Nothing, IFC4} = nothing,
        n_threads::Integer = Threads.nthreads()
    )

    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))
    
    @info "Remapping IFCs to Supercell"
    valid_ifcs_remapped = remap(sc, uc, valid_ifcs...)
    valid_ifcs_remapped_kwargs = LatticeDynamicsToolkit.build_kwargs(valid_ifcs_remapped...)
    
    return _make_energy_dataset(cc_settings, sc, make_calc; valid_ifcs_remapped_kwargs..., n_threads = n_threads)
end

#Assumes IFCs are supercell already
# Comptue true energy given `calc` via AtomsCalculators
function _make_energy_dataset(
    cc_settings::ConfigSettings,
    sc::CrystalStructure,
    make_calc::Function;
    ifc2::IFC2,
    ifc3::Union{Nothing, IFC3} = nothing,
    ifc4::Union{Nothing, IFC4} = nothing,
    n_threads::Integer = Threads.nthreads()
)
    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))

    LatticeDynamicsToolkit.remap_checks(sc, valid_ifcs...)

    dynmat = dynmat_gamma(ifc2, sc)
    freqs_sq, phi = get_modes(dynmat, Val{true}())
    freqs = sqrt.(freqs_sq)  # Will error for negative frequencies which I am ok with

    tep_energies = zeros(SVector{3, Float64}, cc_settings.n_configs)

    f = (config) -> energies(config, ifc2; fc3=ifc3, fc4=ifc4, n_threads=1)

    @info "Building Energy Dataset"
    tep_energies, V = _canonical_configs_V!(
        tep_energies,
        f,
        sc,
        make_calc,
        cc_settings,
        freqs,
        phi,
        sc.m;
        n_threads = n_threads
    )

    return Hartree_to_eV .* tep_energies, V

end

# Uses AtomsCalculators to calculate True energies
function _canonical_configs_V!(
        output, 
        f::Function, 
        sc::CrystalStructure, 
        make_calc::Function,
        CM::ConfigSettings, 
        freqs::AbstractVector,
        phi::AbstractMatrix, 
        atom_masses::AbstractVector;
        n_threads::Int = Threads.nthreads(),
        D::Int = 3
    )
    
    N_atoms = Int(length(freqs) / D)

    freqs_view, phi_view, atom_masses = LatticeDynamicsToolkit.prepare(freqs, phi, D, atom_masses)

    phi_view_T = transpose(phi_view)
    atom_masses_T = transpose(atom_masses)
    mean_amplitude_matrix = LatticeDynamicsToolkit.mean_amplitude.(Ref(CM), freqs_view, atom_masses_T) # D*N_atoms - D x D*N_atoms

    # Pre-scale modes by their amplitudes
    phi_A = phi_view_T .* mean_amplitude_matrix # D*N_atoms x D*N_atoms - D

    V = zeros(Float64, CM.n_configs)

    # LAMMPSCalculator only supports metal units
    x_cart_eq_ang = copy(sc.x_cart) .* bohr_to_A

    # Make LAMMPSCalculator for each thread
    chnl = Channel{LAMMPSCalculator}(n_threads)
    foreach(1:n_threads) do _
        put!(chnl, make_calc(sc))
    end

    p = Progress(CM.n_configs; desc="Calculating Energies", dt = 0.25, color = :magenta)
    @tasks for n in 1:CM.n_configs
        @set begin
            ntasks = n_threads
            scheduler = :static
        end
        @local begin
            tmp = zeros(size(phi_A))
            coord_storage = zeros(D*N_atoms)
            randn_storage = zeros(D*N_atoms - D)
        end

        calc = take!(chnl)

        randn!(randn_storage)
        copy!(tmp, phi_A)

        tmp .*= randn_storage

        # Evalulate user function
        coord_storage .= vec(sum(tmp, dims=1))
        cs = reinterpret(SVector{D, Float64}, coord_storage)
        output[n] = f(cs)

        # Calculate energy with provided calculator
        # all LAMMPSCalculators use metal units
        cs .*= bohr_to_A
        cs .+= x_cart_eq_ang
        V[n] = single_point_potential_energy(cs, calc)

        put!(chnl, calc)
        next!(p)
    end
    finish!(p)

    return output, V
end

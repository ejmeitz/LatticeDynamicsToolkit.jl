export 
    canonical_configs, 
    canonical_configs_and_velocities, 
    canonical_velocities
    
function mean_amplitude(qc::QuantumConfigSettings, freq, mass)
    nᵢ = planck(qc.temperature, freq)
    return sqrt((2*nᵢ + 1) / (2 * mass * freq))
end

function mean_amplitude(cc::ClassicalConfigSettings, freq, mass)
    return sqrt((kB_Hartree * cc.temperature) / mass) / freq
end

function extend_masses(atom_masses, D)
    return collect(Iterators.flatten(Iterators.repeated(el, D) for el in atom_masses))
end

function prepare_freqs_phi(freqs, phi, D::Int)
    idx_rt = sortperm(abs.(ustrip.(freqs)))[1:D]
    if idx_rt != 1:D
        throw(ArgumentError("First D modes should be rigid translation modes, got $(idx_rt). There might be imaginary modes present"))
    end

    # View without rigid translation modes
    @views freqs_view = freqs[D+1:end]
    @views phi_view = phi[:, D+1:end]

    return freqs_view, phi_view
end

function prepare(freqs, phi, D, atom_masses)

    freqs_view, phi_view = prepare_freqs_phi(freqs, phi, D)

    # Extend mass vector so there are D copies of each mass per atom
    # This makes broadcasting easier, could also rehspae phi to be (D*N_atoms - D) x N_atoms x 3
    atom_masses = extend_masses(atom_masses, D)

    return freqs_view, phi_view, atom_masses
end

function canonical_configs(CM::ConfigSettings, freqs::AbstractVector,
                         phi::AbstractMatrix, atom_masses::AbstractVector;
                         n_threads::Int = Threads.nthreads(), D::Int = 3)
    
    N_atoms = Int(length(freqs) / D)

    freqs_view, phi_view, atom_masses = prepare(freqs, phi, D, atom_masses)

    phi_view_T = transpose(phi_view)
    atom_masses_T = transpose(atom_masses)
    mean_amplitude_matrix = mean_amplitude.(Ref(CM), freqs_view, atom_masses_T) # D*N_atoms - D x D*N_atoms

    # Create storage
    configs = zeros(D*N_atoms, CM.n_configs)

    # Pre-scale modes by their amplitudes
    phi_A = phi_view_T .* mean_amplitude_matrix # D*N_atoms x D*N_atoms - D


    p = Progress(CM.n_configs; desc="Generating Disps", dt = 0.1, color = :magenta)
    @tasks for n in 1:CM.n_configs
        @set begin
            ntasks = n_threads
            scheduler = :static
        end
        @local begin
            tmp = zeros(size(phi_A))
            randn_storage = zeros(D*N_atoms - D)
        end

        randn!(randn_storage)
        copy!(tmp, phi_A)

        tmp .*= randn_storage
        configs[:, n] .= vec(sum(tmp, dims=1))

        next!(p)
    end
    finish!(p)

    return configs
end

function canonical_configs_and_velocities(CM::ConfigSettings, freqs::AbstractVector,
                                          phi::AbstractMatrix, atom_masses::AbstractVector;
                                          n_threads::Int = Threads.nthreads(), D::Int = 3 )
    
    N_atoms = Int(length(freqs) / D)

    freqs_view, phi_view, atom_masses = prepare(freqs, phi, D, atom_masses)

    phi_view_T = transpose(phi_view)
    atom_masses_T = transpose(atom_masses)
    mean_amplitude_matrix = mean_amplitude.(Ref(CM), freqs_view, atom_masses_T) # D*N_atoms - D x D*N_atoms

    # Create storage
    configs = zeros(D*N_atoms, CM.n_configs)
    velos = zeros(D*N_atoms, CM.n_configs)

    # Pre-scale modes by their amplitudes
    phi_A = phi_view_T .* mean_amplitude_matrix # D*N_atoms x D*N_atoms - D
    freq_unit = unit(first(freqs_view))
    
    p = Progress(CM.n_configs; desc="Generating Disps & Velos", dt = 0.1, color = :yellow)
    # reinterpret_type = typeof(first(phi_A) * freq_unit)
    @tasks for n in 1:CM.n_configs
        @set begin
            ntasks = n_threads
            scheduler = :static
        end
        @local begin
            tmp = zeros(size(phi_A))
            tmp2 = zeros(size(phi_A)) #* this is only needed with unitful
            randn_storage = zeros(DefaultFloat, D*N_atoms - D)
        end

        randn!(randn_storage)
        copy!(tmp, phi_A)
        
        tmp .*= randn_storage
        configs[:, n] .= vec(sum(tmp, dims=1))

        randn!(randn_storage)
        copy!(tmp, phi_A)
        tmp .*= randn_storage
        tmp2 .= tmp .* freqs_view

        velos[:, n] .= vec(sum(tmp2, dims=1))
        next!(p)
    end
    finish!(p)

    return configs, velos
end

function canonical_velocities(CM::ConfigSettings, freqs::AbstractVector,
    phi::AbstractMatrix, atom_masses::AbstractVector;
    n_threads::Int = Threads.nthreads(), D::Int = 3)

    N_atoms = Int(length(freqs) / D)

    freqs_view, phi_view, atom_masses = prepare(freqs, phi, D, atom_masses)

    phi_view_T = transpose(phi_view)
    atom_masses_T = transpose(atom_masses)
    mean_amplitude_matrix = mean_amplitude.(Ref(CM), freqs_view, atom_masses_T) # D*N_atoms - D x D*N_atoms

    # Create storage
    velos = zeros(D*N_atoms, CM.n_configs)

    # Pre-scale modes by their amplitudes
    phi_A = phi_view_T .* mean_amplitude_matrix # D*N_atoms x D*N_atoms - D
    freq_unit = unit(first(freqs_view))
    
    p = Progress(CM.n_configs; desc="Generating Velos", dt = 0.1, color = :green)
    # reinterpret_type = typeof(first(phi_A) * freq_unit)
    @tasks for n in 1:CM.n_configs
        @set begin
            ntasks = n_threads
            scheduler = :static
        end
        @local begin
            # tmp = zeros(eltype(phi_A), size(phi_A))
            tmp2 = zeros(size(phi_A))
            randn_storage = zeros(D*N_atoms - D)
        end

        randn!(randn_storage)
        # copy!(tmp, phi_A)
        # tmp .*= randn_storage
        tmp2 .= (phi_A .* randn_storage) .* freqs_view

        velos[:, n] .= vec(sum(tmp2, dims=1))
        next!(p)
    end
    finish!(p)

    return velos
end

# in place version that evaluates f on each generated config and stores the result in output
# avoids allocating all configurations in RAM. Expects f(::Vector{SVector{3, Float64}})
function canonical_configs!(output, f::Function, CM::ConfigSettings, freqs::AbstractVector,
                         phi::AbstractMatrix, atom_masses::AbstractVector;
                         n_threads::Int = Threads.nthreads(), D::Int = 3)
    
    N_atoms = Int(length(freqs) / D)

    freqs_view, phi_view, atom_masses = prepare(freqs, phi, D, atom_masses)

    phi_view_T = transpose(phi_view)
    atom_masses_T = transpose(atom_masses)
    mean_amplitude_matrix = mean_amplitude.(Ref(CM), freqs_view, atom_masses_T) # D*N_atoms - D x D*N_atoms

    # Pre-scale modes by their amplitudes
    phi_A = phi_view_T .* mean_amplitude_matrix # D*N_atoms x D*N_atoms - D

    p = Progress(CM.n_configs; desc="Generating Disps", dt = 0.25, color = :magenta)
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

        randn!(randn_storage)
        copy!(tmp, phi_A)

        tmp .*= randn_storage

        coord_storage .= vec(sum(tmp, dims=1))
        output[n] = f(reinterpret(SVector{D, Float64}, coord_storage))

        next!(p)
    end
    finish!(p)

    return output
end


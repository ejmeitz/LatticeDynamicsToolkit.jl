
"""
    sTDEP(
        sys::CrystalStructure,
        calc::LAMMPSCalculator,
        basedir::String,
        niter::Integer,
        nconf::Union{Integer, AbstractVector{Integer}},
        rc2,
        temperature,
        maximum_frequency;
        quantum::Bool = false,
        ncores::Integer = Threads.nthreads(),
        verbose::Bool = false
    )

Calculates force constants self-consistently by sampling configurations from the canonical ensemble,
fitting force constants and iterating until convergence. The harmonic free energy, dispersion and DOS
are computed for each iteration and can be used to asses convergence. If `rc3` and `rc4` are provided,
third and fourth order force constants are fitted and used to compute the harmonic free energy. `rc3` 
can be provided without `rc4`. Results are moved to the `basedir` directory.

This function is written to take in any `CrystalStructure` from the LatticeDynamicsToolkit interface and any force calculator
from the LAMMPSCalculator interface.

# Arguments:
- `sys` : CrystalStructure containing at a minimum the positions and cell definition.
- `calc` : LAMMPSCalculator
- `basedir::String` : Directory where all output files are written
- `niter::Integer` : Number of self-consistent iterations to perform. Includes iter to initialize from freuency.
- `nconf` : If Integer, the number of configurations to generate. If list of integers, the number
     of configuritons to generate at each iteration.
- `rc2` : Force constant cutoff in Angstroms when calculating second-order IFCs from configurations
- `temperature` : Temperature to generate configurations at
- `maximum_frequency` : Estimate of maximum freuqency in THz to generate initial configurations from
- `quantum::Bool = false` : Generate configurations from the quantum canonical ensemble
- `F_H_mesh::Vector{Int} = [30, 30, 30]` : Mesh used to compute the harmonic free energy
- `rc3::Union{Float64, Nothing} = nothing` : Cutoff used to fit third-order IFCs on the last iteration (Angstrom)
- `rc4::Union{Float64, Nothing} = nothing` : Cutoff used to fit fourth-order IFCs on the last iteration (Angstrom)
- `ncores::Integer = Threads.nthreads()` : Number of cores used by IFC calculation
- `verbose::Bool = false` : Enable extra printing
"""
function LatticeDynamicsToolkit.sTDEP(
    sys::CrystalStructure,
    calc::LAMMPSCalculator,
    basedir::String,
    niter::Integer,
    rc2,
    temperature,
    maximum_frequency;
    mix::Bool = true,
    nconf_init::Int = 8,
    max_configs::Int = 512,
    quantum::Bool = false,
    F_H_mesh = [30, 30, 30],
    rc3::Union{Float64, Nothing} = nothing,
    rc4::Union{Float64, Nothing} = nothing,
    ncores::Integer = Threads.nthreads(),
    verbose::Bool = false,
)

    @assert isfile(joinpath(basedir, "infile.ucposcar")) "sTDEP requires an infile.ucposcar in $(basedir)"

    nconf = []
    for i in 1:niter
        nconf_i = nconf_init * 2^(i-1)
        if nconf_i > max_configs
            nconf_i = max_configs
        end
        push!(nconf, nconf_i)
    end

    # Make ssposcar
    TDEPWrapper.write_ssposcar(basedir, sys.L .* bohr_to_A, sys.x_cart .* bohr_to_A, sys.species)

    get_path = (i) -> joinpath(basedir, "iter$(lpad(i,3,'0'))")

    efc = ExtractForceConstants(secondorder_cutoff = rc2)
    if quantum
        pd = PhononDispersionRelations(dos = true, temperature = Float64(temperature))
    else
        pd = PhononDispersionRelations(dos = true)
    end

    init_dir = get_path(0)
    mkdir(init_dir)

    # Generate initial configurations with maximum freuqency
    cp(joinpath(basedir, "infile.ssposcar"), joinpath(init_dir, "infile.ssposcar"))
    cp(joinpath(basedir, "infile.ucposcar"), joinpath(init_dir, "infile.ucposcar"))

    # My config code does not support the max-freuqency option
    # So we will use TDEP to generate the configurations
    cc_init = CanonicalConfiguration(
        temperature = temperature,
        nconf = nconf[1],
        quantum = quantum, 
        maximum_frequency = maximum_frequency
    )

    generate_configs(sys, cc_init, calc, init_dir, verbose)

    F_Hs = zeros(Float64, niter)
    F_Hs[1] = compute_harmonic_free_energy(
        F_H_mesh,
        init_dir,
        temperature,
        ncores,
        quantum;
        ifc_name = "outfile.fakeforceconstant"
    )

    # Generate remaining configurations with IFCs from prior iteration
    p = Progress(niter - 1, desc = "sTDEP IFCs")
    for i in 1:(niter - 1)

        cc = CanonicalConfiguration(
            temperature = temperature,
            nconf = nconf[i+1],
            quantum = quantum, 
        )

        outdir = get_path(i)
        # Move IFCs from last iter to current dir
        prepare_next_dir(get_path(i-1), outdir, mix, i == 1)
        # Generate Configs Given Current IFCs
        nconf_extra = mix ? nconf[i] : 0 # if mixing account for extras
        generate_configs(sys, cc, calc, outdir, verbose; nconf_extra = nconf_extra)
        # Calculate IFCs to Generate Next Set of Configs
        execute(efc, outdir, ncores, verbose)
        # Generate DOS and Dispersion Data
        execute(pd, outdir, ncores, verbose)
        # Compute Harmonic Free Energy
        F_Hs[i+1] = compute_harmonic_free_energy(F_H_mesh, outdir, temperature, ncores, quantum)
        next!(p)
    end
    finish!(p)

    final_iter_dir = get_path(niter - 1)

    if rc4 !== nothing && rc3 === nothing
        rc3 = rc4
        @warn "rc4 provided but not rc3, using rc3 = rc4"
    end

    # Fit Higher-Order IFCs
    efc =
    if rc3 !== nothing && rc4 === nothing
        ExtractForceConstants(
            secondorder_cutoff = rc2,
            thirdorder_cutoff = rc3,
        )
    elseif rc4 !== nothing
        ExtractForceConstants(
            secondorder_cutoff = rc2,
            thirdorder_cutoff = rc3,
            fourthorder_cutoff = rc4,
        )
    else
        nothing
    end

    if efc !== nothing
        t0 = time_ns()
        execute(efc, final_iter_dir, ncores, verbose)
        dt = (time_ns() - t0) / 1e9
        verbose && @info "3rd/4th order IFCs took $(dt) seconds"
    end

    open(joinpath(basedir, "harmonic_free_energy.txt"), "w") do f
        println(f, "# Iteration - F_Harmonic [eV / atom] at T = $(temperature) K")
        for i in 1:niter
            @printf f "%d %.15f\n" i F_Hs[i]*Hartree_to_eV
        end
    end

    # Move results to basedir
    cp(joinpath(final_iter_dir, "outfile.forceconstant"), joinpath(basedir, "infile.forceconstant"))
    rc3 !== nothing && cp(joinpath(final_iter_dir, "outfile.forceconstant_thirdorder"), joinpath(basedir, "infile.forceconstant_thirdorder"))
    rc4 !== nothing && cp(joinpath(final_iter_dir, "outfile.forceconstant_fourthorder"), joinpath(basedir, "infile.forceconstant_fourthorder"))

    return nothing
end

function compute_harmonic_free_energy(
    F_H_mesh::Vector{Int},
    dir::String,
    T::Float64,
    ncores::Integer,
    quantum::Bool;
    ifc_name::String = "outfile.forceconstant"
)

    L = quantum ? Quantum : Classical
    uc = CrystalStructure(joinpath(dir, "infile.ucposcar"))
    ibz = SimpleMesh(uc, F_H_mesh)
    ifc2 = read_ifc2(
        joinpath(dir, ifc_name),
        joinpath(dir, "infile.ucposcar")
    )
    
    F0 = LatticeDynamicsToolkit.sum_over_freqs(
        (ω) -> F_harmonic_single(ω, kB_Hartree*T, L), 
        ibz,
        uc,
        ifc2;
        n_threads = ncores
    )

    return F0 / ibz.n_atoms_prim

end

function prepare_next_dir(current_dir, dest_dir, mix::Bool, init_pass::Bool = false)
    mkdir(dest_dir)
    cp(joinpath(current_dir, "infile.ssposcar"), joinpath(dest_dir, "infile.ssposcar"))
    cp(joinpath(current_dir, "infile.ucposcar"), joinpath(dest_dir, "infile.ucposcar"))

    if init_pass
        cp(joinpath(current_dir, "outfile.fakeforceconstant"), joinpath(dest_dir, "infile.forceconstant"))
    else
        cp(joinpath(current_dir, "outfile.forceconstant"), joinpath(dest_dir, "infile.forceconstant"))
    end

    if mix
        cp(joinpath(current_dir, "infile.positions"), joinpath(dest_dir, "infile.positions"))
        cp(joinpath(current_dir, "infile.forces"), joinpath(dest_dir, "infile.forces"))
        cp(joinpath(current_dir, "infile.stat"), joinpath(dest_dir, "infile.stat"))
    end
end

function generate_configs(
        sys::CrystalStructure,
        cc::CanonicalConfiguration, 
        calc::LAMMPSCalculator,
        outdir::String,
        verbose::Bool;
        nconf_extra::Int = 0
    )

    n_atoms = length(sys)
    cell_ang = sys.L .* bohr_to_A

    execute(cc, outdir, 1, verbose)

    get_filepath = (i) -> joinpath(outdir, "contcar_conf$(lpad(i, 4, '0'))")

    f_buf = zeros(Float64, 3, n_atoms)
    f_zeros_buf = zeros(Float64, 3, n_atoms)

    # Parse coordinates into sys object and calculate forces
    for i in 1:cc.nconf
        filepath = get_filepath(i)
        x_cart, _ = TDEPWrapper.read_poscar_positions(filepath, n_atoms = n_atoms)
        x_frac = LatticeDynamicsToolkit.to_frac_coords.(Ref(cell_ang), x_cart) #! Allocation

        PE, f_buf = single_point_forces_and_energy!(f_buf, f_zeros_buf, x_cart, calc)

        # Add data to infile.forces
        open(joinpath(outdir, "infile.forces"), "a") do ff
            for j in 1:n_atoms
                @printf ff "%.15f %.15f %.15f\n" view(f_buf, :, j)...
            end
        end

        # Add data to infile.positions
        open(joinpath(outdir, "infile.positions"), "a") do pf
            for j in eachindex(x_frac)
                @printf pf "%.15f %.15f %.15f\n" x_frac[j]...
            end
        end

        # Add PE to infile.stat
        write_partial_stat(
            outdir, 
            [PE],
            [0.0],
            [cc.temperature];
            file_mode = "a"
        )
    end

    # Write infile.meta
    write_meta(outdir, ustrip.(cc.temperature), cc.nconf + nconf_extra, 1.0, n_atoms)

end

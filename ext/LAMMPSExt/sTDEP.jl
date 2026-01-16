
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
fitting force constants and iterating until convergence. 

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
    ncores::Integer = Threads.nthreads(),
    verbose::Bool = false
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
    write_ssposcar(basedir, sys.L, sys.x_frac, sys.species)

    get_path = (i) -> joinpath(basedir, "iter$(lpad(i,3,'0'))")

    efc = ExtractForceConstants(secondorder_cutoff = rc2)
    pd = PhononDispersionRelations(dos = true)

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

    # Generate remaining configurations with IFCs from prior iteration
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
        nconf_extra = mix ? nconf[i-1] : 0
        generate_configs(sys, cc, calc, outdir, verbose; nconf_extra = nconf_extra)
        # Calculate IFCs to Generate Next Set of Configs
        execute(efc, outdir, ncores, verbose)
        # Generate DOS and Dispersion Data
        execute(pd, outdir, ncores, verbose)
    end

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
        cp(joinpath(currentdir, "infile.positions"), joinpath(dest_dir, "infile.positions"))
        cp(joinpath(currentdir, "infile.forces"), joinpath(dest_dir, "infile.forces"))
        cp(joinpath(currentdir, "infile.stat"), joinpath(dest_dir, "infile.stat"))
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

    @info "Generating Configurations"
    execute(cc, outdir, 1, verbose)

    get_filepath = (i) -> joinpath(outdir, "contcar_conf$(lpad(i, 4, '0'))")

    # Parse coordinates into sys object and calculate forces
    p = Progress(cc.nconf, desc = "Calculating Forces")
    for i in 1:cc.nconf
        filepath = get_filepath(i)
        x_cart, _ = read_poscar_positions(filepath, n_atoms = length(sys))

        PE, F = single_point_forces_and_energy(x_cart, calc)

        # Add data to infile.forces
        open(joinpath(outdir, "infile.forces"), "a") do ff
            for j in eachindex(F)
                @printf ff "%.15f %.15f %.15f\n" F[j]...
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

        next!(p)
    end
    finish!(p)

    # Write infile.meta
    write_meta(outdir, ustrip.(cc.temperature), cc.nconf + nconf_extra, 1.0, length(sys))

end

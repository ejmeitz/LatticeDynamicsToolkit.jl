
# Will not initialize velocities
function LatticeDynamicsToolkit.LAMMPSCalculator(
        sys::CrystalStructure,
        potential_definition::Union{String, Vector{String}};
        label_type_map::Dict{Symbol, Int} = Dict{Symbol, Int}(),
        logfile_path::String = "none",
        unit_system::String = "metal",
        single_point::Bool = true
    )

    if !LinearAlgebra.isdiag(sys.L)
        error(ArgumentError("LAMMPSCalculator does not support triclinic systems yet. Got non-diagonal cell."))
    end

    if length(potential_definition) == 0
        error(ArgumentError("Cannot pass emptry string as potential definition for LAMMPSCalculator"))
    end

    if unit_system != "metal" 
        error(ArgumentError("Currently only 'metal' unit system is supported by LAMMPSCalculator"))
    end

    LAMMPS.MPI.Init()
    lmp = LMP(["-screen","none"], LAMMPS.MPI.COMM_WORLD)

    all_syms = AtomsBase.atomic_symbol(sys, :)
    unique_syms = unique(all_syms)
    unique_sym_idxs = Dict(sym => findfirst(x -> x == sym, all_syms) for sym in unique_syms)

    if any(unique_syms .== :unknown)
        error(ArgumentError("All atoms must have atomic symbols to use LAMMPSCalculator"))
    end
    
    ids = collect(Int32, 1:length(sys))
    # convert to Angstroms to be comensurate with metal units
    xhi, yhi, zhi = bohr_to_A .* LinearAlgebra.diag(sys.L)

    if length(label_type_map) == 0
        label_type_map = Dict(sym => Int32(i) for (i, sym) in enumerate(unique_syms))
        types = [label_type_map[sym] for sym in all_syms]
    else 
        unique_sym_user = keys(label_type_map)
        if Set(unique_sym_user) != Set(unique_syms)
            error(ArgumentError("You provided a label_type_map with $(unique_sym_user) symbols, but" *
                " the system has $(unique_syms). They must match exactly if you pass label_type_map."))
        end
        types = [Int32(label_type_map[sym]) for sym in all_syms]
    end

    m_lmp = Dict(label_type_map[sym] => emu_to_amu * sys.m[i] for (sym, i) in unique_sym_idxs)

    label_map_cmd = "labelmap atom " * join(["$(i) $(sym)" for (sym,i) in label_type_map], " ") 

    # Figure out a safe skin distance: nearest neighbor distance
    # in crystals atoms should never move more than this distance
    # so we the initial neighbor list will always be valid.
    r_cut = (min(xhi, yhi, zhi) / 2.0) - lo_tol # biggest possible cutoff given cell
    nl = CellListMap.neighborlist(sys.x_cart, r_cut * A_to_bohr; unitcell = sys.L)
    skin_distance = bohr_to_A * minimum(x -> last(x), nl)

    setup_cmd = """
            log $(logfile_path)
            units $(unit_system)
            atom_style atomic
            atom_modify map array sort 0 0
            neighbor $(skin_distance) bin
        """
    
    cell_cmd = """
            boundary p p p
            region cell block 0 $(xhi) 0 $(yhi) 0 $(zhi) units box
            create_box $(length(unique_syms)) cell
            $(label_map_cmd)
        """

    mass_cmd = join(["mass $(type) $(m)" for (type,m) in  m_lmp], "\n")

    command(lmp, setup_cmd)
    command(lmp, cell_cmd)
    command(lmp, mass_cmd)


    LAMMPS.create_atoms(
        lmp,
        reinterpret(reshape, Float64, sys.x_cart .* bohr_to_A),
        ids,
        types
    )   

    try
        command(lmp, potential_definition)
    catch e
        if startswith(e.msg, "Number of element to type mappings does")
            @info "Ensure path to potential definition is wrapped in quotes if there are spaces in path."
        end
        rethrow(e)
    end

    command(lmp, "compute pot_e all pe")

    # Little hack to make single point calculations
    # fast and avoid re-building neighborlist a lot
    if single_point
        command(lmp, "fix hold all nve")
        command(lmp, "fix freeze all setforce 0.0 0.0 0.0")
        command(lmp, "velocity all set 0 0 0")
        command(lmp, "timestep 0.01") # timestep won't matter since we freeze atoms
    end # if not single_point user is responsible for doing things like initializing velos, see TI.jl

    # This allows LAMMPS to register the computes/fixes
    # and build the neighbor list. 
    command(lmp, "run 0 post no")

    return LatticeDynamicsToolkit.LAMMPSCalculator{typeof(lmp), typeof(potential_definition)}(lmp, potential_definition)
end

AtomsCalculators.energy_unit(inter::LAMMPSCalculator) = NoUnits

# Assumes sys.r_cart already in right units 
function AtomsCalculators.potential_energy(sys::CrystalStructure, inter::LAMMPSCalculator; kwargs...)
    scatter!(inter.lmp, "x", reinterpret(reshape, Float64, sys.r_cart))
    command(inter.lmp, "run 1 pre no post yes")
    return extract_compute(inter.lmp, "pot_e", STYLE_GLOBAL, TYPE_SCALAR)[1]
end

# Expect Vector of Vectors or 3 x N Matrix
function single_point_potential_energy(r::AbstractVecOrMat, inter::LAMMPSCalculator)
    scatter!(inter.lmp, "x", reinterpret(reshape, Float64, r))
    command(inter.lmp, "run 1 pre no post yes")
    return extract_compute(inter.lmp, "pot_e", STYLE_GLOBAL, TYPE_SCALAR)[1]
end

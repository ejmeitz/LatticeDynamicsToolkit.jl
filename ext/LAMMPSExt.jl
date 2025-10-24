module LAMMPSExt

using LAMMPS
using TDEPToolkit
using AtomsBase
using AtomsCalculators
using LinearAlgebra
using Unitful

function TDEPToolkit.LAMMPSCalculator(
        sys::CrystalStructure,
        potential_definition::Union{String, Array{String}};
        label_type_map::Dict{Symbol, Int} = Dict{Symbol, Int}(),
        logfile_path::String = "none",
        unit_system::String = "metal"
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

    m_lmp = Dict(label_type_map[sym] => ustrip(sys.mass[i]) for (sym, i) in unique_sym_idxs)

    label_map_cmd = "labelmap atom " * join(["$(i) $(sym)" for (sym,i) in label_type_map], " ") 

    # Figure out a safe skin distance: nearest neighbor distance
    # in crystals atoms should never move more than this distance
    # so we the initial neighbor list will always be valid.
    r_cut = (min(xhi, yhi, zhi) / 2.0) - lo_tol # biggest possible cutoff given cell
    nl = CellListMap.neighborlist(sc.x_cart, r_cut * A_to_bohr; unitcell = sc.L)
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

    m_lmp .*= TDEPToolkit.emu_to_amu
    mass_cmd = join(["mass $(type) $(m)" for (type,m) in  m_lmp], "\n")

    command(lmp, setup_cmd)
    command(lmp, cell_cmd)
    command(lmp, mass_cmd)


    LAMMPS.create_atoms(
        lmp,
        reinterpret(reshape, Float64, sys.position),
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

    # This allows LAMMPS to register the computes/fixes
    # and build the neighbor list. 
    command(lmp, "run 0 post no")

    return TDEP.LAMMPSCalculator{typeof(lmp)}(lmp, -1)
end

AtomsCalculators.energy_unit(inter::LAMMPSCalculator) = NoUnits

# Assumes sys.r_cart already in right units 
function AtomsCalculators.potential_energy(sys::CrystalStructure, inter::LAMMPSCalculator; kwargs...)
    scatter!(inter.lmp, "x", reinterpret(reshape, Float64, sys.r_cart))
    command(inter.lmp, "run 0 pre no post no")
    return extract_compute(inter.lmp, "pot_e", STYLE_GLOBAL, TYPE_SCALAR)[1]
end

# Expect Vector of Vectors or 3 x N Matrix
function single_point_potential_energy(r::AbstractVecOrMat, inter::LAMMPSCalculator)
    scatter!(inter.lmp, "x", reinterpret(reshape, Float64, r))
    command(inter.lmp, "run 0 pre no post no")
    return extract_compute(inter.lmp, "pot_e", STYLE_GLOBAL, TYPE_SCALAR)[1]
end


function TDEPToolkit.make_energy_dataset(
        cc_settings::ConfigSettings,
        uc::CrystalStructure,
        sc::CrystalStructure,
        calc::LAMMPSCalculator;
        ifc2::IFC2, # required, but pass as kwarg
        ifc3::Union{Nothing, IFC3} = nothing,
        ifc4::Union{Nothing, IFC4} = nothing,
        n_threads::Integer = Threads.nthreads()
    )

    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))
    
    @info "Remapping IFCs to Supercell"
    valid_ifcs_remapped = remap(sc, uc, valid_ifcs...)
    valid_ifcs_remapped_kwargs = build_kwargs(valid_ifcs_remapped...)
    
    return _make_energy_dataset(cc_settings, sc, calc; valid_ifcs_remapped_kwargs..., n_threads = n_threads)
end

#Assumes IFCs are supercell already
# Comptue true energy given `calc` via AtomsCalculators
function _make_energy_dataset(
    cc_settings::ConfigSettings,
    sc::CrystalStructure,
    calc::LAMMPSCalculator;
    ifc2::IFC2,
    ifc3::Union{Nothing, IFC3} = nothing,
    ifc4::Union{Nothing, IFC4} = nothing,
    n_threads::Integer = Threads.nthreads()
)
    valid_ifcs = Iterators.filter(!isnothing, (ifc2, ifc3, ifc4))

    remap_checks(sc, valid_ifcs...)

    dynmat = dynmat_gamma(ifc2, sc)
    freqs_sq, phi = get_modes(dynmat)
    freqs = sqrt.(freqs_sq)  # Will error for negative frequencies which I am ok with

    tep_energies = zeros(SVector{3, Float64}, cc_settings.n_configs)

    f = (config) -> energies(config, ifc2; fc3=ifc3, fc4=ifc4, n_threads=1)

    @info "Building Energy Dataset"
    tep_energies, V = canonical_configs_V!(
        tep_energies,
        f,
        sc,
        calc,
        cc_settings,
        freqs,
        phi,
        sc.m;
        n_threads = n_threads
    )

    return Hartree_to_eV .* tep_energies, V .* Hartree_to_eV

end

# Uses AtomsCalculators to calculate True energies
function canonical_configs_V!(
        output, 
        f::Function, 
        sc::CrystalStructure, 
        calc::LAMMPSCalculator,
        CM::ConfigSettings, 
        freqs::AbstractVector,
        phi::AbstractMatrix, 
        atom_masses::AbstractVector;
        n_threads::Int = Threads.nthreads(),
        D::Int = 3
    )
    
    N_atoms = Int(length(freqs) / D)

    freqs_view, phi_view, atom_masses = prepare(freqs, phi, D, atom_masses)

    phi_view_T = transpose(phi_view)
    atom_masses_T = transpose(atom_masses)
    mean_amplitude_matrix = mean_amplitude.(Ref(CM), freqs_view, atom_masses_T) # D*N_atoms - D x D*N_atoms

    # Pre-scale modes by their amplitudes
    phi_A = phi_view_T .* mean_amplitude_matrix # D*N_atoms x D*N_atoms - D

    V = zeros(Float64, CM.n_configs)

    # LAMMPSCalculator only supports metal units
    x_cart_eq_ang = copy(sc.x_cart) .* bohr_to_A

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
            sys = deepcopy(sc)
        end

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

        next!(p)
    end
    finish!(p)

    return output, V
end


end
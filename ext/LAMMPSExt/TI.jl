# just trusts that this is supercell representation of IFCs
function (::Type{<:Matrix})(ifc2::IFC2, sc::CrystalStructure)

    ifc2.na != length(sc) && error(ArgumentError("Cannot convert IFC to Matrix. IFC are built on $(ifc2.na) cell, but supercell has $(length(sc)). Try remapping IFCs to the supercell."))


    n_dof = 3*ifc2.na
    out = zeros(Float64, n_dof, n_dof)
    
    @inbounds for a1 in 1:ifc2.na
        r1 = 3*(a1-1)
        for pair in LatticeDynamicsToolkit.get_interactions(ifc2, a1)
            a2 = pair.idxs[2]
            r2 = 3*(a2-1)
            out[r1+1:r1+3, r2+1:r2+3] .= pair.ifcs 
        end
    end

    # enforce exact symmetry
    @inbounds for j in 1:n_dof, i in j+1:n_dof
        s = 0.5 * (out[i,j] + out[j,i])
        out[i,j] = s; out[j,i] = s
    end

    # enforce acoustic sum rule
    for i in 1:ifc2.na # index of block matrix
        for α in 1:3
            for β in 1:3
                ii = 3*(i-1) + α
                jj = 3*(i-1) + β # i == j because we're on diagonal
                out[ii,jj] = 0.0 # remove existing value
                out[ii,jj] = @views -1*sum(out[ii, β:3:end])
            end
        end
    end

    return out
end

_idx(i) = (3*(i-1) + 1) : (3*i)

function _pack_displacements!(u, ux, uy, uz)
    @inbounds for i in eachindex(ux)
        j = 3*i - 2
        u[j]   = ux[i]
        u[j+1] = uy[i]
        u[j+2] = uz[i]
    end
    return u
end

# Assumes only 1 MPI rank since we have global phi
# could communicate rows of phi but why...
function FixTIExternal(
        ifc2::Matrix{Float64},
        ΔU,
        lmp::LMP,
        n_equil::Integer,
        λ::Float64,
        V₀,
        p::Progress
    )

    # command(lmp, "compute disp all displace/atom")

    function compute_harmonic_ef_ti!(f, u)
        mul!(f, ifc2, u)
        e = 0.5 * dot(u, f)
        @. f = -f
        return (e + V₀), f
    end

    N_atoms = Int(size(ifc2, 1) / 3)

    # Buffers
    u = zeros(Float64, size(ifc2, 1))
    u_lmp = zeros(Float64, 4, N_atoms) # x, y, z, mag
    f_lmp = zeros(Float64, 3, N_atoms) # x, y, z
    f_harm = zeros(Float64, 3*N_atoms)
    count = 1

    # This runs in the "post_force" part of the timestep
    # All other forces are known, final integration 
    # update is not complete yet
    FixExternal(lmp, "ti_julia", "all", 1, 1) do fix::LAMMPS.FixExternal

        atom = LAMMPS.ExtractAtomMultiple{@NamedTuple{id::Int32}}(fix.lmp)

        # These forces are added to existing lammps force buffer
        # automatically (see fix external docs)
        f = vec(fix.f)
        fill!(f, zero(eltype(f)))

        # Get lammps internal force/energy data from other potentials defined before this fix
        gather!(fix.lmp, "f", f_lmp) # ordered by id 
        f_lmp_vec = vec(f_lmp) # re-arrange to 3N x 1

        gather!(fix.lmp, "c_disp", u_lmp) # ordered by id
        # u_lmp = extract_compute(fix.lmp, "disp", STYLE_ATOM, TYPE_ARRAY) #! is this faster?
        @views _pack_displacements!(u, u_lmp[1,:], u_lmp[2,:], u_lmp[3,:])

        e_harm, f_harm = compute_harmonic_ef_ti!(f_harm, u)

        # fix external adds `f` to the existing lammps force buffer.
        # We want the forces to be λ*f_lmp + (1-λ)*f_harm
        # so we add (1-λ)*(f_harm[i] - f_lmp[i])
        # f_lmp + (1-λ)*(f_harm[i] - f_lmp[i]) = λ*f_lmp + (1-λ)*f_harm
        @inbounds for i in 1:N_atoms
            id = atom[i].id
            for j in 1:3
                idx1 = 3*(i-1) + j
                idx2 = 3*(id-1) + j
                f[idx1] = (1.0 - λ) * (f_harm[idx2] - f_lmp_vec[idx2])
            end
        end

        # compute pot_e is defined in LAMMPSCalculator
        e_lmp = extract_compute(fix.lmp, "pot_e", STYLE_GLOBAL, TYPE_SCALAR)[1]

        if fix.timestep > n_equil
            ΔU[count] = e_lmp - e_harm
            count += 1
        end

        # DO NOT UPDATE THE POTENTIAL ENERGY IN LAMMPS
        # extract_compute should only read the energy of the user defined potential

        LAMMPS.check(fix.lmp)
        next!(p)

    end

end

function _TI_step(
        lc::LAMMPSCalculator, 
        s::TISettings,
        ifcs::Matrix{Float64},
        x0,
        λ,
        V₀,
        p::Progress
    )

    command(lc.lmp, "reset_timestep 0")
    command(lc.lmp, "neigh_modify every 1 delay 0 check yes")
    command(lc.lmp, "compute disp all displace/atom")

    scatter!(lc.lmp, "x", reinterpret(reshape, Float64, x0))
    seed = rand(1000:1_000_000)
    command(lc.lmp, "velocity all create $(s.T) $(seed) dist gaussian mom yes")
    
    # This makes energy available inside FixTIExternal
    command(lc.lmp, "thermo 1")
    command(lc.lmp, "thermo_style custom step pe")

    command(lc.lmp, "fix 1 all nvt temp $(s.T) $(s.T) $(s.T_damp)")

    ΔU = zeros(Float64, s.nsteps)
    FixTIExternal(ifcs, ΔU, lc.lmp, s.nsteps_equil, λ, V₀, p)

    command(lc.lmp, "timestep $(s.dt_ps)")
    command(lc.lmp, "run $(s.nsteps + s.nsteps_equil)")

    command(lc.lmp, "unfix 1")
    command(lc.lmp, "uncompute disp") # Defined in FixTIExternal
    command(lc.lmp, "unfix ti_julia") # Defined in FixTIExternal

    return sum(ΔU)/length(ΔU)

end


# IGNORES POLAR INTERACTIONS
# returns free energy in eV
function LatticeDynamicsToolkit.ThermodynmicIntegration(
        ifc2_uc::IFC2,
        sc::CrystalStructure,
        uc::CrystalStructure,
        s::TISettings;
        n_threads::Integer = Threads.nthreads(),
        V₀::Float64 = 0.0,
        kwargs...
    )

    if ifc2_uc.na == length(sc)
        error(ArgumentError("IFCs passed to TI are already remapped to the supercell. Please pass the unitcell IFCs."))
    end

    ifc2_sc = remap(sc, uc, ifc2_uc)[1]

    @assert get(kwargs, :single_point, false) == false "Do not pass 'single_point' kwarg to TI, it is manually set to false."

    ifc2_lmp = forceconstant_2nd_HartreeBohr_to_eVA .* Matrix(ifc2_sc, sc)
    
    x0_ang = bohr_to_A .* sc.x_cart

    # Figure out integration weights
    λs, weights = gausslegendre(s.n_lambda)
    integrands = zeros(Float64, length(λs))

    if !all(-1.0 .<= λs .<= 1.0)
        throw(ArgumentError("Quadrature rule must emit points in [-1, 1]."))
    end

    # Map to [0,1]
    λs = (λs .+ 1) ./ 2
    weights ./= 2

    # Global progress bar updated by each thread
    p = Progress(
        s.n_lambda * (s.nsteps + s.nsteps_equil);
        desc = "Thermodynamic Integration",
        dt = 0.5,
    )

    #! is there a way to ensure each lammps calculator only has one core?
    @tasks for (i,λ) in collect(enumerate(λs))
        @set ntasks = n_threads
        @local lc = LAMMPSCalculator(sc, s.pot_cmds; single_point = false, kwargs...)
        integrands[i] = _TI_step(lc, s, ifc2_lmp, x0_ang, λ, V₀, p)
    end

    finish!(p)

    @info "Calculating harmonic free energy on 30x30x30 mesh"

    # This F0 is per atom
    F₀, _, _, _ = harmonic_properties(
        s.T, 
        uc,
        ifc2_uc,
        [30,30,30],
        Classical;
        n_threads = n_threads
    )

    F₀ *= Hartree_to_eV

    ΔF = dot(weights, integrands) / length(sc)

    return F₀, ΔF # eV / atom
end

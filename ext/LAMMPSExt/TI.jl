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
        dU_dλ,
        lmp::LMP,
        n_equil::Integer,
        λ::Float64 = 0.0
    )

    command(lmp, "compute disp all displace/atom")

    function compute_harmonic_ef_ti!(f, u)
        mul!(f, ifc2, u)
        e = 0.5 * dot(u, f)
        @. f = -f
        return e, f
    end

    N_atoms = Int(size(ifc2, 1) / 3)

    # Buffers
    u = zeros(Float64, size(ifc2, 1))
    u_lmp = zeros(Float64, 4, N_atoms) # x, y, z, mag
    f_lmp = zeros(Float64, 3, N_atoms) # x, y, z
    f_harm_λ = zeros(Float64, 3*N_atoms)
    count = 1

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
        @views _pack_displacements!(u, u_lmp[1,:], u_lmp[2,:], u_lmp[3,:])

        e_harm, f_harm_λ = compute_harmonic_ef_ti!(f_harm_λ, u)

        # fix external adds `f` to the existing lammps force buffer.
        # We want the forces to be λ*f_lmp + (1-λ)*f_harm
        # so we add (1-λ)*(f_harm[i] - f_lmp[i])
        # f_lmp + (1-λ)*(f_harm[i] - f_lmp[i]) = λ*f_lmp + (1-λ)*f_harm
        @inbounds for i in 1:N_atoms
            id = atom[i].id
            for j in 1:3
                idx1 = 3*(i-1) + j
                idx2 = 3*(id-1) + j
                f[idx1] = (1.0 - λ) * (f_harm_λ[idx2] - f_lmp_vec[idx2])
            end
        end

        # compute pot_e is defined in LAMMPSCalculator
        e_lmp = extract_compute(fix.lmp, "pot_e", STYLE_GLOBAL, TYPE_SCALAR)[1]
        if fix.timestep > n_equil
            dU_dλ[count] = e_lmp - e_harm
            count += 1
        end

        # # energy as a result of this fix has similar argument to force
        U_fix = (1.0 - λ) * (e_harm - e_lmp)
        LAMMPS.API.lammps_fix_external_set_energy_global(fix.lmp, fix.name, U_fix)
        LAMMPS.check(fix.lmp)
    end
end

function _TI_step(
        lc::LAMMPSCalculator, 
        λ,
        T,
        x0,
        ifcs::Matrix{Float64},
        nsteps::Integer,
        nsteps_equil::Integer
    )

    command(lc.lmp, "reset_timestep 0")
    command(lc.lmp, "neigh_modify every 1 delay 0 check yes")

    scatter!(lc.lmp, "x", reinterpret(reshape, Float64, x0))
    seed = rand(1000:1_000_000)
    command(lc.lmp, "velocity all create $(Float64(T)) $(seed) dist gaussian mom yes")
    command(lc.lmp, "run 0 post no")
    
    # This makes energy available inside FixTIExternal
    command(lc.lmp, "thermo 1")
    command(lc.lmp, "thermo_style custom step pe")

    command(lc.lmp, "fix 1 all nve")

    dU_dλ = zeros(Float64, nsteps)
    FixTIExternal(ifcs, dU_dλ, lc.lmp, nsteps_equil, λ)

    command(lc.lmp, "run $(nsteps + nsteps_equil)")

    command(lc.lmp, "unfix 1")
    command(lc.lmp, "uncompute disp") # Defined in FixTIExternal
    command(lc.lmp, "unfix ti_julia") # Defined in FixTIExternal

    return sum(dU_dλ)/length(dU_dλ)

end

function LatticeDynamicsToolkit.ThermodynmicIntegration(
    ifc2::IFC2,
    sc::CrystalStructure,
    uc::CrystalStructure,
    T::Real,
    pot_cmds::Union{String, Vector{String}},
    nsteps::Integer,
    nsteps_equil::Integer,
    n_lambda::Integer;
    kwargs...
)

    if ifc2.na == length(sc)
        return ThermodynmicIntegration(
                ifc2,
                sc,
                T,
                pot_cmds,
                nsteps,
                nsteps_equil,
                n_lambda;
                kwargs...
            )
    else
        return ThermodynmicIntegration(
                remap(sc, uc, ifc2)[1],
                sc,
                T,
                pot_cmds,
                nsteps,
                nsteps_equil,
                n_lambda;
                kwargs...
            )
    end

end

# IGNORES POLAR INTERACTIONS
# returns free energy in eV
function LatticeDynamicsToolkit.ThermodynmicIntegration(
        ifc2::IFC2,
        sc::CrystalStructure,
        T::Real,
        pot_cmds::Union{String, Vector{String}},
        nsteps::Integer,
        nsteps_equil,
        n_lambda::Integer;
        kwargs...
    )

    @assert get(kwargs, :single_point, false) == false "Do not pass 'single_point' kwarg to TI, it is manually set to false."

    ifc2.na != length(sc) && error(ArgumentError("In thermodynamic integration, got IFC built on $(ifc2.na) cell, but supercell has $(length(sc)). Try remapping IFCs or passing the unitcell."))

    ifc2_lmp = forceconstant_2nd_HartreeBohr_to_eVA .* Matrix(ifc2, sc)
    
    x0_ang = bohr_to_A .* sc.x_cart

    λs, weights = gausslegendre(n_lambda)
    integrands = zeros(Float64, length(λs))

    lc = LAMMPSCalculator(sc, pot_cmds; single_point = false, kwargs...)
   
    p = Progress(length(λs); desc = "TI", dt = 0.5)
    for (i,λ) in enumerate(λs)
        integrands[i] = _TI_step(lc, λ, T, x0_ang, ifc2_lmp, nsteps, nsteps_equil)
        next!(p)
    end
    finish!(p)

    # free energy estimate in eV
    return sum(weights .* integrands)

end

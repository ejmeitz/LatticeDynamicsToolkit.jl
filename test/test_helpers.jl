# Shared test utilities and paths

const DATA_DIR = abspath(joinpath(@__DIR__, "..", "data"))

"""
    load_ifcs(basepath, ucposcar_path, T, N) -> (ifc2, ifc3, ifc4)

Load IFC2, IFC3, IFC4 from TDEP-style directory layout: `basepath/TK_NUC/infile.forceconstant*`
"""
function load_ifcs(basepath, ucposcar_path, T, N::Integer)
    subdir = joinpath(basepath, "$(Int(T))K_$(N)UC")
    ifc2 = read_ifc2(joinpath(subdir, "infile.forceconstant"), ucposcar_path)
    ifc3 = read_ifc3(joinpath(subdir, "infile.forceconstant_thirdorder"), ucposcar_path)
    ifc4 = read_ifc4(joinpath(subdir, "infile.forceconstant_fourthorder"), ucposcar_path)
    return ifc2, ifc3, ifc4
end

"""True if LAMMPSExt is loaded (LAMMPS available)."""
has_lammps() = Base.get_extension(LatticeDynamicsToolkit, :LAMMPSExt) !== nothing

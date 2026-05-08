# Struct-based amorphous fitting scaffold.
#
# This folder intentionally does not implement numerical solves yet. It only
# defines the public contracts and dispatch points for KrylovKit-backed
# fitting backends.

include("types.jl")
include("validation.jl")
include("contracts.jl")
include("pipeline.jl")
include("ifc2_common.jl")
include("ifc2_sparsecpu.jl")

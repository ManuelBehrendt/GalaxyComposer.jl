__precompile__(true)
module GalaxyComposer

using Parameters
using Printf
using QuadGK
using Roots
using SpecialFunctions
using Unitful, UnitfulAstro




export G, kB, mH
export globalset, generate, overview
export Unitful, UnitfulAstro


include("types.jl")
include("disk.jl")
include("dm.jl")
include("bulge.jl")
include("globalset.jl")
include("generator.jl")
include("overview.jl")


end # module GalaxyComposer

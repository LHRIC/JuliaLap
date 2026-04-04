module LonghornAPI

using Dates
using Genie
using Genie.Renderer.Json
using Genie.Requests
using Genie.Router
using JSON3
using UUIDs

include("kinematics/kin_solver.jl")
include("api/repo_index.jl")
include("api/actions.jl")
include("api/server.jl")

export build_catalog, build_action_registry, serve

end

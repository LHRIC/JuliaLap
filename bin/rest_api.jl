#!/usr/bin/env julia

using Pkg

project_root = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(project_root)

include(joinpath(project_root, "src", "LonghornAPI.jl"))
using .LonghornAPI

host = get(ENV, "LONGHORN_API_HOST", "127.0.0.1")
port = parse(Int, get(ENV, "LONGHORN_API_PORT", "8080"))

println("Starting LonghornAPI on $(host):$(port)")
LonghornAPI.serve(; host=host, port=port)

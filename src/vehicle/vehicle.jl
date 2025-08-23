using StaticArrays
include("rigid_body.jl")
include("simple_shock.jl")

struct Vehicle
    body::RigidBody
    front_shock::SimpleShock
    rear_shock::SimpleShock
end
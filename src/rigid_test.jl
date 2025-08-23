using LinearAlgebra
using StaticArrays
using DifferentialEquations
using Plots
using BenchmarkTools
include("vehicle/rigid_body.jl")
include("rotation/rotation_transformations.jl")


mass = 1.0
center_of_mass = SVector(0.0,0.0,0.0)
inertia_body = @SMatrix [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]
position0 = SVector(0.0,0.0,0.0)
linear_momentum0 = SVector(0.0,0.0,0.0)
rotation0 =SVector(1.0,0.0,0.0,0.0)
angular_momentum0 = SVector(0.0,1.0,0.01)

body = createRigidBody(mass, center_of_mass, inertia_body)

function unforced!(du::Vector{Float64}, u::Vector{Float64}, p::RigidBody, t)
    position = SVector{3, Float64}(u[1:3]) 
    linear_momentum = SVector{3, Float64}(u[4:6])
    rotation = MVector{4, Float64}(u[7:10])
    angular_momentum = SVector{3, Float64}(u[11:13])
    
    normalize!(rotation)
    u[7:10] = rotation
    R_mat = quaternion2matrix(rotation)
    angular_velocity = momentum2inertial(R_mat, p.inertia_body_inv, angular_momentum)

    du[1:3] = linear_momentum/p.mass
    du[4:6] .= 0.0
    du[7:10] = velocity2quaternion_dt(angular_velocity, rotation)
    du[11:13] .= 0.0
    nothing
end

u0 = [position0..., linear_momentum0..., rotation0..., angular_momentum0...]
tspan = (0.0, 50)

p = body

prob = ODEProblem(unforced!, u0, tspan, p)
sol = solve(prob, abstol=1e-12, reltol=1e-12)
eulers = mapslices(quaternion2euler, sol[7:10,:], dims=(1))
norms = mapslices(x-> 1.0-norm(x), sol[7:10,:], dims=(1))
plot(sol.t, eulers[1,:]*360/(2*π))
ylims!(0,180)
plot!(sol.t, eulers[2,:]*360/(2*π))
plot!(sol.t, eulers[3,:]*360/(2*π))

# plot(sol.t, norms[1,:])


# plot(sol.t[:], sol[7,:]/maximum(sol[7,:]))
# plot!(sol.t[:], sol[8,:]/maximum(sol[8,:]))
# plot!(sol.t[:], sol[9,:]/maximum(sol[9,:]))
# plot!(sol.t[:], sol[10,:]/maximum(sol[10,:]))




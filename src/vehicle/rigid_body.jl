using StaticArrays

struct RigidBody
    mass::Float64
    center_of_mass::SVector{3, Float64}

    inertia_body::SMatrix{3, 3, Float64}
    inertia_body_inv::SMatrix{3, 3, Float64}

end

function createRigidBody(mass::Float64, 
    center_of_mass::SVector{3, Float64}, 
    inertia_body::SMatrix{3, 3, Float64}) 

    inertia_body_inv = inv(inertia_body)
    rigid_body = RigidBody(mass, center_of_mass, inertia_body, inertia_body_inv)
    return rigid_body
end

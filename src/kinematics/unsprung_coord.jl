using LinearAlgebra
include("../rotation/rotation_transformations.jl")

function unsprung_coords(float_hdpts, contact_patch)
    LO = 1
    unsprung_coordinate_sys = unsprung_axes(float_hdpts)
    world_to_unsprung = inv(unsprung_coordinate_sys)
    contact_patch_local_wrld = contact_patch - float_hdpts[LO,:]
    contact_patch_local_uns = world_to_unsprung*contact_patch_local_wrld

    function coordinate_transform(float)
        unsprung_coordinates = unsprung_axes(float)
        # world_to_unsprung = inv(unsprung_coordinates)
        contact_patch_local_wrld = unsprung_coordinates * contact_patch_local_uns
        contact_patch_global_wrld = contact_patch_local_wrld + float[LO,:]
        T = [unsprung_coordinates contact_patch_global_wrld ; [0 0 0] 1]
        return T
    end

    return coordinate_transform
end

function unsprung_axes(float_hdpts)
    # Floating Hardpoint Indexes
    LO = 1 # Lower Outboard
    UO = 2 # Upper Outboard
    OT = 3 # Outboard Tie
    IT = 4 # Inboard Tie
    OP = 5 # Outboard Push/Pull
    IP = 6 # Inboard Push/Pull
    SO = 7 # Shock Outboard
    BP = 8 # ARB Bellcrank Pickup
    AP = 9 # ARB Arm Pickup

    axis_1 = float_hdpts[UO,:] - float_hdpts[LO,:]
    normalize!(axis_1)
    axis_2 = cross(axis_1, float_hdpts[IT,:] - float_hdpts[LO,:])
    normalize!(axis_2)
    axis_3 = cross(axis_1, axis_2)
    # println(axis_1)
    # println(axis_2)
    # println(axis_3)
    unsprung_coordinate_sys = [axis_1 axis_2 axis_3]
    return unsprung_coordinate_sys
end
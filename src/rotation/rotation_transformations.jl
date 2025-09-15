using StaticArrays
using LinearAlgebra

function quaternion2matrix(quat)
    w, x, y, z = quat
    R = SMatrix{3,3,Float64}(
        1 - 2*y*y - 2*z*z,  2*x*y - 2*z*w,      2*x*z + 2*y*w,
        2*x*y + 2*z*w,      1 - 2*x*x - 2*z*z,  2*y*z - 2*x*w,
        2*x*z - 2*y*w,      2*y*z + 2*x*w,      1 - 2*x*x - 2*y*y   
        )
    return R
end

function matrix_jac2velocity_jac(T,T_jac)
    type = eltype(T_jac)
    R = T[1:3,1:3]
    R_jac = T_jac[1:3,1:3,:]
    p_jac = T_jac[1:3,4,:]
    omega_jac = Array{type}(undef,(3,size(R_jac)[3]))
    for i = 1:size(R_jac)[3]
        R_dot = R_jac[:,:,i]
        omega_tensor = R_dot * inv(R)
        omega = [-omega_tensor[2,3], omega_tensor[1,3], -omega_tensor[1,2]]
        omega_jac[:,i] = omega
    end
    return vcat(p_jac, omega_jac)
end

function momentum2inertial(R, I_inv, L_inertial)
    return R*I_inv*transpose(R)*L_inertial
end

function quat_multiplication(q1, q2)
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2

    return @SVector [
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
        ]
end

function velocity2quaternion_dt(omega, quaternion) # Fixed reference frame
    return 0.5*quat_multiplication(SVector(0.0, omega...), quaternion)
end
 
function quaternion_dt2velocity(quaternion_dt, quaternion) # Fixed reference frame
    return 2.0*quat_multiplication(quaternion_dt, quaternion)
end

function quaternion2euler(quat)
    w, x, y, z = quat
    roll = atan(2*(w*x + y*z), 1 - 2*(x*x + y*y))
    pitch = asin(2*(w*y - z*x))
    yaw = atan(2*(w*y + x*y), 1 - 2*(y*y + z*z))
    return [roll, pitch, yaw]
end
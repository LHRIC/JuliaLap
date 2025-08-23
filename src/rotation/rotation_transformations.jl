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

function velocity2quaternion_dt(omega, quaternion)
    return 0.5*quat_multiplication(SVector(0.0, omega...), quaternion)
end

function quaternion2euler(quat)
    w, x, y, z = quat
    roll = atan(2*(w*x + y*z), 1 - 2*(x*x + y*y))
    pitch = asin(2*(w*y - z*x))
    yaw = atan(2*(w*y + x*y), 1 - 2*(y*y + z*z))
    return [roll, pitch, yaw]
end
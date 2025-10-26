#Things: 
#w_l (angular velocity for left wheel) 
#w_r (angular velocity for right wheel)
#N (Gear Ratio)
#T_d (drive train torque)
#T_l (torque to left wheel)
#T_r (torque to right wheel)

#Equations:
#w_l + w_r = w_d/N
#(T_l + T_r) * N = T_d
#T_d <= preload torque -> [locked diff]: w_l = w_r
#T_d > preload torque -> [unlocked diff]: the stuff in lsd.py

#Works up until the discs start slipping
#Vehicle State 1 is accelerating, 0 is accelerating

import DifferentialEquations as DE
function lsd(w, T)
    w_l, w_r, w_d = w
    T_l0, T_r0, T_d = T
    out1 = w_d/N - w_r - w_l
    out2 = T_d - N * (T_l + T_r)
    vehicle_state = 0
    if T_d > 0
        vehicle_state = 1 #Car is accelerating
    else
        vehicle_state = -1 #Car is deaccelerating
    end
    if abs(T_d) <= P
        if w_l != w_r
            #find the torque to each wheel that would make the speeds equal, static friction is 0.08 according to gadola, kinetic is 0.12
            w_new = (w_l + w_r)/2
            w_l = w_new
            w_r = w_new
        end
        out3 = w_l - w_r
    else
        delta_C_MAX = (2/3) * ((abs(C_m)/r_ramp) * (1/math.tan(math.radians(sigma_accel)))) * ((R_o^3 - R_i^3)/(R_o^2 - R_i^2)) * n * mu_c
        if T_d >= delta_C_MAX
            T_d = delta_C_MAX
            T_t = abs(T_d - delta_C_MAX)/2
        else
            T_t = 0
        end
        if vehicle_state == 1
            if left_vel <= right_vel
                T_l = (T_d + delta_C_MAX)/2 + T_t
                T_r = (T_d - delta_C_MAX)/2 + T_t
            else
                T_l = (T_d - delta_C_MAX)/2 + T_t
                T_r = (T_d + delta_C_MAX)/2 + T_t
            end
        else
            if left_vel <= right_vel
                T_l = -abs(T_d + delta_C_MAX)/2 - T_t
                T_r = -abs(T_d - delta_C_MAX)/2 - T_t
            else
                T_l = -abs(T_d - delta_C_MAX)/2 - T_t
                T_r = -abs(T_d + delta_C_MAX)/2 - T_t
            end
        end
        out3 = (T_l - T_l0) + (T_r - T_r0)
    end
    return [out1, out2, out3]
end


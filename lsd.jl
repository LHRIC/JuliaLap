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

import DifferentialEquations as DE
function lsd(w, T)
    w_l, w_r, w_d = u
    T_l0, T_r0, T_d = T
    out1 = w_d/N - w_r - w_l
    out2 = T_d - N * (T_l + T_r)
    if abs(T_d) <= P
        out3 = w_l - w_r
    else
        delta_C_MAX = (2/3) * ((abs(C_m)/r_ramp) * (1/math.tan(math.radians(sigma_accel)))) * ((R_o^3 - R_i^3)/(R_o^2 - R_i^2)) * n * mu_c
        if vehicle_state == 1
            if left_vel <= right_vel
                T_l = (T_d + delta_C_MAX)/2
                T_r = (T_d - delta_C_MAX)/2
            else
                T_l = (T_d - delta_C_MAX)/2
                T_r = (T_d + delta_C_MAX)/2
            end
        else
            if left_vel <= right_vel
                T_l = -abs(T_d + delta_C_MAX)/2
                T_r = -abs(T_d - delta_C_MAX)/2
            else
                T_l = -abs(T_d - delta_C_MAX)/2
                T_r = -abs(T_d + delta_C_MAX)/2
            end
        end
        out3 = (T_l - T_l0) + (T_r - T_r0)
    end
    return [out1, out2, out3]
end


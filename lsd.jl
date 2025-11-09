module LSD

using LinearAlgebra
using DifferentialEquations

# ──────────────────────────────────────────────────────────────
# Main function: limited-slip differential behavior
# ──────────────────────────────────────────────────────────────
"""
    lsd(w, T; N, P, C_m, r_ramp, sigma_accel, R_o, R_i, n, mu_c)

Computes the limited-slip differential (LSD) torque and wheel speed behavior.

Arguments:
- `w` = [w_l, w_r, w_d] (wheel and driveshaft angular velocities)
- `T` = [T_l0, T_r0, T_d] (torques to each wheel + driveshaft)
- Keyword parameters define differential geometry and friction constants.

Returns:
A vector `[out1, out2, out3]` representing:
1. Speed constraint residual  
2. Torque balance residual  
3. Wheel torque change term
"""
function lsd(w, T; N, P, C_m, r_ramp, sigma_accel, R_o, R_i, n, mu_c)
    w_l, w_r, w_d = w
    T_l0, T_r0, T_d = T

    # Constraints
    out1 = w_d / N - w_r - w_l
    out2 = T_d - N * (T_l0 + T_r0)

    # Determine vehicle state (accel or decel)
    vehicle_state = T_d > 0 ? 1 : -1

    # Initialize torque outputs
    T_l, T_r = T_l0, T_r0

    # Locked differential (below preload torque)
    if abs(T_d) <= P
        if w_l != w_r
            w_new = (w_l + w_r) / 2
            w_l = w_r = w_new
        end
        out3 = w_l - w_r

    # Unlocked differential (above preload torque)
    else
        delta_C_MAX = (2 / 3) * ((abs(C_m) / r_ramp) * (1 / tan(deg2rad(sigma_accel)))) *
                      ((R_o^3 - R_i^3) / (R_o^2 - R_i^2)) * n * mu_c

        T_d_eff = min(T_d, delta_C_MAX)
        T_t = abs(T_d - delta_C_MAX) / 2

        # Split torque by wheel speed difference
        if vehicle_state == 1
            if w_l <= w_r
                T_l = (T_d_eff + delta_C_MAX) / 2 + T_t
                T_r = (T_d_eff - delta_C_MAX) / 2 + T_t
            else
                T_l = (T_d_eff - delta_C_MAX) / 2 + T_t
                T_r = (T_d_eff + delta_C_MAX) / 2 + T_t
            end
        else
            if w_l <= w_r
                T_l = -abs(T_d_eff + delta_C_MAX) / 2 - T_t
                T_r = -abs(T_d_eff - delta_C_MAX) / 2 - T_t
            else
                T_l = -abs(T_d_eff - delta_C_MAX) / 2 - T_t
                T_r = -abs(T_d_eff + delta_C_MAX) / 2 - T_t
            end
        end

        out3 = (T_l - T_l0) + (T_r - T_r0)
    end

    return [out1, out2, out3]
end

# ──────────────────────────────────────────────────────────────
# Public exports
# ──────────────────────────────────────────────────────────────
export lsd

end # module

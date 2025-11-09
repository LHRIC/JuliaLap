module CompFunc

using ModelingToolkit
using DifferentialEquations

include("DampeningSpringsCode.jl")
include("new_torque_curves.jl")
include("lsd.jl")

using .DampeningSprings
using .TorqueCurves
using .LSD

# Symbolic placeholders
@parameters comp_interp(..)
@parameters reb_interp(..)
@parameters torque_fn(..)

# ────────────────────────────────────────────────
# Build the MTK system with symbolic placeholders
# ────────────────────────────────────────────────
function build_flexible_system(; damper_settings=(2.0,3.0),
                               torque_csv="38normalized.csv",
                               lsd_params=Dict(:N=>3.5, :P=>50.0, :C_m=>1.0,
                                               :r_ramp=>0.1, :sigma_accel=>30.0,
                                               :R_o=>0.1, :R_i=>0.05,
                                               :n=>4, :mu_c=>0.12))

    @parameters t
    @variables x(t) w_l(t) w_r(t) w_d(t)
    @variables T_l0(t) T_r0(t) T_d(t)
    D = Differential(t)

    eqs = [
        D(x) ~ comp_interp(x) + reb_interp(x),
        D(w_l) ~ 0,  # LSD will be applied numerically
        D(w_r) ~ 0,  # LSD will be applied numerically
        D(w_d) ~ torque_fn(w_d)
    ]

    @named sys = ODESystem(eqs, t,
                    [x, w_l, w_r, w_d, T_l0, T_r0, T_d],
                    [comp_interp, reb_interp, torque_fn])

    # Real functions
    comp_func, reb_func, _ = DampeningSprings.get_func(damper_settings)
    torque_real_fn = TorqueCurves.get_torque_interp(torque_csv)

    # Numeric LSD function
    lsd_real_fn = (w,T) -> LSD.lsd(w,T;
                                   N=lsd_params[:N], P=lsd_params[:P], C_m=lsd_params[:C_m],
                                   r_ramp=lsd_params[:r_ramp], sigma_accel=lsd_params[:sigma_accel],
                                   R_o=lsd_params[:R_o], R_i=lsd_params[:R_i],
                                   n=lsd_params[:n], mu_c=lsd_params[:mu_c])

    return sys, Dict(
        comp_interp => comp_func,
        reb_interp  => reb_func,
        torque_fn   => torque_real_fn,
        :lsd_fn     => lsd_real_fn
    )
end

end

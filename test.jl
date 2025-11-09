using DifferentialEquations
include("compfunc.jl")
using .CompFunc

# ─────────────────────────────────────
# Build the symbolic system and functions
# ─────────────────────────────────────
sys, real_funcs = CompFunc.build_flexible_system()
lsd_fn = real_funcs[:lsd_fn]

# Extract functions directly by name
comp_func  = real_funcs[CompFunc.comp_interp]
reb_func   = real_funcs[CompFunc.reb_interp]
torque_func = real_funcs[CompFunc.torque_fn]

# ─────────────────────────────────────
# Define initial conditions and ODE setup
# ─────────────────────────────────────
u0 = [200, 180, 190, 4500, 18, 22, 40]   # x(0-250 mm/s), w_l, w_r, w_d(4000-13000rpm), T_l0, T_r0, T_d (22-69Nm)
tspan = (0.0, 1.0) #timespan

function f!(du, u, p, t)
    x, w_l, w_r, w_d, T_l0, T_r0, T_d = u

    du[1] = comp_func(x) + reb_func(x)
    du[2] = 0.0 #Placeholders so solver can know size of the system?
    du[3] = 0.0
    du[4] = torque_func(w_d)

    # LSD numerical output
    lsd_out = lsd_fn([((2*pi)/60)*w_l,((2*pi)/60)*w_r, ((2*pi)/60)*w_d], [T_l0, T_r0, T_d]) #Converted to rad/sec
    du[5] = lsd_out[1]
    du[6] = lsd_out[2]
    du[7] = lsd_out[3]
end

# ─────────────────────────────────────
# Solve ODE
# ─────────────────────────────────────
prob = ODEProblem(f!, u0, tspan)
sol = solve(prob)

println(sol.t)
println(sol.u)


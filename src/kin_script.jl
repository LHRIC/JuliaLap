using Plots
include("kinematics/kinematics.jl")
using .Kinematics

# @time sol = gen_kin_model("src/parameters/HDPT_Export.xlsx")
# @time sol = gen_kin_model("src/parameters/HDPT_Export.xlsx")
steer_range = range(-20, 20, 20)
f_shock_range = range(-10, 10, 20)
r_shock_range = range(-10, 10, 20)
fl_T_itp, fl_jac_itp, rl_T_itp, rl_jac_itp = gen_kin_models("src/parameters/HDPT_Export.xlsx", steer_range, f_shock_range, r_shock_range)
# @time @profview sol = gen_kin_models("src/parameters/HDPT_Export.xlsx")
# @benchmark sol = gen_kin_models("src/parameters/HDPT_Export.xlsx")
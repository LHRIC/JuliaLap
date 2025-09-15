using Profile
using BenchmarkTools
include("kinematics/kin_solver.jl")
file_name = "src/parameters/HDPT_Export.xlsx"
fl_range = "B8:E25"
rl_range = "B33:E50"
fl_dict, rl_dict = excel2dict(file_name, fl_range, rl_range)
fl_array = dict2cvec(fl_dict, "FL")
fun, ctrl = gen_corner(fl_array)
# @benchmark fun, ctrl = gen_corner(fl_array)
ctrl = [ctrl[1] + 10, ctrl[2]]

# @benchmark pos = fun(ctrl)
# jac = FiniteDiff.finite_difference_jacobian(fl_fun,ctrl1,absstep=1e-12)
# @benchmark jac = ForwardDiff.jacobian(fun,ctrl)
@benchmark jac = FiniteDiff.finite_difference_jacobian(fun,ctrl)

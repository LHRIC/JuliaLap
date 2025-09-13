using Profile
include("kinematics/kin_solver.jl")

# @time sol = gen_kin_model("src/parameters/HDPT_Export.xlsx")
# @time sol = gen_kin_model("src/parameters/HDPT_Export.xlsx")
pos1, jac1 = gen_kin_models("src/parameters/HDPT_Export-9-13.xlsx")
# @time @profview sol = gen_kin_models("src/parameters/HDPT_Export.xlsx")

# @benchmark sol = gen_kin_models("src/parameters/HDPT_Export.xlsx")
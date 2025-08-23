using Profile
include("kinematics/kin_solver.jl")

# @time sol = gen_kin_model("src/parameters/HDPT_Export.xlsx")
# @time sol = gen_kin_model("src/parameters/HDPT_Export.xlsx")
@time sol = gen_kin_model("src/parameters/HDPT_Export.xlsx")
@profview sol = gen_kin_model("src/parameters/HDPT_Export.xlsx")

# @benchmark sol = gen_kin_model("src/parameters/HDPT_Export.xlsx")
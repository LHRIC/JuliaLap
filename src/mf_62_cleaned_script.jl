include("tires/mf_62_cleaned.jl")
include("tires/parse_tire.jl")
using Plots
using ModelingToolkit
plotlyjs()

# model = parse_tir("src/parameters/Round_8_Hoosier_R25B_16x7p5_10_on_7in_12psi_PAC02_UM2.tir")
# model = parse_tir("src/parameters/R20_16x7p5_10_on_7in_12psi_PAC2002.tir")
# model = parse_tir("src/parameters/FSAE_Defaults.tir")
model = parse_tir("src/parameters/R20_16x7p5_10_on_7in_12psi_PAC2002.tir")
# model = parse_tir("src/parameters/Dynamics copy.tir")


model["LMUV"] = 1
data = Dict{String, Any}("vcx" => 1, "vc" => 1)

alpha = -deg2rad(24):deg2rad(0.01):deg2rad(24)
kappa = -1.2:0.01:1.2
# kappa = -0.2:0.01:0.2
gamma = -deg2rad(5):deg2rad(2):deg2rad(5)

fz = 2092.99
out = []
fzs = 1000:100:2092.99
# fzs = 1250
g = 0
# plot(kappa, out, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])
# println(kappa, out)
k = 1
fz = 2092.99
p = model
# MF62.base(1, 1, fz, 0, k, g, 0)
# lfz0, p_i, p_io, lmux, lmuy, lmuv, r0, g, vcx, vc, fz0 = MF62.params
# fun = MF62.fx0
@mtkcompile sys = OptimizationSystem(MF62.fx0, [MF62.params...], [])
# params_dict = Dict(
#     lfz0 => p["LFZO"],
#     fz0  => p["FNOMIN"],
#     lmux => p["LMUX"],
#     lmuy => p["LMUY"],
#     lmuv => p["LMUV"],
#     r0   => p["UNLOADED_RADIUS"],
#     g    => p["GRAVITY"],
#     vcx  => data["vcx"],
#     vc   => data["vc"],
# )
# val = substitute(fun, params_dict)
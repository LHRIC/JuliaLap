include("tires/mf_62_sym.jl")
include("tires/parse_tire.jl")
using Plots
using ModelingToolkit
plotlyjs()

# model = parse_tir("src/parameters/Round_8_Hoosier_R25B_16x7p5_10_on_7in_12psi_PAC02_UM2.tir")
model = parse_tir("src/parameters/R20_16x7p5_10_on_7in_12psi_PAC2002.tir")
# model = parse_tir("src/parameters/FSAE_Defaults.tir")
# model = parse_tir("src/parameters/R20_16x7p5_10_on_7in_12psi_PAC2002.tir")
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
# @mtkcompile sys = OptimizationSystem(MF62.fx0, [MF62.params...], [])

param_map = [
    MF62.lfz0 => p["LFZO"],
    MF62.fz0  => p["FNOMIN"],
    MF62.lmux => p["LMUX"],
    MF62.lmuy => p["LMUY"],
    MF62.lmuv => p["LMUV"],
    MF62.r0   => p["UNLOADED_RADIUS"],
    MF62.g    => p["GRAVITY"],
    MF62.vcx  => 1,
    MF62.vc   => 1,

    MF62.fz0 => p["FNOMIN"],

    MF62.pcx1 => p["PCX1"],
    MF62.pdx1 => p["PDX1"],
    MF62.pdx2 => p["PDX2"],
    MF62.pdx3 => p["PDX3"],
    MF62.pex1 => p["PEX1"],
    MF62.pex2 => p["PEX2"],
    MF62.pex3 => p["PEX3"],
    MF62.pex4 => p["PEX4"],
    MF62.pkx1 => p["PKX1"],
    MF62.pkx2 => p["PKX2"],
    MF62.pkx3 => p["PKX3"],
    MF62.phx1 => p["PHX1"],
    MF62.phx2 => p["PHX2"],
    MF62.pvx1 => p["PVX1"],
    MF62.pvx2 => p["PVX2"],
    MF62.ppx1 => p["PPX1"],
    MF62.ppx2 => p["PPX2"],
    MF62.ppx3 => p["PPX3"],
    MF62.ppx4 => p["PPX4"],

    MF62.lcx => p["LCX"],
    MF62.lhx => p["LHX"],
    MF62.lex => p["LEX"],
    MF62.lvx => p["LVX"],

    MF62.friction_scaling_x => p["friction_scaling_x"],

    MF62.pcy1 => p["PCY1"],
    MF62.pdy1 => p["PDY1"],
    MF62.pdy2 => p["PDY2"],
    MF62.pdy3 => p["PDY3"],
    MF62.pey1 => p["PEY1"],
    MF62.pey2 => p["PEY2"],
    MF62.pey3 => p["PEY3"],
    MF62.pey4 => p["PEY4"],
    MF62.pey5 => p["PEY5"],
    MF62.pky1 => p["PKY1"],
    MF62.pky2 => p["PKY2"],
    MF62.pky3 => p["PKY3"],
    MF62.pky4 => p["PKY4"],
    MF62.pky5 => p["PKY5"],
    MF62.pky6 => p["PKY6"],
    MF62.pky7 => p["PKY7"],
    MF62.phy1 => p["PHY1"],
    MF62.phy2 => p["PHY2"],
    MF62.pvy1 => p["PVY1"],
    MF62.pvy2 => p["PVY2"],
    MF62.pvy3 => p["PVY3"],
    MF62.pvy4 => p["PVY4"],
    MF62.ppy1 => p["PPY1"],
    MF62.ppy2 => p["PPY2"],
    MF62.ppy3 => p["PPY3"],
    MF62.ppy4 => p["PPY4"],
    MF62.ppy5 => p["PPY5"],

    MF62.lcy => p["LCY"],
    MF62.lhy => p["LHY"],
    MF62.ley => p["LEY"],
    MF62.lvy => p["LVY"],
    MF62.lky => p["LKY"],
    MF62.lkyg => p["LKYC"],

    MF62.phy3 => 0.144919,

    MF62.qhz1 => p["QHZ1"],
    MF62.qhz2 => p["QHZ2"],
    MF62.qhz3 => p["QHZ3"],
    MF62.qhz4 => p["QHZ4"],
    MF62.qbz1 => p["QBZ1"],
    MF62.qbz2 => p["QBZ2"],
    MF62.qbz3 => p["QBZ3"],
    MF62.qbz5 => p["QBZ5"],
    MF62.qbz6 => p["QBZ6"],
    MF62.qbz9 => p["QBZ9"],
    MF62.qbz10 => p["QBZ10"],
    MF62.qcz1 => p["QCZ1"],
    MF62.qdz1 => p["QDZ1"],
    MF62.qdz2 => p["QDZ2"],
    MF62.qdz3 => p["QDZ3"],
    MF62.qdz4 => p["QDZ4"],
    MF62.qdz6 => p["QDZ6"],
    MF62.qdz7 => p["QDZ7"],
    MF62.qdz8 => p["QDZ8"],
    MF62.qdz9 => p["QDZ9"],
    MF62.qdz10 => p["QDZ10"],
    MF62.qdz11 => p["QDZ11"],
    MF62.ppz1 => p["PPZ1"],
    MF62.ppz2 => p["PPZ2"],
    MF62.qez1 => p["QEZ1"],
    MF62.qez2 => p["QEZ2"],
    MF62.qez3 => p["QEZ3"],
    MF62.qez4 => p["QEZ4"],
    MF62.qez5 => p["QEZ5"],

    MF62.lky => p["LKY"],
    MF62.ltr => p["LTR"],
    MF62.lres => p["LRES"],
    MF62.lkzc => p["LKZC"],   

    MF62.rbx1 => p["RBX1"],
    MF62.rbx2 => p["RBX2"],
    MF62.rbx3 => p["RBX3"],
    MF62.rcx1 => p["RCX1"],
    MF62.rex1 => p["REX1"],
    MF62.rex2 => p["REX2"],
    MF62.rhx1 => p["RHX1"],

    MF62.lxal => p["LXAL"],

    MF62.rby1 => p["RBY1"],
    MF62.rby2 => p["RBY2"],
    MF62.rby3 => p["RBY3"],
    MF62.rby4 => p["RBY4"],
    MF62.rcy1 => p["RCY1"],
    MF62.rey1 => p["REY1"],
    MF62.rey2 => p["REY2"],
    MF62.rhy1 => p["RHY1"],
    MF62.rhy2 => p["RHY2"],
    MF62.rvy1 => p["RVY1"],
    MF62.rvy2 => p["RVY2"],
    MF62.rvy3 => p["RVY3"],
    MF62.rvy4 => p["RVY4"],
    MF62.rvy5 => p["RVY5"],
    MF62.rvy6 => p["RVY6"],

    MF62.lyk => p["LYKA"],
    MF62.lvyk => p["LVYKA"],

    MF62.qsx1 => p["QSX1"],
    MF62.qsx2 => p["QSX2"],
    MF62.qsx3 => p["QSX3"],
    MF62.qsx4 => p["QSX4"],
    MF62.qsx5 => p["QSX5"],
    MF62.qsx6 => p["QSX6"],
    MF62.qsx7 => p["QSX7"],
    MF62.qsx8 => p["QSX8"],
    MF62.qsx9 => p["QSX9"],
    MF62.qsx10 => p["QSX10"],
    MF62.qsx11 => p["QSX11"],
    MF62.ppmx1 => p["PPMX1"],

    MF62.lvmx => p["LVMX"],
    MF62.lmx => p["LMX"],

    MF62.qsy1 => p["QSY1"],
    MF62.qsy2 => p["QSY2"],
    MF62.qsy3 => p["QSY3"],
    MF62.qsy4 => p["QSY4"],
    MF62.qsy5 => p["QSY5"],
    MF62.qsy6 => p["QSY6"],
    MF62.qsy7 => p["QSY7"],
    MF62.qsy8 => p["QSY8"],

    MF62.lmy => p["LMY"],

    MF62.ssz1 => p["SSZ1"],
    MF62.ssz2 => p["SSZ2"],
    MF62.ssz3 => p["SSZ3"],
    MF62.ssz4 => p["SSZ4"],

]
fx0_fun = substitute(MF62.k_xk, param_map)
fx0_build = eval(build_function(fx0_fun, MF62.fz, MF62.gamma, MF62.kappa, MF62.vx, MF62.vy))
fz_in = range(0,2000,100)
plot(fz_in,fx0_build.(fz_in,0,1.2,1,1))
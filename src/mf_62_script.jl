include("tires/mf_62.jl")
include("tires/parse_tire.jl")
using Plots

model = MF62.MF62Model(parse_tir("src/parameters/FSAE_Defaults.tir"))
kappa = -0.5:0.001:0.5
gamma = -pi/4:pi/8:pi/4

group = []
for i in gamma
    local fx = []
    for j in kappa
        push!(fx, MF62.fx(model, 400.0, j, i))
    end
    push!(group, fx)
end

display(plot(kappa, group, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"]))

group = []
for i in gamma
    local fy = []
    for j in kappa
        push!(fy, MF62.fy(model, 400.0, j, i))
    end
    push!(group, fy)
end

display(plot(kappa, group, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"]))
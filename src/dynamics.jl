include("tires/mf_62.jl")
include("tires/parse_tire.jl")
using Plots
plotlyjs()

model = parse_tir("src/parameters/Dynamics.tir")
# print("model", model)
# default = parse_tir("src/parameters/FSAE_Defaults.tir")
# print("default", default)

alpha = -deg2rad(4):deg2rad(0.01):deg2rad(2)
kappa = -1.2:0.001:1.2
gamma = -deg2rad(3):deg2rad(2):deg2rad(3)

# groupx = []

# for i in gamma
#     local fx = []
#     for j in kappa
#         push!(fx, MF62.fx0(model, 800.0, j, i))
#     end
#     push!(groupx, fx)
# end

# labelsx = reshape([string(round(rad2deg(g)), "°") for g in gamma], 1, :)

# plot(kappa, groupx, 
#     xlabel = "Kappa (slip ratio)",
#     ylabel = "Fx0 (longitudinal force)",
#     title = "Model 800 N 16 in tire", 
#     label = labelsx) #["-10°" "-9°" "-8°" "-7°" "-6°" "-5°" "-4°" "-3°" "-2°" "-1°" "0°" "1°" "2°" "3°" "4°" "5°" "6°" "7°" "8°" "9°" "10°"])



group = []

for i in gamma
    local fy = []
    for j in alpha
        push!(fy, MF62.fy0(model, 1500.0, j, i))
    end
    push!(group, fy)
end

labels = reshape([string(round(rad2deg(g)), "°") for g in gamma], 1, :)

plot(rad2deg.(alpha), group, 
    xlabel = "Alpha (slip angle)",
    ylabel = "Fy0 (lateral force)",
    title = "16x7.5-10_R20_8_HP 1500 N", 
    label = labels)
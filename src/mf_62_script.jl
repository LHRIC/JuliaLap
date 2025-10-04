include("tires/mf_62.jl")
include("tires/parse_tire.jl")
using Plots
plotlyjs()

model = parse_tir("src/parameters/FSAE_Defaults.tir")
alpha = -deg2rad(20):0.005:deg2rad(20)
kappa = -2:0.005:2
gamma = -pi/4:pi/8:pi/4



data = [MF62.fx(model, 400, a, k, gamma[1]) for k in kappa, a in alpha]
plot(alpha, kappa, data;
    seriestype = :surface,
    xlabel = "Alpha (slip angle)",
    ylabel = "Kappa (slip ratio)",
    zlabel = "Fy (Lateral Force)",
    title = "Magic Formula 62 Tire Model",
    legend = false,
    markersize = 4,
    markerstrokewidth = 0,
    markercolor = :blue
)

# for i in gamma 
#     local at = []
#     for j in kappa 
#         push!(at, MF62.at(model, 400.0, j, i))
#     end
#     push!(group, at)
# end
# plot(kappa, group, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])

# for i in gamma
#     local fx0 = []
#     for j in kappa
#         push!(fx0, MF62.fx0(model, 400.0, j, i))
#     end
#     push!(group, fx0)
# end
# plot(kappa, group, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])

# group = []
# for i in gamma
#     local fy0 = []
#     for j in kappa
#         push!(fy0, MF62.fy0(model, 400.0, j, i))
#     end
#     push!(group, fy0)
# end
# plot(kappa, group, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])
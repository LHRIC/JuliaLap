include("tires/mf_62.jl")
include("tires/parse_tire.jl")
using Plots
plotlyjs()

model = parse_tir("src/parameters/FSAE_Defaults.tir")
alpha = -0.5:0.005:0.5
kappa = -0.5:0.005:0.5
gamma = -pi/4:pi/8:pi/4

data = []
# Manually change gamma[n] for the different plots
for k in kappa
    for a in alpha
        push!(data, MF62.fx(model, 400, a, k, gamma[5]))
    end
end
data_mat = reshape(data, (length(alpha), length(kappa)))
# Extract coordinates for plotting

plot(alpha, kappa, data_mat;
    seriestype = :surface,
    xlabel = "Kappa (slip ratio)",
    ylabel = "Alpha (slip angle)",
    zlabel = "Fx (Longitudinal Force)",
    title = "Magic Formula 62 Tire Model",
    legend = false,
    markersize = 4,
    markerstrokewidth = 0,
    markercolor = :blue,
)

# for i in gamma 
#     local at = []
#     for j in kappa 
#         push!(at, MF62.at(model, 400.0, j, i))
#     end
#     push!(group, at)
# end

# display(plot(kappa, group, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"]))

# for i in gamma
#     local fx = []
#     for j in kappa
#         push!(fx, MF62.fx(model, 400.0, j, i))
#     end
#     push!(group, fx)
# end

# display(plot(kappa, group, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"]))

# group = []
# for i in gamma
#     local fy = []
#     for j in kappa
#         push!(fy, MF62.fy(model, 400.0, j, i))
#     end
#     push!(group, fy)
# end

# display(plot(kappa, group, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"]))
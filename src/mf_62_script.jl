include("tires/mf_62.jl")
include("tires/parse_tire.jl")
using Plots

model = parse_tir("src/parameters/FSAE_Defaults.tir")
alpha = -0.5:0.001:0.5
kappa = -0.5:0.001:0.5
gamma = -pi/4:pi/8:pi/4

# group = []

data = []  # collect results as tuples (kappa, alpha, fx)

for i in gamma
    for j in kappa
        for a in alpha
            push!(data, (j, a, MF62.fx(model, 400, j, a, i)))
        end
    end
end

# Extract coordinates for plotting
kappa_vals = [d[1] for d in data]
alpha_vals = [d[2] for d in data]
fx_vals = [d[3] for d in data]

plot3d(kappa_vals, alpha_vals, fx_vals,
    seriestype = :scatter,
    xlabel = "Kappa (slip ratio)",
    ylabel = "Alpha (slip angle)",
    zlabel = "Fx (Longitudinal Force)",
    title = "Magic Formula 62 Tire Model",
    legend = false,
    marker = (:circle, 4, 0.6, :blue)
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
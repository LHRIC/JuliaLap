include("tires/mf_62_cleaned.jl")
include("tires/parse_tire.jl")
using Plots
plotlyjs()

model = parse_tir("src/parameters/FSAE_Defaults.tir")
model["LMUV"] = 1
data = Dict{String, Any}("vcx" => 1, "vc" => 1)

alpha = -deg2rad(24):deg2rad(0.01):deg2rad(24)
kappa = -1.2:0.01:1.2
gamma = -deg2rad(5):deg2rad(2):deg2rad(5)

out = []
for g in gamma
    local at = []
    for a in alpha
        MF62.base(model, data, 1, 1, 800, a, 0, g, 0)
        push!(at, MF62.at(model, data, 800.0, a, 0, g))
    end
    push!(out, at)
end
plot(alpha, out, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])

# out = []
# for g in gamma
#     local rrm = []
#     for k in kappa
#         MF62.base(model, data, 1, 1, 800, 0, k, g, 0)
#         push!(rrm, MF62.rrm(model, data, 1, 800.0, 0, k, g))
#     end
#     push!(out, rrm)
# end
# plot(kappa, out, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])

# out = []
# for g in gamma
#     local oc = []
#     for a in alpha
#         MF62.base(model, data, 1, 1, 800, a, 0, g, 0)
#         push!(oc, MF62.oc(model, data, 800.0, a, 0, g))
#     end
#     push!(out, oc)
# end
# plot(alpha, out, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])

# out = zeros(length(alpha), length(kappa))
# for a in eachindex(alpha)
#     for k in eachindex(kappa)
#         MF62.base(model, data, 1, 1, 400, alpha[a], kappa[k], gamma[1], 0)
#         out[a, k] = MF62.fy(model, data, 400.0, alpha[a], kappa[k], gamma[1])
#     end
# end
# plot(alpha, kappa, out';
#     seriestype = :surface,
#     xlabel = "Alpha (slip angle)",
#     ylabel = "Kappa (slip ratio)",
#     zlabel = "Fy (Lateral Force)",
#     title = "Magic Formula 62 Tire Model",
#     legend = false,
#     markersize = 4,
#     markerstrokewidth = 0,
#     markercolor = :blue
# )

# out = zeros(length(alpha), length(kappa))
# for a in eachindex(alpha)
#     for k in eachindex(kappa)
#         MF62.base(model, data, 1, 1, 400, alpha[a], kappa[k], gamma[1], 0)
#         out[a, k] = MF62.fx(model, data, 400.0, alpha[a], kappa[k], gamma[1])
#     end
# end
# plot(alpha, kappa, out';
#     seriestype = :surface,
#     xlabel = "Alpha (slip angle)",
#     ylabel = "Kappa (slip ratio)",
#     zlabel = "Fx (Longitudinal Force)",
#     title = "Magic Formula 62 Tire Model",
#     legend = false,
#     markersize = 4,
#     markerstrokewidth = 0,
#     markercolor = :blue
# )

# out = []
# for g in gamma
#     local at0 = []
#     for a in alpha
#         MF62.base(model, data, 1, 1, 400, a, 0, g, 0)
#         push!(at0, MF62.at0(model, data, 400.0, a, g))
#     end
#     push!(out, at0)
# end
# plot(alpha, out, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])

# out = []
# for g in gamma
#     local fy0 = []
#     for a in alpha
#         MF62.base(model, data, 1, 1, 400, a, 0, g, 0)
#         push!(fy0, MF62.fy0(model, data, 400.0, a, g))
#     end
#     push!(out, fy0)
# end
# plot(alpha, out, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])

# out = []
# for g in gamma
#     local fx0 = []
#     for k in kappa
#         MF62.base(model, data, 1, 1, 400, 0, k, g, 0)
#         push!(fx0, MF62.fx0(model, data, 400.0, k, g))
#     end
#     push!(out, fx0)
# end
# plot(kappa, out, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])
include("tires/mf_62_cleaned.jl")
include("tires/parse_tire.jl")
using Plots
# plotlyjs()

model = parse_tir("src/parameters/FSAE_Defaults.tir")
model["LMUV"] = 1
data = Dict{String, Any}("vcx" => 1, "vc" => 1)

alpha = -deg2rad(12):deg2rad(0.01):deg2rad(12)
kappa = -1.2:0.001:1.2
gamma = -deg2rad(5):deg2rad(2):deg2rad(5)


group = []
for g in gamma
    local fx0 = []
    for k in kappa
        MF62.base(model, data, 1, 1, 800, 0, k, g, 0)
        push!(fx0, MF62.fx0(model, 400.0, j, i))
    end
    push!(group, fx0)
end
plot(kappa, group, label=["-pi/4" "-pi/8" 0 "pi/8" "pi/4"])


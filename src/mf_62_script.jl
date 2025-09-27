include("tires/mf_62.jl")
include("tires/parse_tire.jl")

model = MF62.MF62Model(parse_tir("src/parameters/FSAE_Defaults.tir"))
println(MF62.fx(model, 400.0, 0.0, 0.0))
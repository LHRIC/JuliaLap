
include(joinpath(@__DIR__, "tires", "parse_ttc.jl"))

df, dict = parse_ttc(joinpath(@__DIR__, "parameters", "B2356raw2.dat"))
df2 = unitful_to_float(df)

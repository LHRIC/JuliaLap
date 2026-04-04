
include(joinpath(@__DIR__, "tires", "parse_ttc.jl"))

filepath = joinpath(@__DIR__, "parameters", "B2356raw2.dat")
df, dict = parse_ttc(filepath)
df2 = unitful_to_float(df)
plot_all_vs_effective_time(df; dataset_label = filepath)

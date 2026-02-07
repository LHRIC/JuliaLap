
include("tires/parse_ttc.jl")

df, dict = parse_ttc("src/parameters/B2356raw2.dat")
df2 = unitful_to_float(df)
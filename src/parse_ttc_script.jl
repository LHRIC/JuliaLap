include("tires/parse_ttc.jl")

# Get the absolute path to your home directory and join it with the folder/file
path = joinpath(homedir(), "lhr", "ttcdata", "B2356raw2.dat")
parse_ttc(path)

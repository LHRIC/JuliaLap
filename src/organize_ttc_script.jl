#find ~/lhr/ttcdata -name "*:Zone.Identifier" -type f -delete

include("tires/organize_ttc.jl")
 
path = joinpath(homedir(), "lhr", "ttcdata")
organize_ttc(path)

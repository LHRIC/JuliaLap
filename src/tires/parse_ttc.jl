
using MAT
using CSV
using DataFrames

# Outputs: (DataFrame df, Dict dict)
# df - DataFrame where each column represents a value and each row represents a trial
# dict - Dictionary consisting of all unused metadata
function parse_ttc(filepath::String, wanted_cols = ["TSTO", "RE", "P", "AMBTMP", 
    "TSTC", "FY", "V", "NFX", "SA", "RST", "N", "ET", "SL", "TSTI", "MX", 
    "FZ", "RUN", "RL", "SR", "MZ", "NFY", "FX", "IA"])

    filetype = split(filepath, ".")[lastindex(split(filepath, "."))]
    df = nothing
    dict = nothing

    if (filetype == "mat")
        dict = matread(filepath)
        
        values = []
        for key in wanted_cols
            push!(values, vec(dict[key]))
            delete!(dict, key)
        end

        df = DataFrame(wanted_cols .=> values)
    elseif (filetype == "dat")

        meta = nothing
        keys = nothing
        units = nothing
        open(filepath, "r") do io
            meta = readline(io)
            keys = readline(io)
            units = readline(io)
        end

        parts = split(meta, ';')
        dict = Dict{String, Any}()
        for part in parts
            part = strip(part)
            isempty(part) && continue

            # Case 1: "Key: Value"
            if occursin(":", part)
                key, value = split(part, ":", limit=2)
                dict[strip(key)] = strip(value)
            
            # Case 2: "Key=Value"
            elseif occursin("=", part)
                key, value = split(part, "=", limit=2)
                dict[strip(key)] = strip(value)

            # Case 3: "ISO False" or "ID YXXXX"
            else
                tokens = split(part)
                if length(tokens) ≥ 2
                    key = tokens[1]
                    value = join(tokens[2:end], " ")
                    dict[key] = value
                end
            end
        end

        keys = split(keys, "\t")
        units = split(units, "\t")
        dict["units"] = Dict(zip(keys, units))
        
        df = CSV.read(filepath, DataFrame; header=2, skipto=4, delim='\t')
    end

    return df, dict
end

# for a list of ttc file
# df - a concatenated list of dataframes; added column for file ID
# dict - fileID => metadata dictionary
function parse_ttc_list(file_list, wanted_cols = ["TSTO", "RE", "P", "AMBTMP", 
    "TSTC", "FY", "V", "NFX", "SA", "RST", "N", "ET", "SL", "TSTI", "MX", 
    "FZ", "RUN", "RL", "SR", "MZ", "NFY", "FX", "IA"])

    df = DataFrame()
    dict = Dict{String, Dict}()

    for filepath in file_list
        temp_df, temp_dict = parse_ttc(filepath, wanted_cols)

        id = temp_dict["ID"]

        temp_df[!, :ID] .= id
        append!(df, temp_df)
        dict[id] = temp_dict
    end

    return df, dict
end
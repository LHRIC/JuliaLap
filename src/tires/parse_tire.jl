
function parse_tir(filepath::String)
    parsed_data = Dict{String, Any}()

    open(filepath, "r") do file
        for line in eachline(file)
            if occursin("=", line)
                # Split key-value pair on "="
                parts = split(line, "=", limit=2)
                key = strip(parts[1])

                # Remove everything after "$"
                value = strip(split(parts[2], "\$", limit=2)[1])

                # Try to convert value to Int or Float
                if tryparse(Int, value) !== nothing
                    value = parse(Int, value)
                elseif tryparse(Float64, value) !== nothing
                    value = parse(Float64, value)
                else
                    value = strip(value, [''', '"'])
                end

                parsed_data[key] = value
            end
        end
    end

    return parsed_data
end
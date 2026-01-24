
function write_tir(filepath::String, data::Dict{String,Dict{String,Any}})
    open(filepath, "w") do io
        for (section, entries) in data
            println(io, "[$section]")

            for (key, value) in entries
                val_str = if value isa String
                    "'$value'"
                elseif value isa Int
                    string(value)
                elseif value isa Float64
                    @sprintf("%.6g", value)
                else
                    string(value)
                end

                @printf(io, "%-25s =    %s\n", key, val_str)
            end

            println(io)  # blank line between sections
        end
    end
end

using Printf

function write_tir(template_path::String, output_path::String, values::Dict{String,Any})

    open(template_path, "r") do infile
        open(output_path, "w") do outfile

            default_value = 0

            for line in eachline(infile)
                stripped = strip(line)

                # Sets the default value
                if startswith("[", stripped) && occursin("COEFFICIENT", line)
                    default_value = 1
                else
                    default_value = 0
                end

                # Only touch lines with "=" that are not comments
                if occursin("=", line) && !startswith(stripped, "\$") && !startswith(stripped, "!")
                    parts = split(line, "=", limit = 2)
                    key = strip(parts[1])

                    # Preserve comment if it exists
                    value_part = parts[2]
                    split_comment = split(value_part, "\$", limit = 2)
                    comment = length(split_comment) == 2 ? " \$" * split_comment[2] : ""

                    # Determine value
                    val = haskey(values, key) ? values[key] : default_value

                    # Format value
                    val_str = if val isa String
                        "'$val'"
                    elseif val isa Int
                        string(val)
                    elseif val isa Float64
                        @sprintf("%.6g", val)
                    else
                        string(val)
                    end

                    @printf(outfile, "%-25s =    %s%s\n", key, val_str, comment)
                else
                    println(outfile, line)
                end
            end
        end
    end
end
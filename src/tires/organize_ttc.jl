using Base.Filesystem

# Include the parser logic (assumed to be in the same folder)
include("parse_ttc.jl")

"""
    organize_ttc(source_dir::String)

Parses .dat and .mat files in `source_dir`. 
- If "Tire_Name" contains "Hoosier", it moves to a folder named after the tire.
- If not, it moves the file to an "Other" folder.
"""
function organize_ttc(source_dir::String)
    if !isdir(source_dir)
        println("Error: Directory not found: $source_dir")
        return
    end

    # Define the path for non-Hoosier tires
    other_dir = joinpath(source_dir, "Other")

    println("Scanning for Hoosier tires in: $source_dir")
    
    files = readdir(source_dir)
    # Filter for data files, but EXCLUDE the "Other" directory itself if it exists
    data_files = filter(f -> (endswith(f, ".dat") || endswith(f, ".mat")) && f != "Other", files)

    processed_count = 0
    other_count = 0

    for file in data_files
        full_path = joinpath(source_dir, file)
        
        # Skip if it's somehow a directory (safety check)
        if isdir(full_path) continue end

        try
            _, meta_dict = parse_ttc(full_path)

            if haskey(meta_dict, "Tire_Name")
                tire_name = meta_dict["Tire_Name"]
                
                if contains(lowercase(tire_name), "hoosier")
                    # Handle Hoosier tires
                    safe_name = replace(tire_name, r"[\\/:*?\"<>|]" => "_")
                    target_dir = joinpath(source_dir, safe_name)
                    mkpath(target_dir)

                    mv(full_path, joinpath(target_dir, file), force=true)
                    println("✓ Hoosier: $file -> $safe_name/")
                    processed_count += 1
                else
                    # Handle "Other" tires
                    mkpath(other_dir)
                    mv(full_path, joinpath(other_dir, file), force=true)
                    println("○ Other:   $file -> Other/")
                    other_count += 1
                end
            else
                # If no Tire_Name metadata exists, move to Other as well
                mkpath(other_dir)
                mv(full_path, joinpath(other_dir, file), force=true)
                println("! No Metadata: Moving $file to Other/")
                other_count += 1
            end
        catch e
            println("!! Failed to process $file: $e")
        end
    end

    println("---")
    println("Summary: Found $processed_count Hoosier files. Moved $other_count files to 'Other'.")
end

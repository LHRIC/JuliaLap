import Pkg;
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("GLMakie")
Pkg.add("Interpolations")
using CSV
using DataFrames
using GLMakie
using Interpolations

function main()

    #Get file path
    print("Enter CSV file path: ")
    file_path = strip(readline())

    #Check if file for low or high speed
    file_type = (occursin("_LS", file_path) || occursin(" LS", file_path)) ? "LS" : "HS"

    #Load file into dataframe
    df = nothing #Needed for scope
    try
        df = CSV.read(file_path, DataFrame)
    catch e
        println("Error loading file; ", e)
        return
    end

    #Extract settings from column headers
    setting_values = Set{Float64}()
    pattern = file_type == "LS" ? r"\((\d+)-4\.3\) V-" : r"\(\d+-([\d.]+)\) V-" #Used chatgpt, apparently works the same as Python
    for col in names(df)
        m = match(pattern, col)
        if !isnothing(m) #isnothing is preferred over !==
            setting = parse(Float64, m.captures[1]) #First group is not the entire string unlike python
            push!(setting_values, setting)
        end
    end

    if isempty(setting_values)
        println("No valid settings found in CSV file")
        return
    end

    available_settings = sort(collect(setting_values)) #collect(set) converts set to array, sorted in ascending order
    println("Available settings: ", available_settings)

    #Get user input
    print("Enter desired setting: ")
    desired_setting = parse(Float64, readline())

    fig = Figure(size = (1000, 700)) #Inches to pixels conversion using 100 dpi
    ax = Axis(fig[1, 1], xlabel="Velocity (mm/sec)", ylabel="Force (N)")

    if desired_setting in available_settings
        if file_type == "LS"
            setting_str = (isapprox(desired_setting, round(desired_setting)) ? 
           string(Int(round(desired_setting))) * "-4.3" : string(Int(round(desired_setting))) * "-4.3" : string(desired_setting) * "-4.3") #isapprox returns true when difference is within x*1e-5
        else

            #Original code rounded to 1 decimal place, can't figure out why so didn't
            setting_str = (isapprox(desired_setting, round(desired_setting)) ? 
            ("0-" * string(Int(round(desired_setting)))) : ("0-" * string(desired_setting))) #'*' works only if both are strings, but is faster than ',' 
        end
    end

    cols = Dict(
        "v_comp" => "($(setting_str)) V-C", #Velocity in compression (mm/sec)
        "f_comp" => "($(setting_str)) C", #Force in compression (N)
        "v_reb"  => "($(setting_str)) V-R", #Velocity in rebound (mm/sec)
        "f_reb"  => "($(setting_str)) R" #Force in rebound (N)
    )

    comp_vel = nothing
    comp_force = nothing
    reb_vel = nothing
    reb_force = nothing
    try
        df_sub = df[!, [cols["v_comp"], cols["f_comp"], cols["v_reb"], cols["f_reb"]]] #Create sub dataframe from need columns
        df_clean = dropmissing(df_sub) #Drop rows with missing values
        comp_vel   = df_clean[!, cols["v_comp"]]
        comp_force = df_clean[!, cols["f_comp"]]
        reb_vel    = df_clean[!, cols["v_reb"]]
        reb_force  = df_clean[!, cols["f_reb"]]

    catch e
        println("Missing data columns for setting $setting_str: ", e)
        return
    end

    # Convert to array and sort data in increasing order, python does automatically
    sorted_index = sortperm(comp_vel) #Short for sorted permuation
    comp_vel = comp_vel[sorted_index]
    comp_force = comp_force[sorted_index]
    sorted_index = sortperm(reb_vel)
    reb_vel = reb_vel[sorted_index]
    reb_force = reb_force[sorted_index]

    #Compression data points plotted in blue, rebound data points in green
    lines!(ax, comp_vel, comp_force, color=:blue, label = "$(setting_str) Compression")
    lines!(ax, reb_vel, reb_force, color=:green, label = "$(setting_str) Rebound")

    last_val = round(comp_vel[end], digits=2)

    print("Enter velocity to query (0.0 - $last_val): ")
    query_velocity = parse(Float64, readline())

    #Interpolation used is linear, for higher degrees replace Linear() with Quadratic(), Cubic(), etc.
    #Syntax is interpolation function = interpolate(tuple of x coordinates, array of y-coordinates, Gridded(degree polynomial))
    comp_interp = extrapolate(interpolate((comp_vel,), comp_force, Gridded(Linear())), Line()) #(comp_vel,) makes comp_vel a tuple, needed to pass vertical line test
    reb_interp  = extrapolate(interpolate((reb_vel,), reb_force, Gridded(Linear())), Line()) #Gridded creates irregular grids based on spacing between x-coordinates

    #Calculate and store interpolated points with user's specified velocity
    comp_query = comp_interp(query_velocity)
    reb_query = reb_interp(query_velocity)

    println("At velocity $query_velocity mm/s:")
    println("  Compression force (N)≈ $comp_query")
    println("  Rebound force (N)    ≈ $reb_query")

    #Scatter! adds points on top of existing axis (plot)
    scatter!(ax, [query_velocity], [comp_query], color=:yellow, markersize=15, label="Queried Compression")
    scatter!(ax, [query_velocity], [reb_query], color=:blue, markersize=15, label="Queried Rebound")
    display(fig)  #shows the figure window


    #Loop for user input
    while true
        println("Enter velocity (0.0 - $last_val) to query, or 'q' to quit: ")
        input = readline();
        if lowercase(input) == "q"
            break
        end
        query_velocity = parse(Float64, input)
        comp_query = comp_interp(query_velocity)
        reb_query = reb_interp(query_velocity)

        println("At velocity $query_velocity mm/s:")
        println("  Compression force (N)≈ $comp_query")
        println("  Rebound force (N)    ≈ $reb_query")

        #Scatter! adds points on top of existing axis (plot)
        scatter!(ax, [query_velocity], [comp_query], color=:yellow, markersize=15, label=nothing)
        scatter!(ax, [query_velocity], [reb_query], color=:blue, markersize=15, label=nothing)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
        main()
end

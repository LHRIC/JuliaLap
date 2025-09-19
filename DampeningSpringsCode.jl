using CSV
using DataFrames
using GLMakie
using Interpolations

function get_settings(file_path::AbstractString) #AbstractString works for both substrings and strings

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
    return available_settings
end

function get_interpolation(file_path::AbstractString, desired_setting::Float64)

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

    available_settings = get_settings(file_path)

    if desired_setting in available_settings
        if file_type == "LS"
            setting_str = (isapprox(desired_setting, round(desired_setting)) ? 
           string(Int(round(desired_setting))) * "-4.3" : string(desired_setting) * "-4.3") #isapprox returns true when difference is within x*1e-5
        else

            #Original code rounded to 1 decimal place, can't figure out why so didn't
            setting_str = (isapprox(desired_setting, round(desired_setting)) ? 
            ("0-" * string(Int(round(desired_setting)))) : ("0-" * string(desired_setting))) #'*' works only if both are strings, but is faster than ',' 
        end
    else
        println("Enter a valid setting: $(available_settings)")
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

    last_val = round(comp_vel[end], digits=2)

    #Interpolation used is linear, for higher degrees replace Linear() with Quadratic(), Cubic(), etc.
    #Syntax is interpolation function = interpolate(tuple of x coordinates, array of y-coordinates, Gridded(degree polynomial))
    comp_interp = extrapolate(interpolate((comp_vel,), comp_force, Gridded(Linear())), Line()) #(comp_vel,) makes comp_vel a tuple, needed to pass vertical line test
    reb_interp  = extrapolate(interpolate((reb_vel,), reb_force, Gridded(Linear())), Line()) #Gridded creates irregular grids based on spacing between x-coordinates

    return comp_interp, reb_interp, last_val
end

function display_curves(comp_func, reb_func, max_vel)
    fig = Figure(size = (1000, 700)) #Inches to pixels conversion using 100 dpi
    ax = Axis(fig[1, 1], xlabel="Velocity (mm/sec)", ylabel="Force (N)")

    #Plot from 0-max_vel (mm/s)
    comp_x = range(0, max_vel, length=100)
    reb_x = range(0, max_vel, length=100)

    #Evaluate with given ranges
    comp_y = comp_func.(comp_x)
    reb_y = reb_func.(reb_x)


    #Compression data points plotted in blue, rebound data points in green
    lines!(ax, comp_x, comp_y, color=:blue, linewidth=2, label="Compression curve")
    lines!(ax, reb_x, reb_y, color=:red, linewidth=2, label="Rebound curve")

    axislegend(ax) #Displays labels
    display(fig) #Shows the figure window
end

function use_funcs(comp_func, reb_func,)
    return query_vel -> (comp_func(query_vel), reb_func(query_vel)) #Uses a closure which is pretty interesting
end

function main() 
    #Returns available settings for specified file in an array
    println(get_settings("DSD_12_LS.csv")) 

    #Returns interpolation functions for compression and rebound, and also max compression velocity in data set, explicitly cast to Float64 to be safe
    comp_func, reb_func, max_vel = get_interpolation("DSD_12_LS.csv", Float64(2)) 

    #Plots and displays damper curves given compression and rebound interpolation functions, from 0-max_vel (mm/sec)
    display_curves(comp_func, reb_func, max_vel)

    #Setup so that only velocity has to be passed into get_force function
    get_force = use_funcs(comp_func, reb_func)

    #Pass velocity in mm/sec, returns compression and rebound force values in N
    comp_val, reb_val = get_force(100.5)
    println("At velocity 100.5 mm/s: Compression ≈ $comp_val, Rebound ≈ $reb_val")

    sleep(10) #Just to show the graph
end

if abspath(PROGRAM_FILE) == @__FILE__
        main()
end
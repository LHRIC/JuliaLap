module DampeningSprings

using CSV
using DataFrames
using Interpolations
using GLMakie

# ──────────────────────────────────────────────────────────────
# Load CSV and extract available damper settings
# ──────────────────────────────────────────────────────────────
function get_settings(file_path::AbstractString)
    file_type = (occursin("_LS", file_path) || occursin(" LS", file_path)) ? "LS" : "HS"

    df = try
        CSV.read(file_path, DataFrame)
    catch e
        println("Error loading file: ", e)
        return nothing
    end

    setting_values = Set{Float64}()
    pattern = file_type == "LS" ? r"\((\d+)-4\.3\) V-" : r"\(\d+-([\d.]+)\) V-"
    for col in names(df)
        m = match(pattern, col)
        if !isnothing(m)
            push!(setting_values, parse(Float64, m.captures[1]))
        end
    end

    isempty(setting_values) && println("No valid settings found in CSV file")
    return sort(collect(setting_values))
end

# ──────────────────────────────────────────────────────────────
# Create interpolations for compression and rebound
# ──────────────────────────────────────────────────────────────
function get_interpolation(file_path::AbstractString, desired_setting::Float64)
    df = try
        CSV.read(file_path, DataFrame)
    catch e
        println("Error loading file: ", e)
        return nothing
    end

    available_settings = get_settings(file_path)
    if !(desired_setting in available_settings)
        println("Enter a valid setting: $(available_settings)")
        return nothing
    end

    file_type = (occursin("_LS", file_path) || occursin(" LS", file_path)) ? "LS" : "HS"

    # Build column names
    if file_type == "LS"
        setting_str = isapprox(desired_setting, round(desired_setting)) ? "$(Int(round(desired_setting)))-4.3" : "$desired_setting-4.3"
    else
        setting_str = isapprox(desired_setting, round(desired_setting)) ? "0-$(Int(round(desired_setting)))" : "0-$desired_setting"
    end

    cols = Dict(
        "v_comp" => "($setting_str) V-C",
        "f_comp" => "($setting_str) C",
        "v_reb"  => "($setting_str) V-R",
        "f_reb"  => "($setting_str) R"
    )

    df_sub = dropmissing(df[:, [cols["v_comp"], cols["f_comp"], cols["v_reb"], cols["f_reb"]]])
    comp_vel, comp_force = df_sub[:, cols["v_comp"]], df_sub[:, cols["f_comp"]]
    reb_vel,  reb_force = df_sub[:, cols["v_reb"]],  df_sub[:, cols["f_reb"]]

    # Sort
    comp_sort = sortperm(comp_vel)
    reb_sort  = sortperm(reb_vel)
    comp_vel, comp_force = comp_vel[comp_sort], comp_force[comp_sort]
    reb_vel, reb_force = reb_vel[reb_sort], reb_force[reb_sort]

    # Linear interpolation
    comp_interp = extrapolate(interpolate((comp_vel,), comp_force, Gridded(Linear())), Line())
    reb_interp  = extrapolate(interpolate((reb_vel,), reb_force, Gridded(Linear())), Line())
    last_val = round(comp_vel[end], digits=2)

    return comp_interp, reb_interp, last_val
end

# ──────────────────────────────────────────────────────────────
# Return the compression and rebound functions for given LS/HS settings
# ──────────────────────────────────────────────────────────────
function get_func(desired_settings::Tuple{Float64, Float64})
    LS_interpolations = [get_interpolation("DSD_12_LS.csv", i) for i in get_settings("DSD_12_LS.csv")]
    HS_interpolations = [get_interpolation("DSD_12_HS.csv", i) for i in get_settings("DSD_12_HS.csv")]

    LS_available = get_settings("DSD_12_LS.csv")
    HS_available = get_settings("DSD_12_HS.csv")
    LS_setting, HS_setting = desired_settings

    if !(LS_setting in LS_available && HS_setting in HS_available)
        println("Invalid LS/HS setting")
        return nothing
    end

    # Mix LS and HS functions as in original logic
    if HS_setting == 4.3
        return LS_interpolations[findfirst(==(LS_setting), LS_available)]
    elseif LS_setting == 0
        return HS_interpolations[findfirst(==(HS_setting), HS_available)]
    else
        diff_comp_func = x -> HS_interpolations[findfirst(==(HS_setting), HS_available)][1](x) +
                             (LS_interpolations[findfirst(==(LS_setting), LS_available)][1](x) -
                              LS_interpolations[1][1](x))
        diff_reb_func = x -> HS_interpolations[findfirst(==(HS_setting), HS_available)][2](x) +
                            (LS_interpolations[findfirst(==(LS_setting), LS_available)][2](x) -
                             LS_interpolations[1][2](x))
        return diff_comp_func, diff_reb_func, HS_interpolations[findfirst(==(HS_setting), HS_available)][3]
    end
end

# ──────────────────────────────────────────────────────────────
# Public exports
# ──────────────────────────────────────────────────────────────
export get_settings, get_interpolation, get_func

end # module

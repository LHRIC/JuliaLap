using CSV
using DataFrames
using Interpolations
using Plots

function get_dataframe(file_path::AbstractString)
    df = CSV.read(file_path, DataFrame)
    return df
    
end

df = get_dataframe("38normalized.csv") #"38normalized.csv", or "26normalized.csv"

x = df[:,1] #RPM
y = df[:,2] #Torque

lin_int = LinearInterpolation(x, y)

rpms = range(first(x), last(x), length=200)
torques = lin_int.(rpms)

plot(
    rpms,
    torques,
    label = "Linearly Interpolated Torque Curve",
    xlabel = "RPM",
    ylabel = "Torque (Nm)",
    linewidth = 2,
    color = :blue,
)





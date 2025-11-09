module TorqueCurves

using CSV
using DataFrames
using Interpolations
using Plots

# ──────────────────────────────────────────────────────────────
# Load and interpolate torque curve data
# ──────────────────────────────────────────────────────────────

"""
    get_torque_interp()

Reads `"38normalized.csv"` and returns a LinearInterpolation
mapping RPM → Torque (Nm).
"""
function get_torque_interp(file_path::AbstractString)
    df = CSV.read(file_path, DataFrame)
    x = df[:, 1]  # RPM
    y = df[:, 2]  # Torque
    lin_int = LinearInterpolation(x, y)
    return lin_int
end

# ──────────────────────────────────────────────────────────────
# Optional visualization helper
# ──────────────────────────────────────────────────────────────
function plot_torque_curve()
    lin_int = get_torque_interp()
    df = CSV.read("38normalized.csv", DataFrame)
    x = df[:, 1]
    rpms = range(first(x), last(x), length=200)
    torques = lin_int.(rpms)

    plt = plot(
        rpms,
        torques,
        label="Linearly Interpolated Torque Curve",
        xlabel="RPM",
        ylabel="Torque (Nm)",
        linewidth=2,
        color=:blue,
        title="Torque vs RPM"
    )
    display(plt)
    return plt
end

# ──────────────────────────────────────────────────────────────
# Public exports
# ──────────────────────────────────────────────────────────────
export get_torque_interp, plot_torque_curve

end # module

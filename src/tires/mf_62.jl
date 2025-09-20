
module MF62

export MF62Model

"""
Magic Formula 6.2 Tire Model
(Pacejka, Tire and Vehicle Dynamics, 3rd edition, Ch. 4.3.2)
"""

mutable struct MF62Model
    params::Dict{String,Any}
    fz0::Float64
    stiffness_tracker::Vector{Vector{Float64}}
end

function MF62Model(params::Dict{String,Any})
    return MF62Model(params, 800.0, Vector{Vector{Float64}}())
end

"""
Compute longitudinal force (pure slip, α = 0).
"""
function fx(model::MF62Model, fz::Float64, kappa::Float64, gamma::Float64)
    p = model.params
    dfz = (fz - model.fz0) / model.fz0

    c_x = p["PCX1"] * 1
    @assert c_x > 0 "c_x must be > 0"

    mu_x = (p["PDX1"] + p["PDX2"] * dfz) *
           (1 + p["PPX3"] * 1 + p["PPX4"] * 1^2) *
           (1 - p["PDX3"] * gamma^2) * 1

    # Saturation to prevent negative mu
    mu_x = mu_x < 1e-2 ? 1e-2 : mu_x

    s_hx = (p["PHX1"] + p["PHX2"] * dfz) * 1
    kappa_x = kappa + s_hx

    e_x = (p["PEX1"] + p["PEX2"] * dfz + p["PEX3"] * dfz^2) *
          (1 - p["PEX4"] * sign(kappa_x)) * 1
    @assert e_x <= 1 "e_x must be <= 1"

    d_x = mu_x * fz * 1
    @assert d_x > 0 "d_x must be > 0"

    k_xk = fz * (p["PKX1"] + p["PKX2"] * dfz) *
           exp(p["PKX3"] * dfz) * (1 + p["PPX1"] * 1 + p["PPX2"] * 1^2)

    # Enforce stiffness saturation
    sp_sat = 12.0
    if (p["PKX1"] + p["PKX2"] * dfz) < sp_sat
        k_xk = fz * (sp_sat * exp(1e-2 * ((p["PKX1"] + p["PKX2"] * dfz) - sp_sat))) *
               exp(p["PKX3"] * dfz) * (1 + p["PPX1"] * 1 + p["PPX2"] * 1^2)
    end

    b_x = k_xk / (c_x * d_x)
    s_vx = fz * (p["PVX1"] + p["PVX2"] * dfz) * 1 * 1 * 1

    fx0 = d_x * sin(c_x * atan(b_x * kappa_x - e_x *
           (b_x * kappa_x - atan(b_x * kappa_x)))) + s_vx

    push!(model.stiffness_tracker, [k_xk, b_x])

    return fx0 * p["friction_scaling_x"]
end

"""
Compute lateral force (pure slip, κ = 0).
"""
function fy(model::MF62Model, fz::Float64, alpha::Float64, gamma::Float64)
    p = model.params
    dfz = (fz - model.fz0) / model.fz0

    c_y = p["PCY1"] * 1
    @assert c_y > 0 "c_y must be > 0"

    mu_y = (p["PDY1"] + p["PDY2"] * dfz) *
           (1 + p["PPY3"] * 1 + p["PPY4"] * 1^2) *
           (1 - p["PDY3"] * gamma^2) * 1

    d_y = mu_y * fz * 1

    k_ya = p["PKY1"] * model.fz0 * (1 + p["PPY1"] * 1) *
           (1 - p["PKY3"] * abs(gamma)) *
           sin(p["PKY4"] * atan(
               (fz / model.fz0) /
               ((p["PKY2"] + p["PKY5"] * gamma^2) * (1 + p["PPY2"] * 1))
           )) * 1 * 1

    b_y = k_ya / (c_y * d_y)

    s_vyy = fz * (p["PVY3"] + p["PVY4"] * dfz) * gamma * 1 * 1 * 1
    s_vy = fz * (p["PVY1"] + p["PVY2"] * dfz) * 1 * 1 * 1 + s_vyy

    k_yy0 = fz * (p["PKY6"] + p["PKY7"] * dfz) * (1 + p["PPY5"] * 1) * 1

    s_hy = (p["PHY1"] + p["PHY2"] * dfz) * 1 + ((k_yy0 * gamma - s_vyy) / (k_ya)) * 1
    alpha_y = alpha + s_hy

    e_y = (p["PEY1"] + p["PEY2"] * dfz) *
          (1 + p["PEY5"] * gamma^2 -
           (p["PEY3"] + p["PEY4"] * gamma) * sign(alpha_y)) * 1
    @assert e_y <= 1 "e_y must be <= 1"

    fy0 = d_y * sin(c_y * atan(b_y * alpha_y - e_y *
          (b_y * alpha_y - atan(b_y * alpha_y)))) + s_vy

    return fy0 * p["friction_scaling_y"]
end

end # module
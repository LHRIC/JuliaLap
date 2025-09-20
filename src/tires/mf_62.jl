
module MF62

export MF62Model

# Magic Formula 6.2 Tire Model (Pacejka, Tire and Vehicle Dynamics, 3rd edition, Ch. 4.3.2)

mutable struct MF62Model
    params::Dict{String,Any}
    fz0::Float64 
    stiffness_tracker::Vector{Vector{Float64}}
end

function MF62Model(params::Dict{String,Any})
    return MF62Model(params, params["FNOMIN"], Vector{Vector{Float64}}())
end

# Compute longitudinal force (pure slip, α = 0)
# Neglecting turn slip, and assuming small camber values (Lamba=1)
# 1"s represent un-used user correction coefficents and un-implemented tire pressure sensitivity

function fx(model::MF62Model, fz::Float64, kappa::Float64, gamma::Float64)

       # initialize variables from parse tire dictionary 
       #old 
       p = model.params
       pcx1 = p["PCX1"]
       pdx1 = p["PDX1"]
       pdx2 = p["PDX2"]
       pdx3 = p["PDX3"]
       pex1 = p["PEX1"]
       pex2 = p["PEX2"]
       pex3 = p["PEX3"]
       pex4 = p["PEX4"]
       pkx1 = p["PKX1"]
       pkx2 = p["PKX2"]
       pkx3 = p["PKX3"]
       phx1 = p["PHX1"]
       phx2 = p["PHX2"]
       pvx1 = p["PVX1"]
       pvx2 = p["PVX2"]
       ppx1 = p["PPX1"]
       ppx2 = p["PPX2"]
       ppx3 = p["PPX3"]
       ppx4 = p["PPX4"]

       # new additions 
       lfz0 = p["LFZO"]                   # Scale factor of nominal (rated) load
       p_i = ["INFLPRES"]                 # Tire inflation pressure
       p_io = ["NOMPRES"]                 # Nominal inflation pressure
       lmux = ["LMUX"]                    # Peak friction coefficient
       r0 = ["UNLOADED_RADIUS"]           # Unloaded tire radius 
       g = ["GRAVITY"]                    # Gravity 
       v0 = sqrt((g*r0))                  # Derived reference velocity 
       lmux = ["LMUX"]                    # Peak friction coefficient
       lmuv = 1                           # Not in tire model: slip speed Vs decaying friction 
       lcx = ["LCX"]                      # Shape factor
       lhx = ["LHX"]                      # Horizontal shift
       lex = ["LEX"]                      # Curvature factor
       lvx = ["LVX"]                      # Vertical shift
       friction_scaling_x = ["friction_scaling_x"]

       # make it instead kappa being input, Vsx and Vcx as input 
       # either input kappa or the two velocities 
    
       fz0p = lfz0 * model.fz0            # (4.E1)

       dfz = (fz - fz0p) / fz0p           # (4.E2a)
       dpi = (p_i - p_io) / p_io          # (4.E2b)

       # come back for 4.E3
       gam_str = sin(gamma)               # (4.E4)
       # come back for 4.E5, 4.E6, 4.E6a
       lmux_str = 1                       # (4.E7)
       A_mu = 10                          # (4.E8)
       lmux_p = (A_mu*lmux_str)/(1+((A_mu - 1)*lmux_str))

       # longitudinal force (alpha = 0)
       c_x = pcx1*lcx                      # (4.E11)
       mux = (pdx1 + (pdx2*dfz)) * (1+(ppx3*dpi) + (ppx4*(dpi^2))) * (1 - (pdx3*(gamma^2))) * lmux_str   # (4.E13)
       genie = 1                          # (i = 0, 1, ..., 8)
       d_x = mux * fz * genie              # (4.E12)
       s_hx = (phx1 + (phx2*dfz))*lhx      # (4.E17)
       kappa_x = kappa + s_hx              # (4.E10)
       e_x = (pex1 + (pex2*dfz) + pex3*(dfz^2)) * (1 - (pex4*sign(kappa_x))) * lex          # (4.E14)
       k_xk = fz * (pkx1 + (pkx2*dfz)) * (exp(pkx3*dfz)) * (1 + (ppx1*dpi) + (ppx2*(dpi^2)))       # (4.E15) note: unknown thing under equation questionable
       epsilon = 0                         # error amount, assume this is perfecto
       b_x = k_xk/((c_x * d_x) + epsilon)  # (4.E16)
       s_vx = fz * (pvx1 + (pvx2*dfz)) * lvx * lmux_p * genie         # (4.E18)
       fx0 = d_x * (sin(c_x * atan(b_x*kappa_x - e_x*(b_x*kappa_x - atan(b_x*kappa_x))))) + s_vx         # (4.E9)

       return fx0*friction_scaling_x
end

# Compute lateral force (pure slip, κ = 0)

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
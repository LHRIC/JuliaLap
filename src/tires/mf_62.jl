
module MF62
# Magic Formula 6.2 Tire Model (Pacejka, Tire and Vehicle Dynamics, 3rd edition, Ch. 4.3.2)

# Compute longitudinal force (pure slip, α = 0)
# Neglecting turn slip, and assuming small camber values (Lamba=1)
# 1"s represent un-used user correction coefficents and un-implemented tire pressure sensitivity
function fx(params::Dict{String,Any}, fz, kappa, gamma)

       # initialize variables from parse tire dictionary 
       #old 
       p = params
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
       p_i = p["INFLPRES"]                 # Tire inflation pressure
       p_io = p["NOMPRES"]                 # Nominal inflation pressure
       lmux = p["LMUX"]                    # Peak friction coefficient
       r0 = p["UNLOADED_RADIUS"]           # Unloaded tire radius 
       g = p["GRAVITY"]                    # Gravity 
       v0 = sqrt(abs(g*r0))                # Derived reference velocity
       lmuv = 1                            # Not in tire model: slip speed Vs decaying friction 
       lcx = p["LCX"]                      # Shape factor
       lhx = p["LHX"]                      # Horizontal shift
       lex = p["LEX"]                      # Curvature factor
       lvx = p["LVX"]                      # Vertical shift
       friction_scaling_x = p["friction_scaling_x"]

       fz0 = p["FNOMIN"]

       # TODO: make it instead kappa being input, Vsx and Vcx as input 
       # either input kappa or the two velocities 
    
       fz0p = lfz0 * fz0            # (4.E1)

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
       zeta = 1                            # (i = 0, 1, ..., 8) 
       d_x = mux * fz * zeta               # (4.E12)
       s_hx = (phx1 + (phx2*dfz))*lhx      # (4.E17)
       kappa_x = kappa + s_hx              # (4.E10)
       e_x = (pex1 + (pex2*dfz) + pex3*(dfz^2)) * (1 - (pex4*sign(kappa_x))) * lex          # (4.E14)
       k_xk = fz * (pkx1 + (pkx2*dfz)) * (exp(pkx3*dfz)) * (1 + (ppx1*dpi) + (ppx2*(dpi^2)))       # (4.E15) note: unknown thing under equation questionable
       epsilon = 0                         # error amount, assume this is perfecto
       b_x = k_xk/((c_x * d_x) + epsilon)  # (4.E16)
       s_vx = fz * (pvx1 + (pvx2*dfz)) * lvx * lmux_p * zeta        # (4.E18)
       fx0 = d_x * (sin(c_x * atan(b_x*kappa_x - e_x*(b_x*kappa_x - atan(b_x*kappa_x))))) + s_vx         # (4.E9)

       return fx0*friction_scaling_x
end

# Compute lateral force (pure slip, κ = 0)
function fy(params::Dict{String,Any}, fz, alpha, gamma)

       # old 
       p = params
       pcy1 = p["PCY1"]
       pdy1 = p["PDY1"]
       pdy2 = p["PDY2"]
       pdy3 = p["PDY3"]
       pey1 = p["PEY1"]
       pey2 = p["PEY2"]
       pey3 = p["PEY3"]
       pey4 = p["PEY4"]
       pey5 = p["PEY5"]
       pky1 = p["PKY1"]
       pky2 = p["PKY2"]
       pky3 = p["PKY3"]
       pky4 = p["PKY4"]
       pky5 = p["PKY5"]
       pky6 = p["PKY6"]
       pky7 = p["PKY7"]
       phy1 = p["PHY1"]
       phy2 = p["PHY2"]
       pvy1 = p["PVY1"]
       pvy2 = p["PVY2"]
       pvy3 = p["PVY3"]
       pvy4 = p["PVY4"]
       ppy1 = p["PPY1"]
       ppy2 = p["PPY2"]
       ppy3 = p["PPY3"]
       ppy4 = p["PPY4"]
       ppy5 = p["PPY5"]

       # new additions 
       lfz0 = p["LFZO"]                   # Scale factor of nominal (rated) load
       p_i = p["INFLPRES"]                 # Tire inflation pressure
       p_io = p["NOMPRES"]                 # Nominal inflation pressure
       lmuy = p["LMUY"]                    # Peak friction coefficient
       r0 = p["UNLOADED_RADIUS"]           # Unloaded tire radius 
       g = p["GRAVITY"]                    # Gravity 
       v0 = sqrt(abs(g*r0))                # Derived reference velocity 
       lmuv = 1                            # Not in tire model: slip speed Vs decaying friction 
       lcy = p["LCY"]                      # Shape factor
       lhy = p["LHY"]                      # Horizontal shift
       ley = p["LEY"]                      # Curvature factor
       lvy = p["LVY"]                      # Vertical shift
       lky = p["LKY"]                      # Cornering stiffness
       lkyg = p["LKYC"]
       friction_scaling_y = p["friction_scaling_y"]

       fz0 = p["FNOMIN"]

       fz0p = lfz0 * fz0            # (4.E1)

       dfz = (fz - fz0p) / fz0p           # (4.E2a)
       dpi = (p_i - p_io) / p_io          # (4.E2b)

       # come back for 4.E3
       gam_str = sin(gamma)               # (4.E4)
       # come back for 4.E5, 4.E6, 4.E6a
       lmuy_str = 1                       # (4.E7)
       A_mu = 10                          # (4.E8)
       lmuy_p = (A_mu*lmuy_str)/(1+((A_mu - 1)*lmuy_str))

       # alpha_str needs velocity, have to implement for the other method, assume it is 1 for now 
       alpha_str = tan(alpha)
       zeta = 1                       # (i = 0, 1, ..., 8)
       k_ya = pky1*fz0p*(1+(ppy1*dpi))*(1-(pky3*abs(gam_str)))*sin(pky4*atan((fz/fz0p)/((pky2+(pky5*(gam_str^2)))*(1+(ppy2*dpi)))))*zeta*lky            # (4.E25)
       k_yg0 = fz*(pky6 + (pky7*dfz))*(1 + (ppy5*dpi))*lkyg          # (4.E30)
       s_vyg = fz*(pvy3+(pvy4*dfz))*gam_str*lkyg*lmuy_p*zeta  # (4.E28)
       epsilon_K = 0                      # Assume for now there is no error in kappa 
       epsilon_y = 0                      # Assume for now there is no error in y 
       s_hy = (phy1 + (phy2*dfz))*lhy + (((k_yg0*gam_str) - s_vyg)/(k_ya + epsilon_K))*zeta + zeta - 1 # (4.E27)
       alpha_y = alpha_str + s_hy                              # (4.E20)
       muy = (pdy1 + (pdy2*dfz))*(1 + (ppy3*dpi) + (ppy4*(dpi^2)))*(1 - (pdy3*(gam_str^2)))*lmuy_str            # (4.E23)
       d_y = muy*fz*zeta                                      # (4.E22)
       c_y = pcy1*lcy                                          # (4.E21)
       @assert c_y > 0 "c_y is less than 0"
       b_y = k_ya/(c_y*d_y + epsilon_y)                        # (4.E26)
       e_y = (pey1 + pey2*dfz)*(1+pey5*(gam_str^2) - (pey3 + pey4*gam_str)*sign(alpha_y))*ley     # (4.E24)
       @assert e_y <= 1 "e_y is greater than 1"
       s_vy = fz*(pvy1 + pvy2*dfz)*lvy*lmuy_p*zeta + s_vyg          # (4.E29)
       fy0 = d_y*sin(c_y*atan(b_y*alpha_y - e_y*(b_y*alpha_y - atan(b_y*alpha_y)))) + s_vy

       return fy0*friction_scaling_y
end

# Aligning Torque (pure slip slip, kappa = 0)
function at(params::Dict{String,Any}, fz, alpha, gamma)
    p = params
    qhz1 = p["QHZ1"]
    qhz2 = p["QHZ2"]
    qhz3 = p["QHZ3"]
    qhz4 = p["QHZ4"]
    qbz1 = p["QBZ1"]
    qbz2 = p["QBZ2"]
    qbz3 = p["QBZ3"]
    qbz5 = p["QBZ5"]
    qbz6 = 0 # p["QBZ6"]
    qbz9 = p["QBZ9"]
    qbz10 = p["QBZ10"]
    qcz1 = p["QCZ1"]
    qdz1 = p["QDZ1"]
    qdz2 = p["QDZ2"]
    qdz3 = p["QDZ3"]
    qdz4 = p["QDZ4"]
    qdz6 = p["QDZ6"]    # 67
    qdz7 = p["QDZ7"]    # 67
    qdz8 = p["QDZ8"]
    qdz9 = p["QDZ9"]
    qdz10 = p["QDZ10"]
    qdz11 = p["QDZ11"]
    ppz1 = p["PPZ1"]
    ppz2 = p["PPZ2"]
    qez1 = p["QEZ1"]
    qez2 = p["QEZ2"]
    qez3 = p["QEZ3"]
    qez4 = p["QEZ4"]
    qez5 = p["QEZ5"]

    # needed fy variables 
    pcy1 = p["PCY1"]
    pdy1 = p["PDY1"]
    pdy2 = p["PDY2"]
    pdy3 = p["PDY3"]
    pey1 = p["PEY1"]
    pey2 = p["PEY2"]
    pey3 = p["PEY3"]
    pey4 = p["PEY4"]
    pey5 = p["PEY5"]
    pky1 = p["PKY1"]
    pky2 = p["PKY2"]
    pky3 = p["PKY3"]
    pky4 = p["PKY4"]
    pky5 = p["PKY5"]
    pky6 = p["PKY6"]
    pky7 = p["PKY7"]
    phy1 = p["PHY1"]
    phy2 = p["PHY2"]
    pvy1 = p["PVY1"]
    pvy2 = p["PVY2"]
    pvy3 = p["PVY3"]
    pvy4 = p["PVY4"]
    ppy1 = p["PPY1"]
    ppy2 = p["PPY2"]
    ppy3 = p["PPY3"]
    ppy4 = p["PPY4"]
    ppy5 = p["PPY5"]

    lky = p["LKY"]
    lmuy = p["LMUY"]
    r0 = p["UNLOADED_RADIUS"]           # Unloaded tire radius 
    ltr = p["LTR"]                      # Pneumatic trail
    lres = p["LRES"]                    # Residual torque
    lkzc = p["LKZC"]                    # Camber torque stiffness
    lhy = p["LHY"]                      # Horizontal shift
    lfz0 = p["LFZO"]                    # Scale factor of nominal (rated) load
    p_i = p["INFLPRES"]                 # Tire inflation pressure
    p_io = p["NOMPRES"]                 # Nominal inflation pressure
    lcy = p["LCY"]                      # Shape factor
    lkyc = p["LKYC"]                    # Camber force stiffness
    lvy = p["LVY"]                      # Vertical shift


    fz0 = p["FNOMIN"]
    
    fz0p = lfz0 * fz0            # (4.E1)

    dfz = (fz - fz0p) / fz0p           # (4.E2a)
    dpi = (p_i - p_io) / p_io          # (4.E2b)

    # come back for 4.E3
    gam_str = sin(gamma)               # (4.E4)
    # come back for 4.E5, 4.E6, 4.E6a
    lmux_str = 1                       # (4.E7)
    A_mu = 10                          # (4.E8)
    lmux_p = (A_mu*lmux_str)/(1+((A_mu - 1)*lmux_str))
    alpha_str = tan(alpha)
    epsilon_K = 0 
    epsilon_y = 0
    lmuy_str = lmuy                 # (4.E7) TODO: velocity 
    zeta = 1

    lmuy_p = (A_mu*lmuy_str)/(1+((A_mu - 1)*lmuy_str))

    d_r = fz*r0*((qdz6+qdz7)*lres*zeta + ((qdz8 + qdz9*dfz)*(1 + ppz2*dpi) + (qdz10 + qdz11*dfz)*abs(gam_str))*gam_str*lkzc*zeta)*lmuy_str*cos(alpha) + zeta - 1 # (4.E47)
    c_r = zeta                      # (4.E46)
    muy = (pdy1 + (pdy2*dfz))*(1 + (ppy3*dpi) + (ppy4*(dpi^2)))*(1 - (pdy3*(gam_str^2)))*lmuy_str            # (4.E23)
    d_y = muy*fz*zeta                                      # (4.E22)
    c_y = pcy1*lcy                                          # (4.E21)
    @assert c_y > 0 "c_y is less than 0"
    k_ya = pky1*fz0p*(1+(ppy1*dpi))*(1-(pky3*abs(gam_str)))*sin(pky4*atan((fz/fz0p)/((pky2+(pky5*(gam_str^2)))*(1+(ppy2*dpi)))))*zeta*lky            # (4.E25)
    b_y = k_ya/(c_y*d_y + epsilon_y)                        # (4.E26)
    b_r = ((qbz9*lky)/(lmuy_str + qbz10*b_y*c_y))*zeta          # (4.E45)
    d_t0 = fz*(r0/fz0p)*(qdz1 + qdz2*dfz)*(1-ppz1*dpi)*ltr      # (4.E42) TODO: velocity sign 
    d_t = d_t0*(1+qdz3*abs(gam_str) + qdz4*gam_str^2)*zeta      # (4.E43)
    c_t = qcz1                    # (4.E41)   
    @assert c_t > 0
    b_t = (qbz1 + qbz2*dfz + qbz3*dfz^2)*(1 + qbz5*abs(gam_str) + qbz6*gam_str^2)*(lky/lmuy_str) # (4.E40)
    k_ya_p = k_ya + epsilon_K                              # (4.E39)
    k_yg0 = fz*(pky6 + (pky7*dfz))*(1 + (ppy5*dpi))*lkyc          # (4.E30)
    s_vyg = fz*(pvy3+(pvy4*dfz))*gam_str*lkyc*lmuy_p*zeta  # (4.E28)
    s_hy = (phy1 + (phy2*dfz))*lhy + (((k_yg0*gam_str) - s_vyg)/(k_ya + epsilon_K))*zeta + zeta - 1 # (4.E27)
    s_vy = fz*(pvy1 + pvy2*dfz)*lvy*lmuy_p*zeta + s_vyg          # (4.E29)
    s_hf = s_hy + s_vy/k_ya_p                              # (4.E38)
    a_r = alpha_str + s_hf                                 # (4.E37)
    m_zr0 = d_r*cos(c_r*atan(b_r*a_r))*cos(alpha)          # (4.E36)
    s_ht = qhz1 + qhz2*dfz + (qhz3 + qhz4*dfz)gam_str      # (4.E35)
    a_t = alpha_str + s_ht      # (4.E34)
    e_t = (qez1 + qez2*dfz + qez3*dfz^2)*(1+(qez4 + qez5*gam_str)*(2/pi)*atan(b_t*c_t*a_t))     # (4.E44)
    @assert e_t <= 1
    t_0 = d_t*cos(c_t*atan(b_t*a_t - e_t*(b_t*a_t-atan(b_t*a_t))))*cos(alpha)       # (4.E33) TODO: cos prime -> velocity implementation 
    f_y0 = fy(p, fz, alpha, gamma)  # function call to fy TODO: want to multiply by friction scalling 
    m_z0_p = -t_0*f_y0          # (4.E32)
    m_z0 = m_z0_p + m_zr0       # (4.E31)
    return m_z0
end

end # module
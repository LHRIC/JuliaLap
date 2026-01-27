module MF62
using ModelingToolkit
# using Symbolics
# Magic Formula 6.2 Tire Model (Pacejka, Tire and Vehicle Dynamics, 3rd edition, Ch. 4.3.2)

params = @variables begin
    lfz0                              
    p_io                              
    lmux                                
    lmuy                              
    lmuv                           
    r0                        
    g                                
    fz0

    pcx1
    pdx1
    pdx2
    pdx3
    pex1
    pex2
    pex3
    pex4
    pkx1
    pkx2
    pkx3
    phx1
    phx2
    pvx1
    pvx2
    ppx1
    ppx2
    ppx3
    ppx4

    lcx
    lhx
    lex
    lvx

    friction_scaling_x

    pcy1
    pdy1
    pdy2
    pdy3
    pey1
    pey2
    pey3
    pey4
    pey5
    pky1
    pky2
    pky3
    pky4
    pky5
    pky6
    pky7
    phy1
    phy2
    phy3
    pvy1
    pvy2
    pvy3
    pvy4
    ppy1
    ppy2
    ppy3
    ppy4
    ppy5

    lcy
    lhy
    ley
    lvy
    lky
    lkyg

    friction_scaling_y

    r0
    zeta
    dfz
    dpi
    gam_str
    lmuy_str
    b_y
    c_y
    alpha_str
    fz0p
    k_ya
    s_hy
    s_vy

    qhz1
    qhz2
    qhz3
    qhz4
    qbz1
    qbz2
    qbz3
    qbz5
    qbz6
    qbz9
    qbz10
    qcz1
    qdz1
    qdz2
    qdz3
    qdz4
    qdz6
    qdz7
    qdz8
    qdz9
    qdz10
    qdz11
    ppz1
    ppz2
    qez1
    qez2
    qez3
    qez4
    qez5

    lky
    ltr
    lres
    lkzc

    rbx1
    rbx2
    rbx3
    rcx1
    rex1
    rex2
    rhx1

    lxal

    rby1
    rby2
    rby3
    rby4
    rcy1
    rey1
    rey2
    rhy1
    rhy2
    rvy1
    rvy2
    rvy3
    rvy4
    rvy5
    rvy6

    lyk
    lvyk
    lvmx
    lmx
    
    qsx1
    qsx2
    qsx3
    qsx4
    qsx5
    qsx6
    qsx7
    qsx8
    qsx9
    qsx10
    qsx11
    ppmx1
 
    qsy1 
    qsy2
    qsy3
    qsy4
    qsy5
    qsy6
    qsy7
    qsy8

    lmy

    ssz1
    ssz2
    ssz3
    ssz4

    ls
end

indep_vars = @parameters begin
    p_i
    fz
    gamma
    alpha
    vs
    Ω
    rl
end

# Pre-processing velocity calculations

vc = Ω*rl
vc_x = vc*cos(alpha)
vc_y = vc*sin(alpha)

vs_x = vs*cos(alpha)
vs_y = vs*sin(alpha)

# Compute longitudinal force (pure slip, α = 0)
# Neglecting turn slip, and assuming small camber values (Lamba=1)
# 1"s represent un-used user correction coefficents and un-implemented tire pressure sensitivity

### fx0

v0 = sqrt(abs(g*r0))                                # Derived reference velocity

fz0p = lfz0 * fz0                                   # (4.E1)
dfz = (fz - fz0p) / fz0p                            # (4.E2a)
dpi = (p_i - p_io) / p_io                           # (4.E2b)
alpha_str = tan(alpha) * sign(vc_x)                 # (4.E3)
gam_str = sin(gamma)                                # (4.E4)
kappa = - vs_x/abs(vc_x)                            # (4.E5)
epsilon = 0                                         # TODO
vcp = vc - epsilon                                  # (4.E6a)
cosalpha_p = vc_x/vcp                               # (4.E6) TODO: above
# cosalpha_p = cos(alpha)
zeta = 1                                            # TODO: check i part

lmux_str = lmux/(1+lmuv*vs/v0)                      # (4.E7)
lmuy_str = lmuy/(1+lmuv*vs/v0)                      # (4.E7)
A_mu = 10                                       
lmux_p = (A_mu*lmux_str)/(1+((A_mu - 1)*lmux_str))  # (4.E8)
lmuy_p = (A_mu*lmuy_str)/(1+((A_mu - 1)*lmuy_str))  # (4.E8)

# longitudinal force (alpha = 0)
s_vx = fz * (pvx1 + (pvx2*dfz)) * lvx * lmux_p * zeta        # (4.E18)
s_hx = (phx1 + (phx2*dfz))*lhx                      # (4.E17)
c_x = pcx1*lcx                                      # (4.E11)
# @assert c_x > 0
mux = (pdx1 + (pdx2*dfz)) * (1+(ppx3*dpi) + (ppx4*(dpi^2))) * (1 - (pdx3*(gamma^2))) * lmux_str   # (4.E13)
d_x = mux * fz * zeta                               # (4.E12)
# @assert d_x > 0
kappa_x = kappa + s_hx                              # (4.E10)
e_x = (pex1 + (pex2*dfz) + pex3*(dfz^2)) * (1 - (pex4*sign(kappa_x))) * lex                     # (4.E14)
# @assert e_x <= 1

k_xk = fz * (pkx1 + (pkx2*dfz)) * (exp(pkx3*dfz)) * (1 + (ppx1*dpi) + (ppx2*(dpi^2)))       # (4.E15) note: unknown thing under equation questionable
epsilon = 0                         # error amount, assume this is perfecto
b_x = k_xk/((c_x * d_x) + epsilon)                  # (4.E16)
fx0 = d_x * (sin(c_x * atan(b_x*kappa_x - e_x*(b_x*kappa_x - atan(b_x*kappa_x))))) + s_vx         # (4.E9)
f_x0 = fx0*friction_scaling_x

### fy0

epsilon_K = 0                      # Assume for now there is no error in kappa 
epsilon_y = 0                      # Assume for now there is no error in y 

k_ya0 = pky1*fz0p*sin(pky4*atan(fz/(pky2*fz0p)))*lkyg
k_yg0 = (phy3*k_ya0 + fz*(pvy3+pvy4*dfz))*lkyg
# k_yg0 = fz*(pky6 + (pky7*dfz))*(1 + (ppy5*dpi))*lkyg          # (4.E30)
s_vyg = fz*(pvy3+(pvy4*dfz))*gam_str*lkyg*lmuy_p*zeta  # (4.E28)
s_vy = fz*(pvy1 + pvy2*dfz)*lvy*lmuy_p*zeta + s_vyg          # (4.E29)
k_ya = pky1*fz0p*(1+(ppy1*dpi))*(1-(pky3*abs(gam_str)))*sin(pky4*atan((fz/fz0p)/((pky2+(pky5*(gam_str^2)))*(1+(ppy2*dpi)))))*zeta*lky            # (4.E25)
s_hy = (phy1 + (phy2*dfz))*lhy + phy3*gam_str*lkyg
# s_hy = (phy1 + (phy2*dfz))*lhy + (((k_yg0*gam_str) - s_vyg)/(k_ya + epsilon_K))*zeta + zeta - 1 # (4.E27)
c_y = pcy1*lcy                                          # (4.E21)
# @assert c_y > 0 "c_y is less than 0"
mu_y = (pdy1 + (pdy2*dfz))*(1 + (ppy3*dpi) + (ppy4*(dpi^2)))*(1 - (pdy3*(gam_str^2)))*lmuy_str            # (4.E23)
d_y = mu_y*fz*zeta                                      # (4.E22)
b_y = k_ya/(c_y*d_y + epsilon_y)                        # (4.E26)
alpha_y = alpha_str + s_hy                              # (4.E20)
e_y = (pey1 + pey2*dfz)*(1+pey5*(gam_str^2) - (pey3 + pey4*gam_str)*sign(alpha_y))*ley     # (4.E24)
# @assert e_y <= 1 "e_y is greater than 1"
fy0 = d_y*sin(c_y*atan(b_y*alpha_y - e_y*(b_y*alpha_y - atan(b_y*alpha_y)))) + s_vy

f_y0 = fy0*friction_scaling_y

# Aligning Torque (pure slip slip, kappa = 0)
epsilon_K = 0

d_r = fz*r0*((qdz6+qdz7)*lres*zeta + ((qdz8 + qdz9*dfz)*(1 + ppz2*dpi) + (qdz10 + qdz11*dfz)*abs(gam_str))*gam_str*lkzc*zeta)*lmuy_str*cos(alpha) + zeta - 1 # (4.E47)
c_r = zeta                      # (4.E46)
b_r = (qbz9*lky/lmuy_str + qbz10*b_y*c_y)*zeta          # (4.E45)
b_t = (qbz1 + qbz2*dfz + qbz3*dfz^2)*(1 + qbz5*abs(gam_str) + qbz6*gam_str^2)*(lky/lmuy_str) # (4.E40)
# @assert b_t > 0
c_t = qcz1                    # (4.E41)   
s_ht = qhz1 + qhz2*dfz + (qhz3 + qhz4*dfz)gam_str      # (4.E35)
a_t = alpha_str + s_ht      # (4.E34)
e_t = (qez1 + qez2*dfz + qez3*dfz^2)*(1+(qez4 + qez5*gam_str)*(2/pi)*atan(b_t*c_t*a_t))     # (4.E44)
# @assert e_t <= 1
d_t0 = fz*(r0/fz0p)*(qdz1 + qdz2*dfz)*(1-ppz1*dpi)*ltr      # (4.E42) TODO: velocity sign 
d_t = d_t0*(1+qdz3*abs(gam_str) + qdz4*gam_str^2)*zeta      # (4.E43)
# @assert c_t > 0
k_ya_p = k_ya + epsilon_K                              # (4.E39)
s_hf = s_hy + s_vy/k_ya_p                              # (4.E38)
a_r = alpha_str + s_hf                                 # (4.E37)
m_zr0 = d_r*cos(c_r*atan(b_r*a_r))*cos(alpha)          # (4.E36)
t_0 = d_t*cos(c_t*atan(b_t*a_t - e_t*(b_t*a_t-atan(b_t*a_t))))*cos(alpha)       # (4.E33) TODO: cos prime -> velocity implementation 
m_z0_p = -t_0*f_y0          # (4.E32)
m_z0 = m_z0_p + m_zr0       # (4.E31)

### fx

# longitudinal force (combined slip)
s_hxa = rhx1                    # (4.E57)
e_xa = rex1 + rex2*dfz          # (4.E56)
# @assert e_xa <= 1
c_xa = rcx1                     # (4.E55)
b_xa = (rbx1 + rbx3*gam_str^2)*cos(atan(rbx2*kappa))*lxal       # (4.E54)
# @assert b_xa > 0 
a_s = alpha_str + s_hxa     # (4.E53)
g_xa0 = cos(c_xa*atan(b_xa*s_hxa - e_xa*(b_xa*s_hxa - atan(b_xa*s_hxa))))  # (4.E52)
g_xa = cos(c_xa*atan(b_xa*a_s - e_xa*(b_xa*a_s - atan(b_xa*a_s))))/g_xa0    # (4.E51)
# @assert g_xa > 0 
f_x = g_xa*f_x0      # (4.E50)

# lateral force (combined slip)

d_vyk = mu_y*fz*(rvy1+rvy2*dfz+rvy3*gam_str)*cos(atan(rvy4*alpha_str)) * zeta
s_vyk = d_vyk*sin(rvy5*atan(rvy6*kappa)) * lvyk
s_hyk = rhy1 + rhy2*dfz
e_yk = rey1 + rey2*dfz
# @assert e_yk <= 1
c_yk = rcy1
b_yk = (rby1 + rby4*gam_str^2)*cos(atan(rby2*(alpha_str - rby3))) * lyk
# @assert b_yk > 0
k_s = kappa + s_hyk
g_yk0 = cos(c_yk*atan(b_yk*s_hyk - e_yk*(b_yk*s_hyk - atan(b_yk*s_hyk))))
g_yk = cos(c_yk*atan(b_yk*k_s - e_yk*(b_yk*k_s - atan(b_yk*k_s))))/g_yk0
# @assert g_yk > 0
f_y = g_yk*f_y0 + s_vyk

# overturning couple
# currently only gives 0, but it makes sense when looking at the inputs

part_1 = qsx1*lvmx - qsx2*gamma*(1+ppmx1*dpi) + qsx3*(f_y/fz0)
part_2_1 = qsx5*(atan(qsx6*(fz/fz0)))^2
part_2_2 = qsx7*gamma + qsx8*atan(qsx9*(f_y/fz0))
part_2 = qsx4*cos(part_2_1)*sin(part_2_2)
part_3 = qsx10*atan(qsx11*(fz/fz0))*gamma

m_x = r0*fz*(part_1 + part_2 + part_3)*lmx

# Rolling Resistance Moment

part_1 = qsy1 + qsy2*f_x/fz0 + qsy3*abs(vs_x/v0) + qsy4*(vs_x/v0)^4 + (qsy5 + qsy6*fz/fz0)*gamma^2
part_2 = (fz/fz0)^qsy7 * (p_i/p_io)^qsy8

m_y = fz*r0*part_1*part_2*lmy

# Aligning Torque (combined slip)
alpha_req = sqrt(a_r^2 + (k_xk/k_ya_p)^2*kappa^2)*sign(a_r)
alpha_teq = sqrt(a_t^2 + (k_xk/k_ya_p)^2*kappa^2)*sign(a_t)
s = r0 * (ssz1 + ssz2*(f_y/fz0p) + (ssz3 + ssz4*dfz)*gam_str)*ls
m_zr = d_r * cos(c_r*atan(b_r*alpha_req))*cosalpha_p
fy_p = g_yk * f_y0
b_alpha_t = b_t*alpha_teq
t = d_t*cos(c_t*atan(b_alpha_t - e_t*(b_alpha_t - atan(b_alpha_t)))) * cosalpha_p
m_zp = -t * fy_p
m_z = m_zp + m_zr + s * f_x

ϵ_c = 1e-8

output_vec = [
    f_x,
    f_y,
    m_x,
    m_y,
    m_z
]

constraints = [
    c_x ≳ 0 + ϵ_c
    d_x ≳ 0 + ϵ_c
    e_x ≲ 1 

    c_y ≳ 0 + ϵ_c
    e_y ≲ 1

    b_t ≳ 0 + ϵ_c
    e_t ≲ 1 
    c_t ≳ 0 + ϵ_c

    e_xa ≲ 1
    b_xa ≳ 0 + ϵ_c
    g_xa ≳ 0 + ϵ_c
    
    e_yk ≲ 1
    b_yk ≳ 0 + ϵ_c
    g_yk ≳ 0 + ϵ_c
]

end
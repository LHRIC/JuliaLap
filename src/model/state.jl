module StateVariables
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables begin
    (x(t))[1:3]
    (q(t))[1:4]
    (z(t))[1:4]
    (θ(t))[1:4]
    Θ(t)

    (v(t))[1:3]
    (q̇(t))[1:4]
    (ż(t))[1:4]
    (ω(t))[1:4]
    Ω(t)
end

export  x, q, z, θ, Θ, v, q̇, ż, ω, Ω

end
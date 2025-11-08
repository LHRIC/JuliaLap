export gen_kin_models
include("hdpt_parsing.jl")
include("kin_structs.jl")
include("kin_params.jl")
include("hdpt_vec.jl")
include("residuals.jl")
include("unsprung_coord.jl")
include("../rotation/rotation_transformations.jl")
using ComponentArrays
using SciMLBase
using DiffEqBase
using NonlinearSolve
using Plots
using ForwardDiff
using BenchmarkTools
using PreallocationTools
using LinearAlgebra
using SparseArrays
using Sparspak
using FiniteDiff
using Interpolations
using GridInterpolations
using StaticArrays
# using SciMLSensitivity

function gen_kin_models(file_name::String, steer_range, f_shock_range, r_shock_range)
    fl_range = "B8:E25"
    rl_range = "B33:E50"
    fl_dict, rl_dict = excel2dict(file_name, fl_range, rl_range)
    fl_array = dict2cvec(fl_dict, "FL")
    rl_array = dict2cvec(rl_dict, "RL")
    fl_fun, fl_ctrl0 = gen_corner(fl_array)
    rl_fun, rl_ctrl0 = gen_corner(rl_array)

    fl_T0 = SMatrix{4,4}(fl_fun(fl_ctrl0))
    fl_jac0 = SMatrix{6,2}(jac_wrapper(fl_fun, fl_ctrl0, fl_T0))
    fl_T_array = zeros(typeof(fl_T0),size(steer_range)...,size(f_shock_range)...)
    fl_jac_array = zeros(typeof(fl_jac0),size(steer_range)...,size(f_shock_range)...)

    rl_T0 = SMatrix{4,4}(fl_fun(fl_ctrl0))
    rl_jac0 = SMatrix{6,2}(jac_wrapper(fl_fun, fl_ctrl0, fl_T0))
    rl_T_array = zeros(typeof(rl_T0),size(r_shock_range)...)
    rl_jac_array = zeros(typeof(rl_jac0),size(r_shock_range)...)


    # TODO: Benchmark ForwardDiff and FiniteDiff jacobian calculations against eachother
    for i = eachindex(steer_range)
        for j = eachindex(f_shock_range)
                ctrl = [fl_ctrl0[1] + f_shock_range[j], fl_ctrl0[2] + steer_range[i]]
                T = fl_fun(ctrl)
                jac = jac_wrapper(fl_fun, fl_ctrl0, T)
                fl_T_array[i,j] = T
                fl_jac_array[i,j] = jac
        end
    end

    for i = eachindex(r_shock_range)
        T = rl_fun(rl_ctrl0)
        jac = jac_wrapper(rl_fun, rl_ctrl0, T)
        rl_T_array[i] = T
        rl_jac_array[i] = jac
    end
    fl_T_itp = Interpolations.interpolate((steer_range, f_shock_range), fl_T_array, Gridded(Linear()))
    fl_jac_itp = Interpolations.interpolate((steer_range, f_shock_range), fl_jac_array, Gridded(Linear()))
    rl_T_itp = Interpolations.interpolate((r_shock_range,), rl_T_array, Gridded(Linear()))
    rl_jac_itp = Interpolations.interpolate((r_shock_range,), rl_jac_array, Gridded(Linear()))
    return fl_T_itp, fl_jac_itp, rl_T_itp, rl_jac_itp
end

function jac_wrapper(fun, ctrl, T)
    # T_jac = FiniteDiff.finite_difference_jacobian(fun,ctrl)
    T_jac = ForwardDiff.jacobian(fun, ctrl)
    T_jac = reshape(T_jac,(size(T)...,2))
    jac = matrix_jac2velocity_jac(T,T_jac)
    return jac
end

function gen_corner(c_array)
    float_hdpts, fixed_hdpts = hdpt_vec(c_array)
    u0 = float_hdpts
    Rvec, Cvec = residual_vec(float_hdpts, fixed_hdpts)
    initial_shock = norm(float_hdpts[7,:] - fixed_hdpts[6,:])
    initial_steer = float_hdpts[4,2]
    ctrl = [initial_shock, initial_steer]
    
    contact_patch = c_array[:CP]

    unsprung_transform = unsprung_coords(float_hdpts, contact_patch)
    f = NonlinearFunction{false}(kin_fun)

    kinematic_problem = NonlinearProblem(f, u0, (fixed_hdpts, Rvec, Cvec, ctrl))

    function kin_NL(y)
        kinematic_problem = remake(kinematic_problem,u0=u0,p=(fixed_hdpts,Rvec,Cvec,y))
        # sol = solve(kinematic_problem, NewtonRaphson(), abstol=1e-12)
        sol = solve(kinematic_problem, NewtonRaphson())

        # println(sol.retcode)
        # println(norm(sol.resid))
        if !SciMLBase.successful_retcode(sol.retcode)
            println("Unsuccessful nonlinear solve with retcode: ", sol.retcode, " at ", y)
        end
        return unsprung_transform(sol.u)
    end
    return kin_NL, ctrl
end

function kin_fun!(R, float_hdpts, (fixed_hdpts, Rvec, Cvec, ctrl))
    for i in eachindex(Rvec)
        residual_fun = Rvec[i]
        R[i] = residual_fun(float_hdpts, fixed_hdpts)
    end
    for i in eachindex(Cvec)
        R[i+length(Rvec)] = Cvec[i](float_hdpts, fixed_hdpts, ctrl[i])
    end
    return
end

function kin_fun(float_hdpts, (fixed_hdpts, Rvec, Cvec, ctrl))
    T = eltype(float_hdpts)
    R = Vector{T}(undef,length(Rvec) + length(Cvec))
    for i in eachindex(Rvec)
        R[i] = Rvec[i](float_hdpts, fixed_hdpts)
    end
    for i in eachindex(Cvec)
        R[i+length(Rvec)] = Cvec[i](float_hdpts, fixed_hdpts, ctrl[i])
    end
    return R
end
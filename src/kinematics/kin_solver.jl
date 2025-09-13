include("hdpt_parsing.jl")
include("kin_structs.jl")
include("kin_params.jl")
include("hdpt_vec.jl")
include("residuals.jl")
include("unsprung_coord.jl")
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
# using SciMLSensitivity

function gen_kin_models(file_name::String)
    fl_range = "B8:E25"
    rl_range = "B33:E50"
    fl_dict, rl_dict = excel2dict(file_name, fl_range, rl_range)
    fl_array = dict2cvec(fl_dict, "FL")
    fl_fun, ctrl0 = gen_corner(fl_array)
    ctrl1 = ctrl0
    ctrl2 = [ctrl0[1] + 10, ctrl0[2]]
    # println(ctrl1)
    # println(ctrl2)
    # println(fl_fun(ctrl1))
    # println(fl_fun(ctrl2))

    pos = fl_fun(ctrl1)
    jac = FiniteDiff.finite_difference_jacobian(fl_fun,ctrl1,absstep=1e-12)
    return pos, jac
    # return pos, jac
end

function gen_corner(c_array)
    float_hdpts, fixed_hdpts = hdpt_vec(c_array)
    u0 = float_hdpts
    Rvec, Cvec = residual_vec(float_hdpts, fixed_hdpts)
    initial_shock = norm(float_hdpts[7,:] - fixed_hdpts[6,:])
    initial_steer = float_hdpts[4,2]
    ctrl = [initial_shock, initial_steer]
    
    contact_patch = c_array[:CP]
    
    #= Some stuff that can be used to check the kinematic model
        R = (zeros(Float64, size(u0)))
        jacobian_function = (R,x) -> kin_fun!(R,x,(fixed_hdpts, Rvec, Cvec, ctrl))
        initial_jac = ForwardDiff.jacobian(jacobian_function, R, u0)
        initial_sparse = sparse(initial_jac)
        rank_of_jac = rank(initial_jac)
        display(initial_sparse)

        sparsity_map = findall(!iszero, initial_jac)
        sparsity_map = initial_jac .!= 0
        
        kin_fun!(R, float_hdpts, (fixed_hdpts, Rvec, Cvec, ctrl))
        println("INITIAL RESIDUAL")
        println(R)
    =#

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




include("hdpt_parsing.jl")
include("kin_structs.jl")
include("kin_params.jl")
include("hdpt_vec.jl")
include("residuals.jl")
using ComponentArrays
using NonlinearSolve
using Plots
using ForwardDiff
using BenchmarkTools
using PreallocationTools
using LinearAlgebra
using SparseArrays
using Sparspak

function gen_kin_models(file_name::String)
    fl_range = "B8:E25"
    rl_range = "B33:E50"
    fl_dict, rl_dict = excel2dict(file_name, fl_range, rl_range)
    fl_array = dict2cvec(fl_dict, "FL")
    idk = gen_corner(fl_array)
    return idk
end

function gen_corner(c_array)
    float_hdpts, fixed_hdpts = hdpt_vec(c_array)
    u0 = float_hdpts
    # u0_cache = DiffCache(float_hdpts)
    Rvec, Cvec = residual_vec(float_hdpts, fixed_hdpts)
    p = (Rvec, Cvec)
    initial_shock = norm(float_hdpts[7,:] - fixed_hdpts[6,:])
    initial_steer = float_hdpts[4,2]
    ctrl = [initial_shock, initial_steer]
    R = (zeros(Float64, size(u0)))
    jacobian_function = (R,x) -> kin_fun!(R,x,(fixed_hdpts, Rvec, Cvec, ctrl))
    # initial_jac = zeros(Float64,(27,27))
    # println(size(initial_jac))
    initial_jac = ForwardDiff.jacobian(jacobian_function, R, u0)
    # println(rank(initial_jac))
    # println(cond(initial_jac))
    initial_sparse = sparse(initial_jac)
    # display(initial_sparse)
    # println(typeof(initial_sparse))
    # println(typeof(initial_jac))
    # sparsity_map = findall(!iszero, initial_jac)
    # sparsity_map = initial_jac .!= 0
    # println(sparsity_map)
    # println(initial_sparse)
    
    # kin_fun!(R, float_hdpts, (fixed_hdpts, Rvec, Cvec, ctrl))
    # println("INITIAL RESIDUAL")
    # println(R)
    
    # println(f)
    # f = NonlinearFunction{true}(kin_fun!; jac_prototype=initial_sparse)
    f = NonlinearFunction{true}(kin_fun!)
    kinematic_problem = NonlinearProblem(f, u0, (fixed_hdpts, Rvec, Cvec, ctrl))
    # kinematic_problem = NonlinearProblem(kin_fun!, u0, (fixed_hdpts, Rvec, Cvec, ctrl))
    solutions = []
    trace = TraceMinimal(1)
    # trace = TraceWithJacobianConditionNumber(100)
    # trace = TraceAll(100)
    # trace = TraceMinimal(;print_frequency=10000000, store_frequency=1)

    for i in 1:1
        ctrl[1] = initial_shock + rand()*20 - 10
        ctrl[2] = initial_steer + rand()*20 - 10
        remake(kinematic_problem,u0=u0,p=(fixed_hdpts,Rvec,Cvec,ctrl))
        # sol = solve(kinematic_problem,NewtonRaphson(),
        #     abstol=5e-4, show_trace=Val(false), trace_level=trace, store_trace=Val(false), maxiters=5000)
        sol = solve(kinematic_problem,NewtonRaphson(),
            abstol=1e-6)
        # push!(solutions,sol)
    end
    # println(size(u0))
    # R = reshape((zeros(Float64, size(u0))),:,1)

    # sol = kin_fun!(R, u0, (fixed_hdpts, Rvec, Cvec, ctrl))
    # sol = solutions[1]
    # residuals = [h.fnorm for h in sol.trace.history]
    # print(residuals)
    # p = plot(1:length(residuals), residuals, xscale=:log10, yscale=:log10)
    # p = plot(1:20, residuals[1:20], xscale=:log10, yscale=:log10)

    # display(p)
    for sol in solutions
        # println(sol.retcode)
    end
    return solutions
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
    R = []
    for i in eachindex(Rvec)
        residual_fun = Rvec[i]
        push!(R,residual_fun(float_hdpts, fixed_hdpts))
    end
    for i in eachindex(Cvec)
        push!(R,Cvec[i](float_hdpts, fixed_hdpts, ctrl[i]))
    end
    return R
end




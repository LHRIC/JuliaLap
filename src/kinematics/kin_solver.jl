include("hdpt_parsing.jl")
include("kin_structs.jl")
include("kin_params.jl")
using ComponentArrays
using NonlinearSolve
using Plots
using ForwardDiff
using BenchmarkTools
using PreallocationTools

function gen_kin_model(file_name::String)
    fl_range = "B8:E25"
    rl_range = "B33:E50"
    fl_dict, rl_dict = excel2dict(file_name, fl_range, rl_range)
    fl_array = dict2cvec(fl_dict, "FL")
    rl_array = dict2cvec(rl_dict, "RL")
    fl_shock_range = collect(LinRange(-10,10,21))
    fl_steer_range = collect(LinRange(-10,10,21))
    # fl_shock_range = 10.0
    # fl_steer_range = 10.0
    fl_sol = kin_solve(rl_array,fl_shock_range,fl_steer_range)
    return fl_sol
end

function kin_solve(c_array::ComponentArray,shock_range,steer_range)
    parameters::kinParameters = kin_params(c_array)

    initial_shock = parameters.shock_ctrl.length
    initial_steer = parameters.hdpts.float_hdpts[:IT][2]
    
    shock_space = shock_range .+ initial_shock
    steer_space = steer_range .+ initial_steer
    parameters.hdpts.float_hdpts.cvec = DiffCache(parameters.hdpts.float_hdpts.cvec)
    kinematic_problem = NonlinearProblem(kin_fun!, parameters.hdpts.float_hdpts.cvec, parameters)
    
    solution_array = [kin_remake!(shock,steer,kinematic_problem,parameters)
    for shock in shock_space, steer in steer_space]
    return solution_array
end

function kin_remake!(shock, steer, kinematic_problem, parameters)
    parameters.shock_ctrl.length = shock
    parameters.steer_ctrl.val = steer
    remake(kinematic_problem, p=parameters)
    sol = solve(kinematic_problem, abstol=1e-6)
    # sol = solve(kinematic_problem,abstol=1e-6, NewtonRaphson(autodiff = AutoFiniteDiff()))

    return sol.retcode
end

function kin_fun!(R, x, p::kinParameters)
    # if typeof(p.hdpts.float_hdpts.cvec) === typeof(x)
    p.hdpts.float_hdpts.cvec .= x
        # By modifying the float hardpoints in place we avoid slow allocations
    # else
        # p.hdpts.float_hdpts.cvec = x 
        #= On the first iteration, it will replace the float hardpoints 
        wrapper with the type specified by NonLinearSolve.jl =#
    # end
    for i in eachindex(p.residual_vector)
        R[i] = residual(p.residual_vector[i])
    end
        return nothing
end
include("../kinematics/hdpt_parsing.jl")
include("../kinematics/kin_structs.jl")
include("../kinematics/kin_params.jl")
include("../kinematics/hdpt_vec.jl")
include("../kinematics/residuals.jl")
include("../kinematics/unsprung_coord.jl")
include("../rotation/rotation_transformations.jl")
using ForwardDiff
using Plots
using SparseArrays
using Printf


file_name = "src/parameters/HDPT_Export-9-13.xlsx"
fl_range = "B8:E25"
rl_range = "B33:E50"
fl_dict, rl_dict = excel2dict(file_name, fl_range, rl_range)
fl_array = dict2cvec(fl_dict, "FL")

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

function hdpt_jac(c_array)
    float_hdpts, fixed_hdpts = hdpt_vec(c_array)
    u0 = float_hdpts
    Rvec, Cvec = residual_vec(float_hdpts, fixed_hdpts)
    initial_shock = norm(float_hdpts[7,:] - fixed_hdpts[6,:])
    initial_steer = float_hdpts[4,2]
    ctrl = [initial_shock, initial_steer]
    contact_patch = c_array[:CP]
    R = (zeros(Float64, size(u0)))
    
    jacobian_function = (x) -> kin_fun(x,(fixed_hdpts, Rvec, Cvec, ctrl))
    initial_jac = ForwardDiff.jacobian(jacobian_function, u0)
    return initial_jac
end

function reorder_xyz(A)
    cols = size(A,2)
    groups = div(cols,3)
    A_xyz = Array{eltype(A)}(undef,size(A))
    for i = 1:groups
        A_xyz[:,3*(i-1) + 1 : 3*(i-1) + 3] = hcat(A[:,i], A[:,groups + i], A[:,2*groups + i])
    end
    # A_xyz = [[A[:,i], A[:,groups+i], A[:,2*groups + i]] for i = 1:groups]
    return A_xyz
end

fl_jac = hdpt_jac(fl_array)
fl_cols = 1:size(fl_jac,2)
fl_jac_xyz = reorder_xyz(fl_jac)

xs = [string("x", i) for i = 1:size(fl_jac)[1]]
ys = reverse([string("r", i) for i = 1:size(fl_jac)[2]])

fl_heat = heatmap(xs, ys, reverse(fl_jac_xyz; dims=1), c = :balance, clims = (-0.01, 0.01))
dim1 = size(fl_jac_xyz,1)
dim2 = size(fl_jac_xyz,2)
for i = 1:dim1
    for j = 1:dim2
        val = fl_jac_xyz[i,j]
        if val == 0.0
            num_str = ""
        else
            num_str = @sprintf("%.1e", fl_jac_xyz[i,j])
        end
        annotate!(j - 0.5, dim2 - i + 0.5, text(num_str, :white, :center, 2))
    end
end

fl_sparse = sparse(fl_jac_xyz)
fl_sparse_original = sparse(fl_jac)
display(fl_sparse)
display(fl_sparse_original)
display(fl_heat)
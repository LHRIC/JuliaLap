using ComponentArrays
using LinearAlgebra

mutable struct hdptWrapper
    cvec::ComponentVector
end
Base.getindex(x::hdptWrapper,sym::Symbol) = getfield(x,1)[sym]
# Allows hdptWrapper to function the same as the component array while also not breaking references to it

struct Linkage
    hdpt1_array::hdptWrapper
    hdpt2_array::hdptWrapper
    hdpt1_symbol::Symbol
    hdpt2_symbol::Symbol
    initial_length::Float64
end

function Linkage(hdpt1_array, hdpt2_array, hdpt1_symbol, hdpt2_symbol)
    initial_length = norm(hdpt1_array[hdpt1_symbol] - hdpt2_array[hdpt2_symbol])
    return Linkage(hdpt1_array, hdpt2_array, hdpt1_symbol, hdpt2_symbol, initial_length)
end

struct Revolute
    origin_array::hdptWrapper
    point_array::hdptWrapper
    origin_symbol::Symbol
    point_symbol::Symbol
    axis::Vector{Float64}
end

mutable struct LinearActuator
    hdpt1_array::hdptWrapper
    hdpt2_array::hdptWrapper
    hdpt1_symbol::Symbol
    hdpt2_symbol::Symbol
    length::Real
end

mutable struct Forced
    hdpt_array::hdptWrapper
    hdtp_symbol::Symbol
    val::Real
    idx::Int
end

function residual(linkage::Linkage)
    pos1 = linkage.hdpt1_array[linkage.hdpt1_symbol]
    pos2 = linkage.hdpt2_array[linkage.hdpt2_symbol]
    diff = linkage.initial_length - norm((pos1 - pos2))
    return diff/linkage.initial_length
end

function residual(revolute::Revolute)
    origin_pos = revolute.origin_array[revolute.origin_symbol]
    point_pos = revolute.point_array[revolute.point_symbol]
    diff_vec = point_pos - origin_pos
    axis = revolute.axis
    return dot(diff_vec, axis)/norm(diff_vec)
end

function residual(linear::LinearActuator)
    pos1 = linear.hdpt1_array[linear.hdpt1_symbol]
    pos2 = linear.hdpt2_array[linear.hdpt2_symbol]
    diff = linear.length - norm((pos1 - pos2))
    return diff/linear.length
end

function residual(forced::Forced)
    return forced.val - forced.hdpt_array[forced.hdtp_symbol][forced.idx]
end

mutable struct kinHdpts
    float_hdpts::hdptWrapper
    fixed_hdpts::hdptWrapper
    # ctrl_hdpts::hdptWrapper
end

mutable struct kinParameters
    hdpts::kinHdpts
    residual_vector::Vector{Union{Linkage, Revolute, LinearActuator, Forced}}
    shock_ctrl::LinearActuator
    steer_ctrl::Forced
end

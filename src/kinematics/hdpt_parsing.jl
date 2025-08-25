using XLSX
using StaticArrays

function excel2dict(file_name, fl_range, rl_range)
    xlsx_file = XLSX.readxlsx(file_name)
    fl_array = xlsx_file["Hardpoint Configurator"][fl_range]
    rl_array = xlsx_file["Hardpoint Configurator"][rl_range]
    fl_dict = Dict{String, SVector{3, Float64}}(zip(fl_array[:,1], eachrow(fl_array[:,2:4])))
    rl_dict = Dict{String, SVector{3, Float64}}(zip(rl_array[:,1], eachrow(rl_array[:,2:4])))
    return fl_dict, rl_dict
end

function dict2cvec(dict, prefix)
    c_vec = ComponentVector(
    LO = dict[prefix * "LO"],
    LIF = dict[prefix * "LIF"],
    LIA = dict[prefix * "LIA"],
    UO = dict[prefix * "UO"],
    UIF = dict[prefix * "UIF"],
    UIA = dict[prefix * "UIA"],
    OT = dict[prefix * "OT"],
    IT = dict[prefix * "IT"],
    CP = dict[prefix * "CP"],
    WC = dict[prefix * "WC"],
    OP = dict[prefix * "OP"],
    IP = dict[prefix * "IP"],
    BC = dict[prefix * "BC"],
    SO = dict[prefix * "SO"],
    SI = dict[prefix * "SI"],
    BP = dict[prefix * "BP"],
    AP = dict[prefix * "AP"],
    BE = dict[prefix * "BE"]
        )
    return c_vec
end

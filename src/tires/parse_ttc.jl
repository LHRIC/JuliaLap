
using MAT
using CSV
using DataFrames
using Statistics
using Unitful
using Unitful: °C, °F

const _HAS_FFTW = let
    try
        @eval using FFTW
        true
    catch
        false
    end
end

const _NOISE_FFT_WARNED = Ref(false)

function unit_from_string(u::AbstractString)
    u = lowercase(strip(u))

    return if u == "kph"
        (1000/3600)u"m/s"
    elseif u == "rpm"
        (2π/60)u"rad/s"
    elseif u == "deg"
        (π/180)u"rad"
    elseif u == "deg c"
        u"°C"
    elseif u == "deg f"
        u"°F"
    elseif u == "cm"
        0.01u"m"
    elseif u == "kpa"
        1000u"Pa"
    elseif u == "n"
        u"N"
    elseif u == "nm" || u == "n-m"
        u"N*m"
    elseif u == "s" || u == "sec"
        u"s"
    elseif u == "kg"
        u"kg"
    else
        nothing
    end
end

_to_float(x::Missing) = nothing
_to_float(x::Unitful.Quantity) = ustrip(x)
_to_float(x::Number) = float(x)
_to_float(x) = nothing

function _numeric_signal(col)
    out = Float64[]
    for v in col
        fv = _to_float(v)
        fv === nothing && continue
        isfinite(fv) || continue
        push!(out, fv)
    end
    return out
end

function _noise_flags(x::Vector{Float64};
    min_samples::Int = 128,
    max_fft_n::Int = 2048,
    eps::Float64 = 1e-12,
    zero_mean_ratio_max::Float64 = 0.30,
    spectral_flatness_min::Float64 = 0.50,
    low_freq_share_max::Float64 = 0.20,
    lag1_autocorr_abs_max::Float64 = 0.15,
    near_constant_std::Float64 = 1e-9)

    n = length(x)
    n < min_samples && return nothing

    μ = mean(x)
    σ = std(x)
    (σ <= near_constant_std) && return (true, μ, σ, 1.0, 1.0, 1.0, "near-constant signal")

    x0 = x .- μ
    if !_HAS_FFTW && !_NOISE_FFT_WARNED[]
        println("[parse_ttc] Note: FFTW is not installed; using slower DFT fallback for noise checks.")
        _NOISE_FFT_WARNED[] = true
    end

    x_fft = n <= max_fft_n ? x0 : @view x0[1:max_fft_n]
    X = _HAS_FFTW ? FFTW.fft(x_fft) : _dft(x_fft)
    n_fft = length(X)
    half_n = fld(n_fft, 2)
    half_n < 2 && return nothing

    power = abs2.(X[2:half_n])
    pmean = mean(power)
    pmean <= eps && return (true, μ, σ, 1.0, 1.0, 1.0, "very low spectral power")

    spectral_flatness = exp(mean(log.(power .+ eps))) / (pmean + eps)
    low_band_end = max(1, floor(Int, 0.10 * length(power)))
    low_freq_share = sum(@view power[1:low_band_end]) / (sum(power) + eps)
    lag1 = cor(view(x0, 1:length(x0)-1), view(x0, 2:length(x0)))
    lag1 = isnan(lag1) ? 0.0 : lag1
    zero_mean_ratio = abs(μ) / (σ + eps)

    looks_white_zero = (zero_mean_ratio <= zero_mean_ratio_max &&
                        spectral_flatness >= spectral_flatness_min &&
                        low_freq_share <= low_freq_share_max &&
                        abs(lag1) <= lag1_autocorr_abs_max)

    reason = looks_white_zero ? "white-noise-like around zero" : ""
    return (looks_white_zero, μ, σ, spectral_flatness, low_freq_share, lag1, reason)
end

function _dft(x::AbstractVector{<:Real})
    n = length(x)
    out = Vector{ComplexF64}(undef, n)
    for k in 0:(n - 1)
        s = 0.0 + 0.0im
        for t in 0:(n - 1)
            s += x[t + 1] * cis(-2π * k * t / n)
        end
        out[k + 1] = s
    end
    return out
end

function check_noise_quality(df::DataFrame; cols = nothing, dataset_label = "dataset")
    cols_to_check = cols === nothing ? names(df) : cols
    suspects = NamedTuple[]

    for col in cols_to_check
        col ∈ names(df) || continue
        x = _numeric_signal(df[!, col])
        flags = _noise_flags(x)
        flags === nothing && continue

        (bad, μ, σ, flatness, low_share, lag1, reason) = flags
        bad || continue

        push!(suspects, (
            col = String(col),
            n = length(x),
            mean = μ,
            std = σ,
            spectral_flatness = flatness,
            low_freq_share = low_share,
            lag1_autocorr = lag1,
            reason = reason,
        ))
    end

    if !isempty(suspects)
        println("[parse_ttc] Warning: possible bad values in $(dataset_label):")
        for s in suspects
            println("  - $(s.col): $(s.reason) (n=$(s.n), mean=$(s.mean), std=$(s.std), flatness=$(s.spectral_flatness), low_freq_share=$(s.low_freq_share), lag1=$(s.lag1_autocorr))")
        end
    end

    return suspects
end

# Outputs: (DataFrame df, Dict dict)
# df - DataFrame where each column represents a value and each row represents a trial
# dict - Dictionary consisting of all unused metadata
function parse_ttc(filepath::String, wanted_cols = ["TSTO", "RE", "P", "AMBTMP", 
    "TSTC", "FY", "V", "NFX", "SA", "RST", "N", "ET", "SL", "TSTI", "MX", 
    "FZ", "RUN", "RL", "SR", "MZ", "NFY", "FX", "IA"]; check_noise::Bool = true, noise_cols = nothing)

    filetype = split(filepath, ".")[lastindex(split(filepath, "."))]
    df = nothing
    dict = nothing

    if (filetype == "mat")
        dict = matread(filepath)
        
        values = Vector{Vector}()
        for key in wanted_cols
            push!(values, vec(dict[key]))
            delete!(dict, key)
        end

        df = DataFrame(wanted_cols .=> values)
    elseif (filetype == "dat")

        meta = nothing
        keys = nothing
        units = nothing
        open(filepath, "r") do io
            meta = readline(io)
            keys = readline(io)
            units = readline(io)
        end

        parts = split(meta, ';')
        dict = Dict{String, Any}()
        for part in parts
            part = strip(part)
            isempty(part) && continue

            # Case 1: "Key: Value"
            if occursin(":", part)
                key, value = split(part, ":", limit=2)
                dict[strip(key)] = strip(value)
            
            # Case 2: "Key=Value"
            elseif occursin("=", part)
                key, value = split(part, "=", limit=2)
                dict[strip(key)] = strip(value)

            # Case 3: "ISO False" or "ID YXXXX"
            else
                tokens = split(part)
                if length(tokens) ≥ 2
                    key = tokens[1]
                    value = join(tokens[2:end], " ")
                    dict[key] = value
                end
            end
        end

        df = CSV.read(filepath, DataFrame; header=2, skipto=4, delim='\t')
        keys = split(keys, "\t")
        units = split(units, "\t")

        for (col, unit_str) in zip(keys, units)
            if col ∈ names(df)
                u = unit_from_string(unit_str)
                if u !== nothing
                    df[!, col] = df[!, col] .* u
                end
            end
        end
    end

    if check_noise && df !== nothing
        check_noise_quality(df; cols = noise_cols, dataset_label = filepath)
    end

    return df, dict
end

# for a list of ttc file
# df - a concatenated list of dataframes; added column for file ID
# dict - fileID => metadata dictionary
function parse_ttc_list(file_list, wanted_cols = ["TSTO", "RE", "P", "AMBTMP", 
    "TSTC", "FY", "V", "NFX", "SA", "RST", "N", "ET", "SL", "TSTI", "MX", 
    "FZ", "RUN", "RL", "SR", "MZ", "NFY", "FX", "IA"]; check_noise::Bool = true, noise_cols = nothing)

    df = DataFrame()
    dict = Dict{String, Dict}()

    for filepath in file_list
        temp_df, temp_dict = parse_ttc(filepath, wanted_cols; check_noise = check_noise, noise_cols = noise_cols)

        id = temp_dict["ID"]

        temp_df[!, :ID] = fill(id, nrow(temp_df))
        append!(df, temp_df)
        dict[id] = temp_dict
    end

    return df, dict
end

function unitful_to_float(df::DataFrame)
    out = deepcopy(df)

    for col in names(out)
        if eltype(out[!, col]) <: Unitful.Quantity
            out[!, col] = ustrip.(out[!, col])
        end
    end

    return out
end

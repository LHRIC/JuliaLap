const CATALOG_IGNORE_DIRS = Set([
    ".git",
    ".ipynb_checkpoints",
    ".venv",
    "__pycache__",
    "node_modules",
])

const JULIA_FUNCTION_PATTERN = r"^\s*function\s+([A-Za-z_][A-Za-z0-9_!]*)"
const PYTHON_FUNCTION_PATTERN = r"^\s*def\s+([A-Za-z_][A-Za-z0-9_]*)\s*\("

function repository_root()
    return normpath(joinpath(@__DIR__, "..", "..", ".."))
end

function relative_repo_path(path::AbstractString)
    return relpath(path, repository_root())
end

function detect_language(path::AbstractString)
    extension = lowercase(splitext(path)[2])
    if extension == ".jl"
        return "julia"
    elseif extension == ".py"
        return "python"
    elseif extension == ".ipynb"
        return "jupyter"
    else
        return "other"
    end
end

function make_script_record(path::AbstractString)
    language = detect_language(path)
    kind = language == "jupyter" ? "notebook" : "script"
    return Dict(
        "id" => string(uuid4()),
        "kind" => kind,
        "language" => language,
        "name" => basename(path),
        "path" => relative_repo_path(path),
        "line" => nothing,
    )
end

function parse_symbol_records(path::AbstractString)
    language = detect_language(path)
    pattern = language == "julia" ? JULIA_FUNCTION_PATTERN :
        language == "python" ? PYTHON_FUNCTION_PATTERN : nothing

    isnothing(pattern) && return Dict{String, Any}[]

    records = Dict{String, Any}[]
    for (line_number, line) in enumerate(eachline(path))
        match_result = match(pattern, line)
        isnothing(match_result) && continue
        push!(records, Dict(
            "id" => string(uuid4()),
            "kind" => "function",
            "language" => language,
            "name" => match_result.captures[1],
            "path" => relative_repo_path(path),
            "line" => line_number,
        ))
    end
    return records
end

function scan_repository(root::AbstractString=repository_root())
    records = Dict{String, Any}[]

    for (current_root, dirs, files) in walkdir(root)
        filter!(dir -> !(dir in CATALOG_IGNORE_DIRS), dirs)

        for file in sort(files)
            path = joinpath(current_root, file)
            extension = lowercase(splitext(file)[2])
            extension in (".jl", ".py", ".ipynb") || continue

            push!(records, make_script_record(path))
            append!(records, parse_symbol_records(path))
        end
    end

    sort!(records; by=record -> (
        get(record, "path", ""),
        get(record, "line", 0) === nothing ? 0 : get(record, "line", 0),
        get(record, "kind", ""),
        get(record, "name", ""),
    ))
    return records
end

function attach_action_metadata!(records, actions)
    for record in records
        record["runnable"] = false
        record["action_id"] = nothing
    end

    by_path = Dict{String, String}()
    by_function = Dict{Tuple{String, String}, String}()

    for (action_id, action) in actions
        path = get(action, "path", nothing)
        symbol_name = get(action, "symbol", nothing)

        if !isnothing(path) && isnothing(symbol_name)
            by_path[path] = action_id
        elseif !isnothing(path) && !isnothing(symbol_name)
            by_function[(path, symbol_name)] = action_id
        end
    end

    for record in records
        action_id = nothing
        if record["kind"] == "function"
            action_id = get(by_function, (record["path"], record["name"]), nothing)
        elseif record["kind"] in ("script", "notebook")
            action_id = get(by_path, record["path"], nothing)
        end

        if !isnothing(action_id)
            record["runnable"] = true
            record["action_id"] = action_id
        end
    end

    return records
end

function build_catalog(actions)
    records = scan_repository()
    return attach_action_metadata!(records, actions)
end

function filter_catalog(records; q=nothing, kind=nothing, language=nothing, runnable=nothing)
    lowered_query = isnothing(q) ? nothing : lowercase(String(q))

    return filter(records) do record
        if !isnothing(kind) && get(record, "kind", nothing) != kind
            return false
        end
        if !isnothing(language) && get(record, "language", nothing) != language
            return false
        end
        if !isnothing(runnable) && get(record, "runnable", nothing) != runnable
            return false
        end
        if isnothing(lowered_query)
            return true
        end

        haystack = lowercase(string(record["name"], " ", record["path"]))
        return occursin(lowered_query, haystack)
    end
end

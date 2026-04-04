const SCRIPT_RUNNERS = Dict(
    "julia" => `julia --project=.`,
    "python" => `python3`,
)

function flatten_array(value)
    return collect(vec(Array(value)))
end

function matrix_to_rows(value)
    matrix = Array(value)
    return [collect(matrix[row, :]) for row in axes(matrix, 1)]
end

function normalize_input_path(raw_path::AbstractString)
    root = repository_root()
    absolute = normpath(joinpath(root, raw_path))
    startswith(absolute, root) || error("Path escapes repository root.")
    isfile(absolute) || error("File does not exist: $(raw_path)")
    return absolute
end

function run_command_capture(command::Cmd; dir::AbstractString)
    return cd(dir) do
        stdout_buffer = IOBuffer()
        stderr_buffer = IOBuffer()
        process = run(pipeline(command, stdout=stdout_buffer, stderr=stderr_buffer); wait=false)
        wait(process)

        return Dict(
            "exit_code" => process.exitcode,
            "stdout" => String(take!(stdout_buffer)),
            "stderr" => String(take!(stderr_buffer)),
            "success" => process.exitcode == 0,
        )
    end
end

function run_kinematics_model(params)
    file = get(params, "file", "JuliaLap/src/parameters/HDPT_Export-9-13.xlsx")
    absolute = normalize_input_path(file)
    transform, jacobian = gen_kin_models(absolute)

    return Dict(
        "input_file" => relative_repo_path(absolute),
        "transform" => matrix_to_rows(transform),
        "jacobian" => matrix_to_rows(jacobian),
        "transform_flat" => flatten_array(transform),
        "jacobian_flat" => flatten_array(jacobian),
    )
end

function run_allowlisted_script(path::AbstractString)
    absolute = normalize_input_path(path)
    language = detect_language(absolute)
    haskey(SCRIPT_RUNNERS, language) || error("No runner configured for $(language) files.")

    runner = SCRIPT_RUNNERS[language]
    workdir = language == "julia" ? normpath(joinpath(repository_root(), "JuliaLap")) : dirname(absolute)
    script_path = relpath(absolute, workdir)
    command = `$runner $script_path`
    return run_command_capture(command; dir=workdir)
end

function build_action_registry()
    return Dict(
        "kinematics.front_model" => Dict(
            "id" => "kinematics.front_model",
            "kind" => "function",
            "runtime" => "julia",
            "description" => "Compute the front-corner kinematic transform and velocity Jacobian from an Excel hardpoint file.",
            "path" => "JuliaLap/src/kinematics/kin_solver.jl",
            "symbol" => "gen_kin_models",
            "parameters" => Dict(
                "file" => "Optional repository-relative path to an Excel hardpoint file.",
            ),
            "handler" => params -> run_kinematics_model(params),
        ),
        "script.julialap.kin_script" => Dict(
            "id" => "script.julialap.kin_script",
            "kind" => "script",
            "runtime" => "julia",
            "description" => "Run the existing Julia kinematics script in the JuliaLap project environment.",
            "path" => "JuliaLap/src/kin_script.jl",
            "parameters" => Dict(),
            "handler" => _ -> run_allowlisted_script("JuliaLap/src/kin_script.jl"),
        ),
    )
end

function describe_actions(actions)
    descriptions = Dict{String, Any}[]
    for action_id in sort(collect(keys(actions)))
        action = actions[action_id]
        push!(descriptions, Dict(
            "id" => action["id"],
            "kind" => action["kind"],
            "runtime" => action["runtime"],
            "description" => action["description"],
            "path" => get(action, "path", nothing),
            "symbol" => get(action, "symbol", nothing),
            "parameters" => action["parameters"],
        ))
    end
    return descriptions
end

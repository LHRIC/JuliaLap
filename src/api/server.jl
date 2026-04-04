function parse_bool(value)
    if isnothing(value)
        return nothing
    end
    lowered = lowercase(String(value))
    if lowered in ("1", "true", "yes")
        return true
    elseif lowered in ("0", "false", "no")
        return false
    else
        error("Invalid boolean value: $(value)")
    end
end

function json_response(status::Integer, payload)
    return Genie.Renderer.Json.json(payload; status=status)
end

function error_response(status::Integer, message::AbstractString)
    return json_response(status, Dict("ok" => false, "error" => message))
end

function parse_json_body()
    payload = Genie.Requests.jsonpayload()
    isnothing(payload) && return Dict{String, Any}()
    return Dict{String, Any}(pairs(payload))
end

function query_param(name::AbstractString)
    return Genie.Router.query(Symbol(name), nothing)
end

function register_routes(actions)
    catalog = build_catalog(actions)

    Genie.Router.route("/health", method=Genie.Router.GET) do
        try
            return json_response(200, Dict(
                "ok" => true,
                "service" => "LonghornAPI",
                "framework" => "Genie",
                "catalog_entries" => length(catalog),
                "action_count" => length(actions),
                "generated_at" => string(now()),
            ))
        catch error
            return error_response(500, sprint(showerror, error))
        end
    end

    Genie.Router.route("/api/v1/catalog", method=Genie.Router.GET) do
        try
            filtered = filter_catalog(
                catalog;
                q=query_param("q"),
                kind=query_param("kind"),
                language=query_param("language"),
                runnable=parse_bool(query_param("runnable")),
            )
            return json_response(200, Dict("ok" => true, "items" => filtered))
        catch error
            return error_response(500, sprint(showerror, error))
        end
    end

    Genie.Router.route("/api/v1/actions", method=Genie.Router.GET) do
        try
            return json_response(200, Dict("ok" => true, "actions" => describe_actions(actions)))
        catch error
            return error_response(500, sprint(showerror, error))
        end
    end

    Genie.Router.route("/api/v1/actions/:action_id", method=Genie.Router.GET) do
        try
            action_id = string(Genie.Router.params(:action_id))
            haskey(actions, action_id) || return error_response(404, "Unknown action: $(action_id)")
            action = actions[action_id]
            return json_response(200, Dict("ok" => true, "action" => first(describe_actions(Dict(action_id => action)))))
        catch error
            return error_response(500, sprint(showerror, error))
        end
    end

    Genie.Router.route("/api/v1/actions/:action_id/run", method=Genie.Router.POST) do
        try
            action_id = string(Genie.Router.params(:action_id))
            haskey(actions, action_id) || return error_response(404, "Unknown action: $(action_id)")
            action = actions[action_id]
            params = parse_json_body()
            result = action["handler"](params)
            return json_response(200, Dict("ok" => true, "action_id" => action_id, "result" => result))
        catch error
            return error_response(500, sprint(showerror, error))
        end
    end
end

function serve(; host::AbstractString="127.0.0.1", port::Integer=8080)
    actions = build_action_registry()
    register_routes(actions)
    Genie.config.run_as_server = true
    Genie.config.server_host = String(host)
    Genie.config.server_port = Int(port)
    return Genie.up(Int(port), String(host); open_browser=false, async=false)
end

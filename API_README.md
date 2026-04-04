# LonghornAPI

This adds a Genie-based Julia REST API on top of the repository so a web server or frontend can:

- query the repository for Julia and Python scripts/functions
- search for runnable entry points
- execute a small allowlist of approved actions

## Why an allowlist

Running arbitrary files from an HTTP request would make this machine a remote code execution target. The current server indexes the whole repo for discovery, but only runs actions you explicitly register in `src/api/actions.jl`.

## Start the server

```bash
julia --project=JuliaLap JuliaLap/bin/rest_api.jl
```

Optional environment variables:

- `LONGHORN_API_HOST`
- `LONGHORN_API_PORT`

## Endpoints

### `GET /health`

Basic status check.

### `GET /api/v1/catalog`

Returns indexed repository items.

Supported query params:

- `q`
- `kind` with values like `script`, `function`, `notebook`
- `language` with values like `julia`, `python`
- `runnable` with values `true` or `false`

### `GET /api/v1/actions`

Returns the allowlisted runnable actions.

### `GET /api/v1/actions/{id}`

Returns metadata for one action.

### `POST /api/v1/actions/{id}/run`

Runs an action with a JSON body.

Example:

```bash
curl -X POST http://127.0.0.1:8080/api/v1/actions/kinematics.front_model/run \
  -H 'Content-Type: application/json' \
  -d '{"file":"JuliaLap/src/parameters/HDPT_Export-9-13.xlsx"}'
```

## Current runnable actions

- `kinematics.front_model`
- `script.julialap.kin_script`

## Next steps

- add more wrappers around real simulation functions
- add authentication before exposing this outside localhost
- add async jobs and timeouts for long-running solves
- split legacy scripts into cleaner parameterized functions

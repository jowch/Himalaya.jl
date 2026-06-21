using HTTP, JSON3, DBInterface, Tables, Oxygen

# ─────────────────────────────────────────────────────────────────────────────
# Filesystem probes for the ingest funnel (spec §6.1). Read-only — no DB
# writes, no with_idempotency, no apply_event!, no SSE, no client_op_id,
# no pendingDeferreds participation. Queue-orthogonal by design.
# ─────────────────────────────────────────────────────────────────────────────

"""
    register_fs_routes!()

Filesystem probes for the ingest funnel (spec §6.1). Read-only, no DB writes.
- GET /api/fs/suggest?prefix=  → directory autocomplete for the picker.
- GET /api/fs/validate?path=   → cheap picker gate (exists + not-already-an-experiment).
- GET /api/fs/manifest?path=&{image,metadata,integration}_pattern=  → phase-1 file manifest.

`validate` and `manifest` are added in Tasks 1.4 / 1.5.
"""
function register_fs_routes!()
    @get "/api/fs/suggest" function(req::HTTP.Request)
        q      = HTTP.queryparams(HTTP.URI(req.target))
        prefix = get(q, "prefix", "")
        isempty(prefix) && return _json(200, Dict(:suggestions => String[]))
        dir  = isdir(prefix) ? prefix : dirname(prefix)
        base = isdir(prefix) ? "" : basename(prefix)
        isdir(dir) || return _json(200, Dict(:suggestions => String[]))
        kids = String[]
        for name in readdir(dir; sort = true)
            startswith(name, base) || continue
            full = joinpath(dir, name)
            isdir(full) && push!(kids, full)
            length(kids) >= 20 && break
        end
        _json(200, Dict(:suggestions => kids))
    end
    @get "/api/fs/validate" function(req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = get(q, "path", "")
        if isempty(path) || !isdir(path)
            return _json(200, Dict(:ok => false, :matched => 0, :scanned => 0,
                               :message => "path does not exist or is not a directory"))
        end
        dup = !isempty(DBInterface.execute(current_db(),
            "SELECT 1 FROM experiments WHERE data_dir = ? LIMIT 1", [path]) |> Tables.rowtable)
        if dup
            return _json(200, Dict(:ok => false, :matched => 0, :scanned => 0,
                               :message => "an experiment already uses this directory"))
        end
        scanned = count(!startswith("."), readdir(path))   # cheap; rich count is /manifest
        _json(200, Dict(:ok => true, :matched => scanned, :scanned => scanned, :message => nothing))
    end
    # manifest added in Task 1.5
end

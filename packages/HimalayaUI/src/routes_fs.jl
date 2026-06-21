using HTTP, JSON3, Oxygen

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
    # validate + manifest added in Tasks 1.4 / 1.5
end

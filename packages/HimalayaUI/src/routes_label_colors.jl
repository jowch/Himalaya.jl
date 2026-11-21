using HTTP, JSON3, DBInterface, Tables, Oxygen

const LABEL_COLORS = ("rose", "plum", "violet", "blue", "teal", "green", "ochre", "stone")
const LABEL_COLOR_ENTITY_ID = 0

function label_colors_list(db)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT label, color FROM label_colors ORDER BY label COLLATE NOCASE"))
    [Dict(:label => String(r.label), :color => String(r.color)) for r in rows]
end

function _label_color(db, label::String)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT color FROM label_colors WHERE label = ? COLLATE NOCASE", [label]))
    isempty(rows) ? nothing : String(rows[1].color)
end

function set_label_color!(db, req, label::String, color::Union{String,Nothing})
    _label_color(db, label) == color && return
    result = apply_event!(InTransaction(), db, req;
        kind = "tag_color_set", entity_type = "tag", entity_id = LABEL_COLOR_ENTITY_ID,
        payload = Dict(:label => label, :color => color))
    _enqueue_broadcast_from_result!(result, "tag_color_set", "tag", LABEL_COLOR_ENTITY_ID)
end

function register_label_color_routes!()
    @get "/api/label-colors" function(req::HTTP.Request)
        _json_ok(label_colors_list(current_db()))
    end

    @put "/api/label-colors" function(req::HTTP.Request)
        db = current_db()
        body = try
            json(req)
        catch
            return _json_error(400, "body must be JSON")
        end
        body isa AbstractDict || return _json_error(400, "body must be an object")
        label = normalize_label(get(body, "label", nothing))
        label === nothing && return _json_error(400, "label must be 1 to $MAX_LABEL_LENGTH characters on one line")
        any(r -> lowercase(r) == lowercase(label), ROLE_LABELS) &&
            return _json_error(400, "$(join(ROLE_LABELS, ", ")) have fixed colours")
        haskey(body, "color") || return _json_error(400, "color is required; send null to reset")
        raw = get(body, "color", nothing)
        (raw === nothing || (raw isa AbstractString && raw in LABEL_COLORS)) ||
            return _json_error(400, "color must be one of $(join(LABEL_COLORS, ", ")) or null")
        color = raw === nothing ? nothing : String(raw)

        return with_idempotency(db, req) do
            set_label_color!(db, req, label, color)
            _json_ok(label_colors_list(db))
        end
    end
end

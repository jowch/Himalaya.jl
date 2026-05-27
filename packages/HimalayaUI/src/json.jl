"""
    row_to_json(row; bool_keys = ())

Convert a NamedTuple row to a Symbol-keyed Dict suitable for JSON3 encoding.
Converts `missing` to `nothing`. Coerces any keys listed in `bool_keys` to Bool.
"""
function row_to_json(row::NamedTuple; bool_keys::Tuple = ())
    out = Dict{Symbol, Any}()
    for k in propertynames(row)
        v = getproperty(row, k)
        if v isa Missing
            out[k] = nothing
        elseif k in bool_keys
            out[k] = v != 0
        else
            out[k] = v
        end
    end
    out
end

"""
    rows_to_json(rows; bool_keys = ())

Convert a collection of rows (typically the output of `Tables.rowtable`).
"""
function rows_to_json(rows; bool_keys::Tuple = ())
    [row_to_json(r; bool_keys) for r in rows]
end

"""
    _json_error(status, msg; extra...)

Single-shot JSON error response with `Content-Type: application/json`.

Shared route helper — relocated here from `routes_comparisons.jl` when the
Compare routes were retired (I3.6 #177), since `routes_series.jl` reuses it.
"""
function _json_error(status::Int, msg::AbstractString; extra...)
    body = Dict{Symbol, Any}(:error => msg)
    for (k, v) in extra
        body[k] = v
    end
    HTTP.Response(status, ["Content-Type" => "application/json"],
                  JSON3.write(body))
end

"""
    _view_fields_error(body) -> Union{HTTP.Response, Nothing}

Type-guard the optional view-choice fields. Returns a 400 `HTTP.Response`
when a field is present with the wrong type, else `nothing`.
`view_show_peak_ticks` / `view_show_peak_labels` land in INTEGER columns
read back through `Bool(...)`; a non-boolean would otherwise throw
`InexactError` on a later GET (a 500 instead of a clean 400 at write
time). A present-but-null value is allowed — it resets to the per-tab
default (spec §6.4).

Shared route helper — relocated here from `routes_comparisons.jl` (I3.6 #177);
`routes_series.jl` reuses it.
"""
function _view_fields_error(body)
    if haskey(body, :view_grouping_mode) && body.view_grouping_mode !== nothing &&
            !(body.view_grouping_mode isa AbstractString)
        return _json_error(400, "view_grouping_mode must be a string")
    end
    for k in (:view_show_peak_ticks, :view_show_peak_labels)
        if haskey(body, k) && body[k] !== nothing && !(body[k] isa Bool)
            return _json_error(400, "$(k) must be a boolean")
        end
    end
    nothing
end

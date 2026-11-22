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
"""
function _json_error(status::Int, msg::AbstractString; extra...)
    body = Dict{Symbol, Any}(:error => msg)
    for (k, v) in extra
        body[k] = v
    end
    HTTP.Response(status, ["Content-Type" => "application/json"],
                  JSON3.write(body))
end

_json(status::Int, body) =
    HTTP.Response(status, ["Content-Type" => "application/json"], JSON3.write(body))

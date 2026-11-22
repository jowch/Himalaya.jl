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

function rows_to_json(rows; bool_keys::Tuple = ())
    [row_to_json(r; bool_keys) for r in rows]
end

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

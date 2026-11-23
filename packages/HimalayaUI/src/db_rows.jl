using SQLite, Tables

function _declared_column_type(decl::Ptr)
    decl == C_NULL && return Any
    t = uppercase(unsafe_string(decl))
    occursin("INT", t) && return Union{Missing, Int64}
    any(k -> occursin(k, t), ("CHAR", "CLOB", "TEXT")) && return Union{Missing, String}
    any(k -> occursin(k, t), ("REAL", "FLOA", "DOUB")) && return Union{Missing, Float64}
    Any
end

function sqlrows(q::SQLite.Query)
    names = q.names
    handle = SQLite._get_stmt_handle(q.stmt)
    declared = Type[_declared_column_type(SQLite.C.sqlite3_column_decltype(handle, i - 1)) for i in eachindex(names)]
    values = Vector{Any}[]
    for row in q
        push!(values, Any[Tables.getcolumn(row, i) for i in eachindex(names)])
    end
    types = Type[all(v -> v[i] isa T, values) ? T : Any for (i, T) in enumerate(declared)]
    NT = NamedTuple{Tuple(names), Tuple{types...}}
    rows = Vector{NT}(undef, length(values))
    for (k, v) in enumerate(values)
        rows[k] = NT(Tuple(v))
    end
    rows
end

import JSON

function to_json(x)::String
    io = IOBuffer()
    write_json(io, x)
    String(take!(io))
end

write_json(io::IOBuffer, x) = (write(io, JSON.json(x; sort_keys = false)); nothing)
write_json(io::IOBuffer, ::Nothing) = (write(io, "null"); nothing)
write_json(io::IOBuffer, ::Missing) = (write(io, "null"); nothing)
write_json(io::IOBuffer, x::Bool) = (write(io, x ? "true" : "false"); nothing)
write_json(io::IOBuffer, x::Integer) = (print(io, x); nothing)
write_json(io::IOBuffer, x::Base.IEEEFloat) = (write(io, _finite_json_number(x)); nothing)
write_json(io::IOBuffer, x::AbstractFloat) = write_json(io, Float64(x))
write_json(io::IOBuffer, x::Symbol) = _write_json_string(io, String(x))
write_json(io::IOBuffer, x::AbstractString) = _write_json_string(io, String(x))
write_json(io::IOBuffer, x::AbstractChar) = _write_json_string(io, string(x))

function _finite_json_number(x::Base.IEEEFloat)::String
    isfinite(x) || throw(ArgumentError("$x not allowed to be written in JSON spec"))
    Base.Ryu.writeshortest(x)
end

function write_json(io::IOBuffer, x::AbstractDict)
    write(io, UInt8('{'))
    first = true
    for (k, v) in x
        first || write(io, UInt8(','))
        first = false
        _write_json_key(io, k)
        write(io, UInt8(':'))
        write_json(io, v)
    end
    write(io, UInt8('}'))
    nothing
end

function write_json(io::IOBuffer, x::NamedTuple)
    write(io, UInt8('{'))
    for (i, k) in enumerate(keys(x))
        i == 1 || write(io, UInt8(','))
        _write_json_string(io, String(k))
        write(io, UInt8(':'))
        write_json(io, getfield(x, k))
    end
    write(io, UInt8('}'))
    nothing
end

function write_json(io::IOBuffer, x::Union{AbstractVector, Tuple, AbstractSet})
    write(io, UInt8('['))
    first = true
    for v in x
        first || write(io, UInt8(','))
        first = false
        write_json(io, v)
    end
    write(io, UInt8(']'))
    nothing
end

_write_json_key(io::IOBuffer, k::Symbol) = _write_json_string(io, String(k))
_write_json_key(io::IOBuffer, k::AbstractString) = _write_json_string(io, String(k))
_write_json_key(io::IOBuffer, k) = _write_json_string(io, string(k)::String)

const _JSON_SHORT_ESCAPES = Dict{UInt8, UInt8}(
    0x08 => UInt8('b'), 0x09 => UInt8('t'), 0x0a => UInt8('n'), 0x0c => UInt8('f'), 0x0d => UInt8('r'))

_needs_json_escape(b::UInt8) = b < 0x20 || b == 0x7f || b == UInt8('"') || b == UInt8('\\')

function _write_json_string(io::IOBuffer, s::String)
    write(io, UInt8('"'))
    for b in codeunits(s)
        _needs_json_escape(b) ? _write_json_escape(io, b) : write(io, b)
    end
    write(io, UInt8('"'))
    nothing
end

function _write_json_escape(io::IOBuffer, b::UInt8)
    write(io, UInt8('\\'))
    if b == UInt8('"') || b == UInt8('\\')
        write(io, b)
    elseif haskey(_JSON_SHORT_ESCAPES, b)
        write(io, _JSON_SHORT_ESCAPES[b])
    else
        write(io, "u00", string(b; base = 16, pad = 2))
    end
    nothing
end

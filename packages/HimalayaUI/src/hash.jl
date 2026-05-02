using SHA

"""
    hash_trace_file(path) -> String

SHA-256 over the bytes of the .dat file. Stable across Julia versions and
re-reads of the same file. Used to detect "is the trace input unchanged
since findpeaks last ran?"
"""
function hash_trace_file(path::AbstractString)::String
    bytes2hex(open(SHA.sha256, path))
end

"""
    hash_peak_set(eff::NamedTuple) -> String

Content hash of an effective peak set. Inputs are sorted by q, encoded as
consecutive (Float64, Float64) tuples in native byte order, and SHA-256ed.
Used to detect "is the indexpeaks input unchanged since indices were
computed?" Order-independent (sort key is q).
"""
function hash_peak_set(eff::NamedTuple)::String
    n = length(eff.q)
    buf = Vector{UInt8}(undef, 16n)
    perm = sortperm(eff.q)
    for (i, k) in enumerate(perm)
        q  = Float64(eff.q[k])
        sh = Float64(eff.sharpness[k])
        q_bytes = reinterpret(UInt8, [q])
        sh_bytes = reinterpret(UInt8, [sh])
        copyto!(buf, 16(i-1) + 1, q_bytes)
        copyto!(buf, 16(i-1) + 9, sh_bytes)
    end
    bytes2hex(SHA.sha256(buf))
end

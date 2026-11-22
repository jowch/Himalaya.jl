using SHA

function hash_trace_file(path::AbstractString)::String
    bytes2hex(open(SHA.sha256, path))
end

const PEAK_SET_HASH_CACHE_CAP = 1024
const _PEAK_SET_HASH_CACHE = Dict{Tuple{Vector{Float64}, Vector{Float64}}, String}()
const _PEAK_SET_HASH_CACHE_ORDER = Vector{Tuple{Vector{Float64}, Vector{Float64}}}()
const _PEAK_SET_HASH_CACHE_LOCK = ReentrantLock()

function _clear_peak_set_hash_cache!()
    lock(_PEAK_SET_HASH_CACHE_LOCK) do
        empty!(_PEAK_SET_HASH_CACHE)
        empty!(_PEAK_SET_HASH_CACHE_ORDER)
    end
end

function hash_peak_set(eff::NamedTuple)::String
    fast = eff.q isa Vector{Float64} && eff.sharpness isa Vector{Float64}
    qs = fast ? eff.q        : collect(Float64, eff.q)
    ss = fast ? eff.sharpness : collect(Float64, eff.sharpness)
    lookup_key = (qs, ss)

    lock(_PEAK_SET_HASH_CACHE_LOCK) do
        cached = get(_PEAK_SET_HASH_CACHE, lookup_key, nothing)
        cached !== nothing && return cached

        n = length(qs)
        buf = Vector{UInt8}(undef, 16n)
        perm = sortperm(qs)
        for (i, k) in enumerate(perm)
            q_bytes  = reinterpret(UInt8, [qs[k]])
            sh_bytes = reinterpret(UInt8, [ss[k]])
            copyto!(buf, 16(i-1) + 1, q_bytes)
            copyto!(buf, 16(i-1) + 9, sh_bytes)
        end
        result = bytes2hex(SHA.sha256(buf))

        store_key = fast ? (copy(qs), copy(ss)) : lookup_key

        if length(_PEAK_SET_HASH_CACHE_ORDER) >= PEAK_SET_HASH_CACHE_CAP
            oldest = popfirst!(_PEAK_SET_HASH_CACHE_ORDER)
            delete!(_PEAK_SET_HASH_CACHE, oldest)
        end
        _PEAK_SET_HASH_CACHE[store_key] = result
        push!(_PEAK_SET_HASH_CACHE_ORDER, store_key)
        result
    end
end

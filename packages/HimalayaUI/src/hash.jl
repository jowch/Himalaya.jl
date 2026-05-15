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

Process-local memoized: replay-as-rerun (frontend mutation queue) repeatedly
recomputes the same hash on every no-op replay. Cache key compares (q, sharpness)
with `==`, so collisions are impossible. Bounded to `PEAK_SET_HASH_CACHE_CAP`
entries via FIFO eviction.
"""
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
    # Fast path: when eff.q and eff.sharpness are already Vector{Float64}, the
    # Dict lookup compares tuples by `==` (content-equality on each Vector), so
    # we can probe with the caller-owned arrays without allocating. Only on a
    # miss do we copy — the cache then owns its keys, so external mutation of
    # eff.q can't corrupt them.
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

        # Store under a cache-owned key so external mutation can't corrupt it.
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

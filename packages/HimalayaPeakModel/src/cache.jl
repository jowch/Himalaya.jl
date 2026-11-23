const STALE_TMP_SECONDS = 3600
const CACHE_MAGIC = b"HPMFEAT1"
const CACHE_HEADER_BYTES = 8 + 3 * 4

"""
    DiskFeatureCache(dir, version, key)

Candidate features cached as raw files in `dir`. `key(e)` names an exposure's
trace (its content hash), or returns `nothing` to skip caching it; `version`
is the feature set's `feature_version`. A file name holds both, so a changed
trace or feature version is a new file and old files are never read.
Features are stored as Float32 and `build_dataset` rounds fresh features the
same way, so a cold and a warm run build identical datasets.
"""
struct DiskFeatureCache
    dir::String
    version::String
    key::Function
end

_version_tag(c::DiskFeatureCache) = bytes2hex(sha1(c.version))[1:12]

function cache_path(c::DiskFeatureCache, e)
    k = c.key(e)
    k === nothing ? nothing : joinpath(c.dir, string(k, "-", _version_tag(c), ".f32"))
end

"""File layout: magic, then Int32 `n` candidates, `p` features, `h` hybrid peaks, then `n·p` Float32 (column-major), `n` Int32 candidate indices, `h` Int32 hybrid indices."""
function write_features(path::AbstractString, idx::AbstractVector{<:Integer}, X::AbstractMatrix,
                        hybrid::AbstractVector{<:Integer})
    n, p = size(X)
    n == length(idx) || throw(DimensionMismatch("one feature row per candidate"))
    mkpath(dirname(path))
    tmp = string(path, ".", getpid(), ".", Threads.threadid(), ".tmp")
    open(tmp, "w") do io
        write(io, CACHE_MAGIC, Int32(n), Int32(p), Int32(length(hybrid)))
        write(io, Float32.(X), Int32.(idx), Int32.(hybrid))
    end
    mv(tmp, path; force = true)
    path
end

"""A cache is an optimisation: a directory that cannot be written leaves the file out."""
function try_write_features(path, idx, X, hybrid)
    try
        write_features(path, idx, X, hybrid)
    catch err
        err isa Union{SystemError, Base.IOError} || rethrow()
        nothing
    end
end

"""`(idx, X, hybrid)` from a cache file, or `nothing` if it is missing, malformed or for another feature count."""
function read_features(path::AbstractString, n_features::Integer)
    isfile(path) || return nothing
    try
        open(path) do io
            filesize(path) >= CACHE_HEADER_BYTES || return nothing
            read(io, 8) == CACHE_MAGIC || return nothing
            n, p, h = Int.(read!(io, Vector{Int32}(undef, 3)))
            (n >= 0 && h >= 0 && p == n_features) || return nothing
            filesize(path) == CACHE_HEADER_BYTES + 4 * (n * p + n + h) || return nothing
            X = read!(io, Matrix{Float32}(undef, n, p))
            idx = Int.(read!(io, Vector{Int32}(undef, n)))
            hyb = Int.(read!(io, Vector{Int32}(undef, h)))
            idx, Float64.(X), hyb
        end
    catch err
        err isa Union{SystemError, EOFError} || rethrow()
        nothing
    end
end

"""
Deletes the cache's own `.f32` files except those of this feature version for `keep` keys, and
`.tmp` files left by an interrupted write once they are over an hour old. Other files are left alone.
A file that cannot be removed (unwritable directory, or another process removed it first) is skipped.
Returns how many it removed.
"""
function prune!(c::DiskFeatureCache, keep)
    isdir(c.dir) || return 0
    wanted = Set(string(k, "-", _version_tag(c), ".f32") for k in keep)
    n = 0
    for f in readdir(c.dir)
        path = joinpath(c.dir, f)
        try
            stale = endswith(f, ".f32") ? !(f in wanted) :
                    endswith(f, ".tmp") && time() - mtime(path) > STALE_TMP_SECONDS
            stale && isfile(path) || continue
            rm(path); n += 1
        catch err
            err isa Union{SystemError, Base.IOError} || rethrow()
        end
    end
    n
end

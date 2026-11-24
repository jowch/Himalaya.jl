const UV_DIRENT_UNKNOWN = Cint(0)
const UV_DIRENT_DIR     = Cint(2)
const UV_DIRENT_LINK    = Cint(3)
const DIR_BATCH = 64

const _DirEntryT = Base.Filesystem.uv_dirent_t

function _with_fs_req(f)
    req = Libc.malloc(Base.Filesystem._sizeof_uv_fs)
    try
        f(req)
    finally
        Libc.free(req)
    end
end

function _scandir(dir::AbstractString)::Vector{Tuple{String,Cint}}
    _with_fs_req() do req
        err = @ccall gc_safe=true uv_fs_scandir(C_NULL::Ptr{Cvoid}, req::Ptr{Cvoid}, String(dir)::Cstring,
                                                0::Cint, C_NULL::Ptr{Cvoid})::Cint
        try
            err < 0 && return Tuple{String,Cint}[]
            ent = Ref{_DirEntryT}()
            out = Tuple{String,Cint}[]
            while ccall(:uv_fs_scandir_next, Cint, (Ptr{Cvoid}, Ptr{_DirEntryT}), req, ent) != Base.UV_EOF
                push!(out, (unsafe_string(ent[].name), ent[].typ))
            end
            out
        finally
            Base.Filesystem.uv_fs_req_cleanup(req)
        end
    end
end

list_names(dir::AbstractString)::Vector{String} = sort!(first.(_scandir(dir)))

_is_folder(dir, name, typ) = typ == UV_DIRENT_DIR ||
    ((typ == UV_DIRENT_LINK || typ == UV_DIRENT_UNKNOWN) && isdir(joinpath(dir, name)))

_folders_in(dir, entries) =
    sort!(String[name for (name, typ) in entries if !startswith(name, ".") && _is_folder(dir, name, typ)])

folder_names(dir::AbstractString)::Vector{String} = _folders_in(dir, _scandir(dir))

function folder_listing(dir::AbstractString, holds_image)
    entries = _scandir(dir)
    (folders = _folders_in(dir, entries), holds_image = any(e -> holds_image(first(e)), entries))
end

function _open_dir(req, dir::AbstractString)::Ptr{Cvoid}
    err = @ccall gc_safe=true uv_fs_opendir(C_NULL::Ptr{Cvoid}, req::Ptr{Cvoid}, String(dir)::Cstring,
                                            C_NULL::Ptr{Cvoid})::Cint
    handle = err < 0 ? C_NULL : ccall(:uv_fs_get_ptr, Ptr{Cvoid}, (Ptr{Cvoid},), req)
    Base.Filesystem.uv_fs_req_cleanup(req)
    handle
end

function _close_dir(req, handle::Ptr{Cvoid})
    @ccall gc_safe=true uv_fs_closedir(C_NULL::Ptr{Cvoid}, req::Ptr{Cvoid}, handle::Ptr{Cvoid},
                                       C_NULL::Ptr{Cvoid})::Cint
    Base.Filesystem.uv_fs_req_cleanup(req)
end

function _read_batch!(req, handle::Ptr{Cvoid}, buffer::Vector{_DirEntryT})::Vector{String}
    n = @ccall gc_safe=true uv_fs_readdir(C_NULL::Ptr{Cvoid}, req::Ptr{Cvoid}, handle::Ptr{Cvoid},
                                          C_NULL::Ptr{Cvoid})::Cint
    names = String[unsafe_string(buffer[i].name) for i in 1:max(n, 0)]
    Base.Filesystem.uv_fs_req_cleanup(req)
    names
end

function any_name(pred, dir::AbstractString)::Bool
    _with_fs_req() do req
        handle = _open_dir(req, dir)
        handle == C_NULL && return false
        buffer = Vector{_DirEntryT}(undef, DIR_BATCH)
        try
            GC.@preserve buffer begin
                unsafe_store!(Ptr{Ptr{_DirEntryT}}(handle), pointer(buffer))
                unsafe_store!(Ptr{Csize_t}(handle + sizeof(Ptr{Cvoid})), Csize_t(DIR_BATCH))
                while true
                    names = _read_batch!(req, handle, buffer)
                    isempty(names) && return false
                    any(pred, names) && return true
                end
            end
        finally
            _close_dir(req, handle)
        end
    end
end

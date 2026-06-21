using Test, HTTP, JSON3, SQLite
using HimalayaUI

# Filesystem probe routes (spec §6.1). Read-only — no DB writes, no event log,
# no SSE. The DB argument to with_inproc_routes is required by the harness but
# unused by /api/fs/* (these routes never touch the DB).

@testset "GET /api/fs/suggest lists child dirs of a prefix" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            root = mktempdir()
            mkpath(joinpath(root, "alpha"))
            mkpath(joinpath(root, "beta"))
            touch(joinpath(root, "afile.txt"))
            resp = call("GET", "/api/fs/suggest?prefix=$(HTTP.escapeuri(joinpath(root, "a")))")
            @test resp.status == 200
            s = JSON3.read(resp.body).suggestions
            @test any(endswith.(s, "alpha"))
            @test !any(endswith.(s, "afile.txt"))   # files excluded
            @test !any(endswith.(s, "beta"))         # prefix-filtered
        end
    end
end

@testset "GET /api/fs/suggest empty prefix returns empty list" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            resp = call("GET", "/api/fs/suggest")
            @test resp.status == 200
            s = JSON3.read(resp.body).suggestions
            @test length(s) == 0
        end
    end
end

@testset "GET /api/fs/suggest exact dir lists its subdirs" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            root = mktempdir()
            mkpath(joinpath(root, "one"))
            mkpath(joinpath(root, "two"))
            # prefix is the directory itself with trailing slash → isdir branch
            resp = call("GET", "/api/fs/suggest?prefix=$(HTTP.escapeuri(root * "/"))")
            @test resp.status == 200
            s = JSON3.read(resp.body).suggestions
            @test any(endswith.(s, "one"))
            @test any(endswith.(s, "two"))
        end
    end
end

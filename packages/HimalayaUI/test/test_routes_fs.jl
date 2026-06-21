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

@testset "GET /api/fs/validate gates the picker" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            good = mktempdir()
            touch(joinpath(good, "x.tif"))

            # existing dir → ok
            r1 = call("GET", "/api/fs/validate?path=$(HTTP.escapeuri(good))")
            @test r1.status == 200
            body1 = JSON3.read(r1.body)
            @test body1.ok == true
            @test body1.message === nothing

            # nonexistent path → not ok + message
            r2 = call("GET", "/api/fs/validate?path=$(HTTP.escapeuri(good * "_nope"))")
            @test r2.status == 200
            body2 = JSON3.read(r2.body)
            @test body2.ok == false
            @test body2.message !== nothing

            # already-an-experiment → not ok
            HimalayaUI.create_experiment!(db; name="e", path=good,
                data_dir=good, analysis_dir=good)
            r3 = call("GET", "/api/fs/validate?path=$(HTTP.escapeuri(good))")
            @test r3.status == 200
            @test JSON3.read(r3.body).ok == false
        end
    end
end

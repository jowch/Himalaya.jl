using Test, HTTP
using HimalayaUI

const _REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _DEPLOY_SWAP = joinpath(_REPO_ROOT, "scripts", "deploy_swap.sh")

function _fake_privileged_path(tmp)
    bin = mkpath(joinpath(tmp, "fakebin"))
    marker = joinpath(tmp, "privileged-calls")
    for name in ("sudo", "systemctl")
        tool = joinpath(bin, name)
        write(tool, "#!/usr/bin/env bash\necho \"$name \$*\" >> '$marker'\n")
        chmod(tool, 0o755)
    end
    (bin, marker)
end

function _stage_release(root, stamp, body)
    release = mkpath(joinpath(root, "dist.releases", stamp))
    write(joinpath(release, "index.html"), body)
    next = mkpath(joinpath(root, "build", "next"))
    write(joinpath(next, "himalaya.so"), "so-$stamp")
    write(joinpath(next, "julia_version"), "1.13.4")
    write(joinpath(next, "depot_path"), "/depot")
    release
end

function _run_swap(root, release; path_prefix)
    out = IOBuffer()
    err = IOBuffer()
    cmd = addenv(`bash $_DEPLOY_SWAP $(joinpath(root, "dist")) $release $(joinpath(root, "build", "next")) $(joinpath(root, "build")) himalaya`,
                 "PATH" => "$path_prefix:$(ENV["PATH"])")
    ok = success(pipeline(cmd; stdout=out, stderr=err))
    (ok = ok, out = String(take!(out)), err = String(take!(err)))
end

@testset "make deploy swaps without sudo and leaves the restart to the user" begin
    mktempdir() do tmp
        bin, marker = _fake_privileged_path(tmp)
        root = mkpath(joinpath(tmp, "checkout"))
        dist = joinpath(root, "dist")
        mkpath(dist)
        write(joinpath(dist, "index.html"), "old-ui")
        mkpath(joinpath(root, "build"))
        write(joinpath(root, "build", "himalaya.so"), "so-old")

        running_sysimage = open(joinpath(root, "build", "himalaya.so"))
        release = _stage_release(root, "20261124-010000", "new-ui")
        r = _run_swap(root, release; path_prefix=bin)

        @test r.ok
        @test !isfile(marker)
        @test occursin("sudo systemctl restart himalaya", r.out)
        @test islink(dist)
        @test readlink(dist) == joinpath("dist.releases", "20261124-010000")
        @test read(joinpath(dist, "index.html"), String) == "new-ui"
        @test read(joinpath(root, "dist.releases", "00000000-000000-before-20261124-010000", "index.html"), String) == "old-ui"
        @test read(joinpath(root, "build", "himalaya.so"), String) == "so-20261124-010000"
        @test read(running_sysimage, String) == "so-old"
        close(running_sysimage)
        @test read(joinpath(root, "build", "julia_version"), String) == "1.13.4"
        @test !isfile(joinpath(root, "build", "next", "himalaya.so"))
    end
end

@testset "make deploy keeps the newest frontend releases and the live one" begin
    mktempdir() do tmp
        bin, _ = _fake_privileged_path(tmp)
        root = mkpath(joinpath(tmp, "checkout"))
        mkpath(joinpath(root, "build"))
        stamps = ["20261124-0$(i)0000" for i in 1:5]
        for s in stamps
            @test _run_swap(root, _stage_release(root, s, "ui-$s"); path_prefix=bin).ok
        end
        @test sort(readdir(joinpath(root, "dist.releases"))) == stamps[3:5]
        @test read(joinpath(root, "dist", "index.html"), String) == "ui-$(stamps[5])"
    end
end

@testset "make deploy points the frontend back when the sysimage cannot move" begin
    mktempdir() do tmp
        bin, _ = _fake_privileged_path(tmp)
        root = mkpath(joinpath(tmp, "checkout"))
        mkpath(joinpath(root, "build"))
        @test _run_swap(root, _stage_release(root, "20261124-010000", "ui-1"); path_prefix=bin).ok
        release = _stage_release(root, "20261124-020000", "ui-2")
        rm(joinpath(root, "build", "next", "himalaya.so"))
        r = _run_swap(root, release; path_prefix=bin)
        @test !r.ok
        @test occursin("the old frontend is back", r.err)
        @test read(joinpath(root, "dist", "index.html"), String) == "ui-1"
        @test !isdir(release)
    end
end

@testset "make deploy never prunes the release that was live before it, even past failed deploys" begin
    mktempdir() do tmp
        bin, _ = _fake_privileged_path(tmp)
        root = mkpath(joinpath(tmp, "checkout"))
        mkpath(joinpath(root, "build"))
        @test _run_swap(root, _stage_release(root, "20261124-010000", "ui-1"); path_prefix=bin).ok
        for orphan in ("20261124-020000", "20261124-030000")
            mkpath(joinpath(root, "dist.releases", orphan))
        end
        @test _run_swap(root, _stage_release(root, "20261124-040000", "ui-4"); path_prefix=bin).ok
        @test read(joinpath(root, "dist.releases", "20261124-010000", "index.html"), String) == "ui-1"
        @test read(joinpath(root, "dist", "index.html"), String) == "ui-4"
    end
end

@testset "make deploy refuses to convert a plain dist folder over an earlier conversion" begin
    mktempdir() do tmp
        bin, _ = _fake_privileged_path(tmp)
        root = mkpath(joinpath(tmp, "checkout"))
        mkpath(joinpath(root, "build"))
        mkpath(joinpath(root, "dist"))
        write(joinpath(root, "dist", "index.html"), "plain")
        release = _stage_release(root, "20261124-010000", "new-ui")
        mkpath(joinpath(root, "dist.releases", "00000000-000000-before-20261124-010000"))
        r = _run_swap(root, release; path_prefix=bin)
        @test !r.ok
        @test !islink(joinpath(root, "dist"))
        @test read(joinpath(root, "dist", "index.html"), String) == "plain"
    end
end

@testset "make deploy's recipe calls neither sudo nor systemctl" begin
    recipe = read(addenv(`make -n -s -C $_REPO_ROOT deploy`, "HIMALAYA_ENV_FILE" => "/nonexistent"), String)
    @test occursin("deploy_swap.sh", recipe)
    @test !occursin(r"\bsudo\b", recipe)
    @test !occursin(r"\bsystemctl\b", recipe)
end

@testset "the server keeps serving the frontend release it started with" begin
    mktempdir() do tmp
        releases = mkpath(joinpath(tmp, "dist.releases"))
        for (name, body) in (("a", "release-a"), ("b", "release-b"))
            mkpath(joinpath(releases, name))
            write(joinpath(releases, name, "index.html"), "<!doctype html><body>$body</body>")
            write(joinpath(releases, name, "asset.js"), "asset-$body")
        end
        dist = joinpath(tmp, "dist")
        symlink(joinpath("dist.releases", "a"), dist)

        prev = get(ENV, "HIMALAYA_FRONTEND_DIST", nothing)
        ENV["HIMALAYA_FRONTEND_DIST"] = dist
        try
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                symlink(joinpath("dist.releases", "b"), dist * ".link")
                run(`mv -T $(dist * ".link") $dist`)
                @test read(joinpath(dist, "index.html"), String) == "<!doctype html><body>release-b</body>"
                @test occursin("release-a", String(HTTP.get("$base/workspace").body))
                @test String(HTTP.get("$base/asset.js").body) == "asset-release-a"
            end
        finally
            prev === nothing ? delete!(ENV, "HIMALAYA_FRONTEND_DIST") : (ENV["HIMALAYA_FRONTEND_DIST"] = prev)
        end
    end
end

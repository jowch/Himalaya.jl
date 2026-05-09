using Test, HTTP, SQLite, DBInterface
using HimalayaUI

# Helpers are available via runtests.jl include order — see test_routes_resolve.jl.

# Spec §3.2 — SPA catch-all serves index.html for any non-API, non-asset
# path so deep-link URLs like /inspect/exp/sample reach the frontend.
# /api/* always returns 404 to prevent unregistered API typos from being
# masked as 200 HTML responses.

@testset "SPA fallback" begin
    mktempdir() do tmp
        # Synthesize a minimal dist directory.
        dist = joinpath(tmp, "dist")
        mkpath(dist)
        index_html = "<!doctype html><html><body>shell</body></html>"
        write(joinpath(dist, "index.html"), index_html)
        write(joinpath(dist, "asset.png"), b"\x89PNG\r\n\x1a\nfake")

        # Point the server at the synthetic dist.
        prev = get(ENV, "HIMALAYA_FRONTEND_DIST", nothing)
        ENV["HIMALAYA_FRONTEND_DIST"] = dist
        try
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                @testset "single-segment unknown path returns index.html" begin
                    r = HTTP.get("$base/foo"; status_exception=false)
                    @test r.status == 200
                    @test occursin("shell", String(r.body))
                    @test HTTP.header(r, "Cache-Control") == "no-store"
                    @test occursin("text/html", HTTP.header(r, "Content-Type"))
                end

                @testset "multi-segment unknown path returns index.html (pins /** syntax)" begin
                    r = HTTP.get("$base/inspect/exp/sample"; status_exception=false)
                    @test r.status == 200
                    @test occursin("shell", String(r.body))
                end

                @testset "/api/* unregistered returns 404 (does not fall through)" begin
                    r = HTTP.get("$base/api/typo-not-registered"; status_exception=false)
                    @test r.status == 404
                end

                @testset "/api with no subpath returns 404 (regression for bare-api guard hole)" begin
                    r = HTTP.get("$base/api"; status_exception=false)
                    @test r.status == 404
                end

                @testset "asset path served by dynamicfiles, not catch-all" begin
                    r = HTTP.get("$base/asset.png"; status_exception=false)
                    @test r.status == 200
                    # dynamicfiles doesn't add no-store; the catch-all does.
                    @test HTTP.header(r, "Cache-Control") != "no-store"
                end
            end
        finally
            if prev === nothing
                delete!(ENV, "HIMALAYA_FRONTEND_DIST")
            else
                ENV["HIMALAYA_FRONTEND_DIST"] = prev
            end
        end
    end
end

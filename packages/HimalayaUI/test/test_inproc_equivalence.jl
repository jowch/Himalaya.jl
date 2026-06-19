# Differential equivalence harness: proves Oxygen.internalrequest (in-process)
# produces the SAME HTTP.Response as a real socket request, for the shapes most
# likely to diverge. This is the single load-bearing equivalence proof (spec
# §5.1). Keep it green permanently.
using Test, HTTP, SQLite, JSON3, DBInterface
using FileIO, ImageCore, TiffImages   # for the binary-PNG contract row
using HimalayaUI

if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

# Drop only transport-only headers; assert every other header matches so
# contract headers (Content-Type, Cache-Control, X-Image-*) cannot slip.
const _TRANSPORT_HEADERS = Set(lowercase.(
    ["Date", "Server", "Transfer-Encoding", "Connection", "Content-Length"]))

_sig_headers(r::HTTP.Response) = sort([
    lowercase(k) => v for (k, v) in r.headers
    if !(lowercase(k) in _TRANSPORT_HEADERS)])

# Normalize body to a String regardless of representation: the wire path returns
# bytes (Vector{UInt8}) parsed off the socket; the in-process path may leave the
# body as the raw String set by format_response! (misc.jl:402). copy(::String)
# would throw, so branch on the type.
_body(r::HTTP.Response) = r.body isa AbstractString ? String(r.body) : String(copy(r.body))

"Assert a wire response and an in-process response are equivalent."
function _assert_equiv(wire::HTTP.Response, inproc::HTTP.Response; label="")
    @test wire.status == inproc.status
    @test _body(wire) == _body(inproc)
    @test _sig_headers(wire) == _sig_headers(inproc)
end

# A fully-seeded fixture DB (its own temp dir) that both transports run against.
# NOTE: idempotency-replay rows below build their OWN per-transport DB so the
# second transport can't replay the first's cached row.
function _seed(dir)
    analysis_dir = joinpath(dir, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db = HimalayaUI.open_db(joinpath(dir, "h.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=dir,
        data_dir=joinpath(dir, "data"), analysis_dir=analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    (; db, exp_id, s_id, e_id)
end

@testset "in-process ≡ wire equivalence" begin

    @testset "GET 200 JSON (list)" begin
        mktempdir() do d
            fx = _seed(d)
            w = with_test_server(fx.db) do port, base
                HTTP.get("$base/api/experiments/$(fx.exp_id)/samples")
            end
            ip = with_inproc_routes(fx.db) do call
                call("GET", "/api/experiments/$(fx.exp_id)/samples")
            end
            _assert_equiv(w, ip; label="list")
        end
    end

    @testset "GET binary PNG image — exercises the X-Image-*/Cache-Control contract (B-2)" begin
        # Must seed a REAL image_path or the route 404s on BOTH transports and the
        # 200 PNG branch (the only code that emits image/png + X-Image-*/Cache-Control)
        # is never reached. TIFF-seeding mirrors test_routes_image.jl:4-19. Use the
        # FULL image (not thumb) so X-Image-Width/Height are present, and so the
        # response is recomputed deterministically each call (same source ⇒ identical
        # PNG bytes across transports).
        mktempdir() do d
            tiff = joinpath(d, "img.tiff")
            save(tiff, Gray.(rand(Float32, 512, 384)))
            db = HimalayaUI.open_db(joinpath(d, "h.db"))   # M0.4 predates M2; use open_db
            exp_id = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, image_path=tiff)
            path = "/api/exposures/$e_id/image?thumb=0"   # explicit full branch; exercises query parsing
            w  = with_test_server(db) do port, base; HTTP.get("$base$path"; status_exception=false) end
            ip = with_inproc_routes(db) do call; call("GET", path) end
            @test w.status == 200                          # guard: the contract branch was actually hit
            @test HTTP.header(w, "Content-Type") == "image/png"
            @test HTTP.header(w, "X-Image-Width") != ""    # full branch only
            _assert_equiv(w, ip; label="image-png")        # body bytes + X-Image-*/Cache-Control parity
        end
    end

    @testset "uncaught throw → 500 (same status both paths)" begin
        # Drive a route that raises (FK violation → SQLiteException). serialize=true
        # must wrap it into the same status the wire produces.
        mktempdir() do d
            fx = _seed(d)
            hdrs = ["Content-Type" => "application/json", "X-Username" => "alice"]
            bodyj = JSON3.write(Dict(:body => "orphan"))
            w = with_test_server(fx.db) do port, base
                HTTP.post("$base/api/samples/99999/messages"; body=bodyj, headers=hdrs, status_exception=false)
            end
            ip = with_inproc_routes(fx.db) do call
                call("POST", "/api/samples/99999/messages";
                     headers=hdrs, body=Vector{UInt8}(bodyj))
            end
            @test w.status >= 400
            _assert_equiv(w, ip; label="fk-throw")
        end
    end

    @testset "DELETE verb parity (404 on missing id, no body)" begin
        # Exercises the DELETE method + a no-body request + an error-response
        # body, asserting parity. Uses a REAL delete route (/api/peaks/{id},
        # routes_peaks.jl:348) with a nonexistent id so no mutation/setup is
        # needed and both transports see the identical not-found response.
        mktempdir() do d
            fx = _seed(d)
            hdrs = ["X-Username" => "alice"]
            w  = with_test_server(fx.db) do port, base
                HTTP.request("DELETE", "$base/api/peaks/99999", hdrs; status_exception=false)
            end
            ip = with_inproc_routes(fx.db) do call
                call("DELETE", "/api/peaks/99999"; headers=hdrs)
            end
            _assert_equiv(w, ip; label="delete-404")
        end
    end

    @testset "idempotency replay: identical body, per-transport DB" begin
        # Each transport gets its OWN seeded DB + its OWN op_id, so neither
        # replays the other's cached row (idempotency.jl:108-110).
        function _replay(call_with, opid)
            hdrs = ["Content-Type"=>"application/json", "X-Username"=>"alice",
                    "X-Client-Id"=>"tab-1", "X-Client-Op-Id"=>opid]
            bodyj = JSON3.write(Dict(:body => "hello"))
            call_with(hdrs, bodyj)
        end
        wire_bodies = mktempdir() do d
            fx = _seed(d)
            with_test_server(fx.db) do port, base
                r1 = _replay((h,b)->HTTP.post("$base/api/samples/$(fx.s_id)/messages"; body=b, headers=h), "op-wire")
                r2 = _replay((h,b)->HTTP.post("$base/api/samples/$(fx.s_id)/messages"; body=b, headers=h), "op-wire")
                (r1.status, _body(r1), _body(r2))
            end
        end
        inproc_bodies = mktempdir() do d
            fx = _seed(d)
            with_inproc_routes(fx.db) do call
                r1 = _replay((h,b)->call("POST", "/api/samples/$(fx.s_id)/messages"; headers=h, body=Vector{UInt8}(b)), "op-ip")
                r2 = _replay((h,b)->call("POST", "/api/samples/$(fx.s_id)/messages"; headers=h, body=Vector{UInt8}(b)), "op-ip")
                (r1.status, _body(r1), _body(r2))
            end
        end
        @test wire_bodies[1] == inproc_bodies[1]          # same status
        @test wire_bodies[2] == wire_bodies[3]            # wire: replay identical
        @test inproc_bodies[2] == inproc_bodies[3]        # in-process: replay identical
    end

    @testset "GET numeric-array JSON body (trace) — serialization parity" begin
        # The trace route returns numeric arrays (q/I) — the shape most prone to
        # JSON serialization divergence (float precision, key order). routes_trace.jl:4.
        mktempdir() do d
            fx = _seed(d)
            w  = with_test_server(fx.db) do port, base; HTTP.get("$base/api/exposures/$(fx.e_id)/trace"; status_exception=false) end
            ip = with_inproc_routes(fx.db) do call; call("GET", "/api/exposures/$(fx.e_id)/trace") end
            @test w.status == 200
            _assert_equiv(w, ip; label="trace")
        end
    end

    @testset "SSE broadcast fans out identically in-process" begin
        # The wire path and in-process path must both flush the post-commit
        # broadcast to a fake subscriber. Push a (pending=Channel,) sub directly
        # (mirrors test_idempotency_replay_invariant.jl::_capture_sse_during) and
        # assert exactly one frame on each transport. NOTE: in-process flush is
        # synchronous (events.jl _flush_post_commit_broadcasts! runs in-task), so
        # no sleep is needed; the wire path needs a short drain.
        function _count_frames(seed_op)
            pending = Channel{String}(64)
            sub = (pending = pending,)
            lock(HimalayaUI.SSE_LOCK) do; push!(HimalayaUI.SSE_SUBSCRIBERS[], sub); end
            try
                seed_op()
            finally
                lock(HimalayaUI.SSE_LOCK) do
                    filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
                end
                close(pending)
            end
            count(f -> !startswith(f, ":") && occursin("post_message", f), collect(pending))
        end
        hdrs = ["Content-Type"=>"application/json", "X-Username"=>"alice"]
        bodyj = JSON3.write(Dict(:body => "hi"))
        wire_n = mktempdir() do d
            fx = _seed(d)
            _count_frames() do
                with_test_server(fx.db) do port, base
                    HTTP.post("$base/api/samples/$(fx.s_id)/messages"; body=bodyj, headers=hdrs)
                    sleep(0.3)   # wire broadcast fires off-task; allow the drain
                end
            end
        end
        inproc_n = mktempdir() do d
            fx = _seed(d)
            _count_frames() do
                with_inproc_routes(fx.db) do call
                    call("POST", "/api/samples/$(fx.s_id)/messages"; headers=hdrs, body=Vector{UInt8}(bodyj))
                end  # synchronous flush — frame already on the channel
            end
        end
        @test wire_n == 1
        @test inproc_n == 1
    end
end

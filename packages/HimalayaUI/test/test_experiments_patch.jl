# packages/HimalayaUI/test/test_experiments_patch.jl
using Test, JSON3, HimalayaUI
if !isdefined(@__MODULE__, :with_inproc_routes)
    include("test_http.jl")
end
using HimalayaUI: with_inproc_routes, open_prepared_clone

@testset "PATCH /api/experiments/:id — widened fields" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            # seed a minimal experiment row
            body_create = Vector{UInt8}(JSON3.write(Dict("path" => tmp, "name" => "test")))
            r = call("POST", "/api/experiments"; headers=["Content-Type"=>"application/json","X-Username"=>"alice"], body=body_create)
            exp_id = JSON3.read(r.body)["id"]

            # name is patchable (not in _READONLY_FIELDS)
            r2 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("name" => "renamed"))))
            @test r2.status == 200
            r2b = call("GET", "/api/experiments/$exp_id")
            @test JSON3.read(r2b.body)["name"] == "renamed"

            # description is patchable
            r3 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("description" => "AgBe SAXS run"))))
            @test r3.status == 200
            r3b = call("GET", "/api/experiments/$exp_id")
            @test JSON3.read(r3b.body)["description"] == "AgBe SAXS run"

            # image_pattern is patchable AND resets scan_signature
            r4 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("image_pattern" => "*.tif"))))
            @test r4.status == 200
            r4b = call("GET", "/api/experiments/$exp_id")
            exp4 = JSON3.read(r4b.body)
            @test exp4["image_pattern"] == "*.tif"
            @test ismissing(exp4["scan_signature"]) || exp4["scan_signature"] === nothing

            # data_dir stays read-only (400)
            r5 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("data_dir" => "/other"))))
            @test r5.status == 400
        end
    end
end

# packages/HimalayaUI/test/test_experiments_patch.jl
using Test, JSON3, HimalayaUI, Tables, DBInterface
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

            r3c = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}("{\"description\":null}"))
            @test r3c.status == 200
            @test JSON3.read(r3c.body)["description"] === nothing
            r3d = call("GET", "/api/experiments/$exp_id")
            @test JSON3.read(r3d.body)["description"] === nothing

            r4 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("image_pattern" => "*.tif"))))
            @test r4.status == 200
            r4b = call("GET", "/api/experiments/$exp_id")
            exp4 = JSON3.read(r4b.body)
            @test exp4["image_pattern"] == "*.tif"
            @test !haskey(exp4, "scan_signature")

            # data_dir stays read-only (400)
            r5 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("data_dir" => "/other"))))
            @test r5.status == 400
        end
    end
end

@testset "PATCH /api/experiments/:id logs and broadcasts the change" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        sub = (pending = Channel{String}(64),)
        hdrs = ["Content-Type" => "application/json", "X-Username" => "alice"]
        with_inproc_routes(db) do call
            r = call("POST", "/api/experiments"; headers=hdrs, body=Vector{UInt8}(JSON3.write(Dict("path" => tmp, "name" => "test"))))
            exp_id = JSON3.read(r.body)["id"]
            lock(HimalayaUI.SSE_LOCK) do
                push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
            end
            @test call("PATCH", "/api/experiments/$exp_id"; headers=hdrs,
                body=Vector{UInt8}(JSON3.write(Dict("name" => "  renamed  ")))).status == 200
            @test call("PATCH", "/api/experiments/$exp_id"; headers=hdrs,
                body=Vector{UInt8}("{\"description\":null}")).status == 200
            for bad in ("{\"name\":\"\"}", "{\"name\":\"   \"}", "{\"name\":3}")
                @test call("PATCH", "/api/experiments/$exp_id"; headers=hdrs, body=Vector{UInt8}(bad)).status == 400
            end
            @test JSON3.read(call("GET", "/api/experiments/$exp_id").body)["name"] == "renamed"

            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT entity_type, entity_id, payload FROM user_actions WHERE action = 'experiment_updated' ORDER BY id"))
            @test length(rows) == 2
            @test all(r -> r.entity_type == "experiment" && r.entity_id == exp_id, rows)
            @test JSON3.read(rows[1].payload)["name"] == "renamed"
            @test JSON3.read(rows[2].payload)["description"] === nothing
        end
        frames = String[]
        while isready(sub.pending)
            push!(frames, take!(sub.pending))
        end
        updated = filter(f -> occursin("experiment_updated", f), frames)
        @test length(updated) == 2
        @test occursin("\"name\":\"renamed\"", updated[1])
        @test occursin("\"description\":null", updated[2])
    end
end

using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

const _COLOR_HEADERS = ["Content-Type" => "application/json", "X-Username" => "alice"]

_set_color(call, label, color) = call("PUT", "/api/label-colors"; headers=_COLOR_HEADERS,
    body=Vector{UInt8}(JSON3.write(Dict(:label => label, :color => color))))

function _colors(call)
    r = call("GET", "/api/label-colors")
    @test r.status == 200
    Dict(String(x.label) => String(x.color) for x in JSON3.read(String(r.body)))
end

_color_rows(db) = Dict(String(r.label) => String(r.color) for r in
    Tables.rowtable(DBInterface.execute(db, "SELECT label, color FROM label_colors")))

@testset "label colors" begin
    @testset "a colour belongs to the label name everywhere, ignoring case" begin
        db = open_prepared_clone(mktempdir())
        with_inproc_routes(db) do call
            @test _colors(call) == Dict()
            r = _set_color(call, "37C", "teal")
            @test r.status == 200
            @test Dict(String(x.label) => String(x.color) for x in JSON3.read(String(r.body))) == Dict("37C" => "teal")
            _set_color(call, "37c", "rose")
            @test _colors(call) == Dict("37C" => "rose")
        end
    end

    @testset "null resets a label to its default colour" begin
        db = open_prepared_clone(mktempdir())
        with_inproc_routes(db) do call
            _set_color(call, "Buffer", "rose")
            @test _set_color(call, "buffer", nothing).status == 200
            @test _colors(call) == Dict()
        end
    end

    @testset "the eight colours match the frontend's list, in order" begin
        api = read(joinpath(@__DIR__, "..", "frontend", "src", "api.ts"), String)
        m = match(r"export const LABEL_COLORS = \[([^\]]*)\]", api)
        @test m !== nothing
        @test Tuple(String(x.captures[1]) for x in eachmatch(r"\"([a-z]+)\"", m.captures[1])) == HimalayaUI.LABEL_COLORS
        for c in HimalayaUI.LABEL_COLORS
            db = open_prepared_clone(mktempdir())
            with_inproc_routes(db) do call
                @test _set_color(call, "37C", c).status == 200
            end
        end
    end

    @testset "invalid requests are refused" begin
        db = open_prepared_clone(mktempdir())
        with_inproc_routes(db) do call
            @test _set_color(call, "37C", "magenta").status == 400
            @test _set_color(call, "37C", 3).status == 400
            @test _set_color(call, "   ", "teal").status == 400
            @test call("PUT", "/api/label-colors"; headers=_COLOR_HEADERS, body=Vector{UInt8}("[]")).status == 400
            @test call("PUT", "/api/label-colors"; headers=_COLOR_HEADERS, body=Vector{UInt8}("{\"label\":\"37C\"}")).status == 400
        end
        @test _color_rows(db) == Dict()
    end

    @testset "setting the same colour again writes no event" begin
        db = open_prepared_clone(mktempdir())
        with_inproc_routes(db) do call
            _set_color(call, "37C", "teal")
            _set_color(call, "37C", "teal")
            _set_color(call, "never-set", nothing)
        end
        @test length(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE action = 'tag_color_set'"))) == 1
    end

    @testset "colours rebuild from the log" begin
        db = open_prepared_clone(mktempdir())
        with_inproc_routes(db) do call
            _set_color(call, "37C", "teal")
            _set_color(call, "Buffer", "rose")
            _set_color(call, "retake", "rose")
            _set_color(call, "retake", nothing)
            _set_color(call, "37C", "stone")
        end
        live = _color_rows(db)
        @test live == Dict("37C" => "stone", "Buffer" => "rose")
        DBInterface.execute(db, "DELETE FROM label_colors")
        HimalayaUI.rebuild_views_from_log!(db, HimalayaUI.LABEL_COLOR_ENTITY_ID; entity_type = "tag")
        @test _color_rows(db) == live
    end

    @testset "a colour change broadcasts one SSE frame" begin
        db = open_prepared_clone(mktempdir())
        sub = (pending = Channel{String}(64),)
        with_inproc_routes(db) do call
            lock(HimalayaUI.SSE_LOCK) do
                push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
            end
            _set_color(call, "37C", "teal")
        end
        frames = String[]
        while isready(sub.pending)
            push!(frames, take!(sub.pending))
        end
        set = filter(f -> occursin("tag_color_set", f), frames)
        @test length(set) == 1
        @test occursin("\"color\":\"teal\"", only(set))
    end
    @testset "open_db renames an old tag_colors table and keeps its colours" begin
        path = joinpath(mktempdir(), "old.db")
        db = HimalayaUI.open_db(path)
        DBInterface.execute(db, "DROP TABLE label_colors")
        DBInterface.execute(db, "CREATE TABLE tag_colors (label TEXT NOT NULL, color TEXT NOT NULL)")
        DBInterface.execute(db, "CREATE UNIQUE INDEX idx_tag_colors_label ON tag_colors(label COLLATE NOCASE)")
        DBInterface.execute(db, "INSERT INTO tag_colors VALUES ('37C', 'teal')")
        close(db)
        db = HimalayaUI.open_db(path)
        names(t) = [String(r.name) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type = ?", [t]))]
        @test "tag_colors" ∉ names("table")
        @test "idx_tag_colors_label" ∉ names("index")
        @test "idx_label_colors_label" ∈ names("index")
        @test _color_rows(db) == Dict("37C" => "teal")
        close(db)
    end
end

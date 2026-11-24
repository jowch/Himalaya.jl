using Test
import JSON
using HimalayaUI: row_to_json, rows_to_json

@testset "row_to_json" begin
    nt = (id = 1, name = "foo", active = 1, notes = missing)
    j  = row_to_json(nt)
    @test j[:id]     == 1
    @test j[:name]   == "foo"
    @test j[:active] == 1
    @test j[:notes]  === nothing

    j2 = row_to_json(nt; bool_keys = (:active,))
    @test j2[:active] === true
end

@testset "rows_to_json" begin
    rows = [(id = 1, q = 0.1), (id = 2, q = 0.2)]
    @test rows_to_json(rows) == [Dict(:id => 1, :q => 0.1),
                                  Dict(:id => 2, :q => 0.2)]
end

@testset "to_json writes the bytes JSON3 wrote" begin
    v = (big = 2^62, ok = true, none = nothing, n = Int64[], id = 7, phase = :Pn3m,
         name = "a\"b\\c\n\t\u0001\u007f é/", row = (peak_id = 3, q = 0.1f0), sub = SubString("xyz", 2, 3),
         q = [0.1, 1/3, 1.0e-7, 1.23456789e8, -0.0, 1.0e22])
    @test HimalayaUI.to_json(v) ==
          """{"big":4611686018427387904,"ok":true,"none":null,"n":[],"id":7,"phase":"Pn3m","name":"a\\"b\\\\c\\n\\t\\u0001\\u007f é/","row":{"peak_id":3,"q":0.1},"sub":"yz","q":[0.1,0.3333333333333333,1.0e-7,1.23456789e8,-0.0,1.0e22]}"""
    @test HimalayaUI.to_json([(id = 1, q = 0.5, excluded = false, note = missing)]) ==
          """[{"id":1,"q":0.5,"excluded":false,"note":null}]"""
end

@testset "to_json agrees with JSON.json in iteration order" begin
    every_char = join(Char(c) for c in 0:0x3000 if isvalid(Char, c))
    values = Any[
        every_char,
        Dict{Symbol, Any}(:a => 1, :b => [1.5, nothing], :c => Dict{Symbol, Any}(:d => "e")),
        Dict{String, Any}("x" => Float32(1.0e20), "y" => Float16(0.5)),
        Dict(1 => 2, 3 => 4),
        Set([1]),
        (1, "two", 3.0),
        JSON.parse("""{"z":1,"a":[true,false,null],"m":{"k":-2.5e-12}}"""),
        big(2)^70,
        'q',
        UInt8(5),
    ]
    for v in values
        @test HimalayaUI.to_json(v) == JSON.json(v; sort_keys = false)
    end
end

@testset "to_json refuses values JSON cannot hold" begin
    for x in (NaN, Inf, -Inf, Float32(NaN))
        @test_throws ArgumentError HimalayaUI.to_json(Dict(:r2 => x))
    end
end

using Test
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

using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "samples routes" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t", data_dir="/t/d",
                                             analysis_dir="/t/a")
    s1 = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    s2 = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D2", name="UX2")

    with_test_server(db) do port, base
        # List
        r = HTTP.get("$base/api/experiments/$exp_id/samples")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].label == "D1"
        @test list[1].tags  == []

        # PATCH
        r = HTTP.patch("$base/api/samples/$s1";
            body = JSON3.write(Dict(:notes => "hello")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).notes == "hello"

        # POST tag
        r = HTTP.post("$base/api/samples/$s1/tags";
            body = JSON3.write(Dict(:key => "lipid", :value => "DOPC")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 201
        tag = JSON3.read(String(r.body))
        @test tag.key   == "lipid"
        @test tag.value == "DOPC"
        tag_id = tag.id

        # Samples list now shows the tag
        r    = HTTP.get("$base/api/experiments/$exp_id/samples")
        list = JSON3.read(String(r.body))
        @test length(list[1].tags) == 1
        @test list[1].tags[1].key == "lipid"

        # DELETE tag
        r = HTTP.delete("$base/api/samples/$s1/tags/$tag_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        r    = HTTP.get("$base/api/experiments/$exp_id/samples")
        list = JSON3.read(String(r.body))
        @test list[1].tags == []
    end
end

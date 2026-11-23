using Test, HTTP, JSON3, SQLite, DBInterface

@testset "combine_traces" begin
    q = [0.01, 0.02, 0.03]
    a = (q, [10.0, 20.0, 30.0], [3.0, 4.0, 5.0])
    b = (q, [30.0, 40.0, 50.0], [4.0, 3.0, 12.0])
    cq, cI, cσ = HimalayaUI.combine_traces([a, b])
    @test cq == q
    @test cI ≈ [20.0, 30.0, 40.0]
    @test cσ ≈ [2.5, 2.5, 6.5]
    @test_throws ArgumentError HimalayaUI.combine_traces([a, ([0.01, 0.02], [1.0, 1.0], [1.0, 1.0])])
    @test_throws ArgumentError HimalayaUI.combine_traces([a, ([0.01, 0.02, 0.04], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0])])
    @test_throws ArgumentError HimalayaUI.combine_traces([])
end

function _combine_fixture()
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis")
    mkpath(analysis_dir)
    data = joinpath(@__DIR__, "..", "..", "..", "test", "data")
    cp(joinpath(data, "cubic_tot.dat"), joinpath(analysis_dir, "a_tot.dat"))
    cp(joinpath(data, "form-factor_tot.dat"), joinpath(analysis_dir, "b_tot.dat"))
    cp(joinpath(data, "cubic_tot.dat"), joinpath(analysis_dir, "c_tot.dat"))
    db = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path = tmp,
        data_dir = joinpath(tmp, "data"), analysis_dir = analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "S")
    ids = [HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = f,
               exposure_time = t, transmitted_counts = 100.0 * t, counting_time = t)
           for (f, t) in (("a_tot", 5.0), ("b_tot", 10.0), ("c_tot", 15.0))]
    for id in ids
        HimalayaUI.analyze_exposure!(db, id)
    end
    (; db, tmp, analysis_dir, exp_id, s_id, ids)
end

_sources(db, cid) = HimalayaUI.combined_source_ids(db, cid)
_kind_count(db, s_id) = Int(first(HimalayaUI.sqlrows(DBInterface.execute(db,
    "SELECT COUNT(*) AS n FROM exposures WHERE sample_id = ? AND kind = 'combined'", [s_id]))).n)

@testset "refresh_combined! builds and analyses one combined image per sample" begin
    fx = _combine_fixture()
    db = fx.db
    @test HimalayaUI.combined_exposure_id(db, fx.s_id) === nothing
    cid = HimalayaUI.refresh_combined!(db, fx.s_id)
    @test cid isa Int
    @test _sources(db, cid) == fx.ids
    row = first(HimalayaUI.sqlrows(DBInterface.execute(db,
        "SELECT kind, status, filename, exposure_time, trace_hash, analysis_version, analysis_error FROM exposures WHERE id = ?", [cid])))
    @test row.kind == "combined"
    @test ismissing(row.filename)
    @test row.exposure_time ≈ 30.0
    @test row.trace_hash == HimalayaUI.combined_trace_hash(db, cid)
    @test row.analysis_version == HimalayaUI.ANALYSIS_VERSION
    @test ismissing(row.analysis_error)
    @test HimalayaUI.count_auto_peaks(db, cid) > 0

    @test HimalayaUI.refresh_combined!(db, fx.s_id) == cid
    @test _kind_count(db, fx.s_id) == 1
    @test !HimalayaUI.combined_due(db, cid)

    a = HimalayaUI.load_dat(joinpath(fx.analysis_dir, "a_tot.dat"))
    b = HimalayaUI.load_dat(joinpath(fx.analysis_dir, "b_tot.dat"))
    _, I, σ = HimalayaUI.load_combined_trace(db, cid)
    @test I ≈ (2 .* a[2] .+ b[2]) ./ 3
    @test σ ≈ sqrt.(2 .* a[3] .^ 2 .+ b[3] .^ 2) ./ 3
end

@testset "a sample with one image gets no combined image" begin
    fx = _combine_fixture()
    db = fx.db
    DBInterface.execute(db, "UPDATE exposures SET status = 'rejected' WHERE id IN (?, ?)", fx.ids[1:2])
    @test HimalayaUI.refresh_combined!(db, fx.s_id) === nothing
    @test _kind_count(db, fx.s_id) == 0
end

@testset "images without a trace yet are left out of the combination" begin
    fx = _combine_fixture()
    db = fx.db
    DBInterface.execute(db, "UPDATE exposures SET trace_hash = NULL WHERE id = ?", [fx.ids[3]])
    cid = HimalayaUI.refresh_combined!(db, fx.s_id)
    @test _sources(db, cid) == fx.ids[1:2]
end

@testset "drop, restore and include routes keep the combined image in step" begin
    fx = _combine_fixture()
    db = fx.db
    cid = HimalayaUI.refresh_combined!(db, fx.s_id)
    first_hash = HimalayaUI.read_trace_hash(db, cid)
    json_hdr = ["Content-Type" => "application/json"]
    rep() = only(filter(s -> s.id == fx.s_id, JSON3.read(String(
        HimalayaUI.sqlrows === nothing ? "" : ""))))
    with_inproc_routes(db) do call
        representative() = only(filter(s -> s.id == fx.s_id,
            JSON3.read(String(call("GET", "/api/samples?experiment_id=$(fx.exp_id)").body)))).representative_exposure_id
        @test representative() == cid

        listing = JSON3.read(String(call("GET", "/api/samples/$(fx.s_id)/exposures").body))
        combined = only(filter(e -> e.id == cid, listing))
        @test combined.kind == "combined"
        @test sort([Int(s.source_exposure_id) for s in combined.sources]) == fx.ids
        @test all(e -> e.included === true, listing)
        @test !any(e -> haskey(e, :combine_excluded), listing)

        tr = JSON3.read(String(call("GET", "/api/exposures/$cid/trace").body))
        @test sort(collect(tr.sources)) == fx.ids
        @test tr.monitors.transmitted_counts ≈ 1000.0
        @test tr.monitors.exposure_time ≈ 1.0 && tr.monitors.counting_time ≈ 1.0

        r = call("PATCH", "/api/exposures/$(fx.ids[2])/status"; body = JSON3.write(Dict(:status => "rejected")), headers = json_hdr)
        @test r.status == 200
        @test _sources(db, cid) == [fx.ids[1], fx.ids[3]]
        @test HimalayaUI.read_trace_hash(db, cid) != first_hash
        @test HimalayaUI.read_trace_hash(db, cid) == HimalayaUI.combined_trace_hash(db, cid)

        r = call("PATCH", "/api/exposures/$(fx.ids[3])/included"; body = JSON3.write(Dict(:included => false)), headers = json_hdr)
        @test r.status == 200
        @test _sources(db, cid) == [fx.ids[1]]
        @test representative() == fx.ids[3]
        @test call("GET", "/api/exposures/$cid/trace").status == 404
        listing = JSON3.read(String(call("GET", "/api/samples/$(fx.s_id)/exposures").body))
        @test only(filter(e -> e.id == fx.ids[3], listing)).included === false

        call("PATCH", "/api/exposures/$(fx.ids[2])/status"; body = JSON3.write(Dict(:status => "accepted")), headers = json_hdr)
        call("PATCH", "/api/exposures/$(fx.ids[3])/included"; body = JSON3.write(Dict(:included => true)), headers = json_hdr)
        @test _sources(db, cid) == fx.ids
        @test HimalayaUI.read_trace_hash(db, cid) == first_hash
        @test representative() == cid
        @test _kind_count(db, fx.s_id) == 1

        @test call("PATCH", "/api/exposures/$cid/status"; body = JSON3.write(Dict(:status => "rejected")), headers = json_hdr).status == 422
        @test call("PATCH", "/api/exposures/$cid/included"; body = JSON3.write(Dict(:included => false)), headers = json_hdr).status == 422
        @test call("PATCH", "/api/exposures/$(fx.ids[1])/included"; body = JSON3.write(Dict(:included => "no")), headers = json_hdr).status == 422
        @test call("PATCH", "/api/exposures/99999/included"; body = JSON3.write(Dict(:included => true)), headers = json_hdr).status == 404
    end
end

@testset "a combined image whose sources disagree on q records the error once" begin
    fx = _combine_fixture()
    db = fx.db
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(fx.analysis_dir, "c_tot.dat"); force = true)
    HimalayaUI.analyze_exposure!(db, fx.ids[3])
    cid = HimalayaUI.refresh_combined!(db, fx.s_id)
    err = first(HimalayaUI.sqlrows(DBInterface.execute(db,
        "SELECT analysis_error FROM exposures WHERE id = ?", [cid]))).analysis_error
    @test occursin("q grids", err)
    @test !HimalayaUI.combined_due(db, cid)
end

@testset "migration adds combine_excluded to an old exposures table" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.migrate_schema!(db)
    DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN combine_excluded")
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_EXPOSURES_COMBINE_EXCLUDED])
    HimalayaUI.migrate_schema!(db)
    @test "combine_excluded" in HimalayaUI.cols_of(db, "exposures")
    HimalayaUI.migrate_schema!(db)
    @test "combine_excluded" in HimalayaUI.cols_of(db, "exposures")
end

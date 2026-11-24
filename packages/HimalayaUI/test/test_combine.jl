using Test, HTTP, JSON3, SQLite, DBInterface, Tables

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

_select!(db, ids) = DBInterface.execute(db,
    "UPDATE exposures SET deselected = CASE WHEN id IN ($(join(ids, ","))) THEN 0 ELSE 1 END WHERE kind = 'file'")

@testset "refresh_selection! averages every selected image, one curve per set" begin
    fx = _combine_fixture()
    db = fx.db
    @test HimalayaUI.selected_image_ids(db, fx.s_id) == fx.ids
    cid = HimalayaUI.refresh_selection!(db, fx.s_id)
    @test cid isa Int
    @test HimalayaUI.analyzed_exposure_id(db, fx.s_id) == cid
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
    run = first(HimalayaUI.sqlrows(DBInterface.execute(db,
        "SELECT payload FROM user_actions WHERE action = 'analyze_run' AND entity_id = ? ORDER BY id DESC LIMIT 1", [cid])))
    @test JSON3.read(run.payload).average_of_sample == fx.s_id
    single = first(HimalayaUI.sqlrows(DBInterface.execute(db,
        "SELECT payload FROM user_actions WHERE action = 'analyze_run' AND entity_id = ? ORDER BY id DESC LIMIT 1", [fx.ids[1]])))
    @test !haskey(JSON3.read(single.payload), :average_of_sample)
    @test HimalayaUI.refresh_selection!(db, fx.s_id) == cid
    @test _kind_count(db, fx.s_id) == 1
    @test !HimalayaUI.combined_due(db, cid)

    a = HimalayaUI.load_dat(joinpath(fx.analysis_dir, "a_tot.dat"))
    b = HimalayaUI.load_dat(joinpath(fx.analysis_dir, "b_tot.dat"))
    _, I, σ = HimalayaUI.load_combined_trace(db, cid)
    @test I ≈ (2 .* a[2] .+ b[2]) ./ 3
    @test σ ≈ sqrt.(2 .* a[3] .^ 2 .+ b[3] .^ 2) ./ 3

    _select!(db, fx.ids[1:2])
    pair = HimalayaUI.refresh_selection!(db, fx.s_id)
    @test pair != cid
    @test _sources(db, pair) == fx.ids[1:2]
    @test _sources(db, cid) == fx.ids
    @test _kind_count(db, fx.s_id) == 2

    _select!(db, fx.ids)
    @test HimalayaUI.refresh_selection!(db, fx.s_id) == cid
    @test _kind_count(db, fx.s_id) == 2
end

@testset "one selected image is analyzed as itself, none selected analyzes nothing" begin
    fx = _combine_fixture()
    db = fx.db
    _select!(db, [fx.ids[2]])
    @test HimalayaUI.refresh_selection!(db, fx.s_id) == fx.ids[2]
    @test HimalayaUI.analyzed_exposure_id(db, fx.s_id) == fx.ids[2]
    @test _kind_count(db, fx.s_id) == 0
    _select!(db, Int[0])
    @test HimalayaUI.refresh_selection!(db, fx.s_id) === nothing
    @test HimalayaUI.analyzed_exposure_id(db, fx.s_id) === nothing
end

@testset "a set with an image not yet analyzed gets no average yet" begin
    fx = _combine_fixture()
    db = fx.db
    DBInterface.execute(db, "UPDATE exposures SET trace_hash = NULL WHERE id = ?", [fx.ids[3]])
    @test HimalayaUI.refresh_selection!(db, fx.s_id) === nothing
    @test _kind_count(db, fx.s_id) == 0
end

@testset "selection, drop and restore routes" begin
    fx = _combine_fixture()
    db = fx.db
    cid = HimalayaUI.refresh_selection!(db, fx.s_id)
    json_hdr = ["Content-Type" => "application/json"]
    put(ids) = call -> call("PUT", "/api/samples/$(fx.s_id)/selection";
        body = JSON3.write(Dict(:exposure_ids => ids)), headers = json_hdr)
    status(call, id, st) = call("PATCH", "/api/exposures/$id/status";
        body = JSON3.write(Dict(:status => st)), headers = json_hdr)
    with_inproc_routes(db) do call
        analyzed() = only(filter(s -> s.id == fx.s_id,
            JSON3.read(String(call("GET", "/api/samples?experiment_id=$(fx.exp_id)").body)))).analyzed_exposure_id
        listing() = JSON3.read(String(call("GET", "/api/samples/$(fx.s_id)/exposures").body))
        gap() = only(filter(s -> s.id == fx.s_id,
            JSON3.read(String(call("GET", "/api/samples?experiment_id=$(fx.exp_id)").body)))).analysis_gap
        @test analyzed() == cid
        @test gap() === nothing

        combined = only(filter(e -> e.id == cid, listing()))
        @test combined.kind == "combined"
        @test !haskey(combined, :selected)
        @test sort([Int(s.source_exposure_id) for s in combined.sources]) == fx.ids
        @test all(e -> e.selected === true, filter(e -> e.kind == "file", listing()))
        @test !any(e -> haskey(e, :deselected), listing())

        tr = JSON3.read(String(call("GET", "/api/exposures/$cid/trace").body))
        @test sort(collect(tr.sources)) == fx.ids
        @test tr.monitors.transmitted_counts ≈ 1000.0
        @test tr.monitors.i2_collapsed == false
        @test haskey(tr.monitors, :mode)
        DBInterface.execute(db, "UPDATE exposures SET incident_monitor = 1.5 WHERE kind = 'file'")
        original = only(Tables.rowtable(DBInterface.execute(db, "SELECT transmitted_counts FROM exposures WHERE id = ?", [fx.ids[2]]))).transmitted_counts
        DBInterface.execute(db, "UPDATE exposures SET transmitted_counts = 1.0 WHERE id = ?", [fx.ids[2]])
        @test JSON3.read(String(call("GET", "/api/exposures/$cid/trace").body)).monitors.i2_collapsed == true
        DBInterface.execute(db, "UPDATE exposures SET transmitted_counts = ? WHERE id = ?", [original, fx.ids[2]])

        r = put([fx.ids[3], fx.ids[1]])(call)
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test collect(body.exposure_ids) == [fx.ids[1], fx.ids[3]]
        pair = body.analyzed_exposure_id
        @test _sources(db, pair) == [fx.ids[1], fx.ids[3]]
        @test analyzed() == pair
        @test only(filter(e -> e.id == fx.ids[2], listing())).selected === false

        r = put([fx.ids[2]])(call)
        @test JSON3.read(String(r.body)).analyzed_exposure_id == fx.ids[2]
        @test analyzed() == fx.ids[2]

        r = put(Int[])(call)
        @test r.status == 200
        @test JSON3.read(String(r.body)).analyzed_exposure_id === nothing
        @test analyzed() === nothing
        @test gap() == "none_selected"

        put(fx.ids)(call)
        @test analyzed() == cid
        @test status(call, fx.ids[2], "rejected").status == 200
        @test HimalayaUI.selected_image_ids(db, fx.s_id) == [fx.ids[1], fx.ids[3]]
        @test analyzed() == HimalayaUI.combined_for_set(db, fx.s_id, [fx.ids[1], fx.ids[3]])

        @test status(call, fx.ids[2], "accepted").status == 200
        @test only(filter(e -> e.id == fx.ids[2], listing())).selected === false
        @test status(call, fx.ids[2], "rejected").status == 200
        undo = call("PATCH", "/api/exposures/$(fx.ids[2])/status"; body = JSON3.write(Dict(:status => nothing, :selected => true)), headers = json_hdr)
        @test undo.status == 200
        @test HimalayaUI.selected_image_ids(db, fx.s_id) == fx.ids
        @test call("PATCH", "/api/exposures/$(fx.ids[2])/status"; body = JSON3.write(Dict(:status => nothing, :selected => "yes")), headers = json_hdr).status == 422
        @test put([fx.ids[1], fx.ids[3]])(call).status == 200
        @test analyzed() == HimalayaUI.combined_for_set(db, fx.s_id, [fx.ids[1], fx.ids[3]])

        @test put([fx.ids[2]])(call).status == 200
        @test status(call, fx.ids[1], "rejected").status == 200
        @test status(call, fx.ids[3], "rejected").status == 200
        @test status(call, fx.ids[2], "rejected").status == 200
        @test isempty(HimalayaUI.selected_image_ids(db, fx.s_id))
        @test analyzed() === nothing
        @test gap() == "all_dropped"
        @test put(Int[])(call).status == 200
        @test status(call, fx.ids[1], "accepted").status == 200
        @test HimalayaUI.selected_image_ids(db, fx.s_id) == [fx.ids[1]]

        @test status(call, cid, "rejected").status == 422
        @test put([fx.ids[3]])(call).status == 422
        @test put([cid])(call).status == 422
        @test call("PUT", "/api/samples/$(fx.s_id)/selection"; body = JSON3.write(Dict(:exposure_ids => "all")), headers = json_hdr).status == 422
        @test call("PUT", "/api/samples/99999/selection"; body = JSON3.write(Dict(:exposure_ids => Int[])), headers = json_hdr).status == 404
    end
end

@testset "the selection route refuses a sample of an autosampler setup" begin
    fx = _combine_fixture()
    db = fx.db
    auto = HimalayaUI.create_setup!(db; experiment_id = fx.exp_id, folder = "solution",
        root_dir = fx.tmp, data_dir = fx.tmp, analysis_dir = fx.analysis_dir, mode = "autosampler")
    water = HimalayaUI.create_sample!(db; experiment_id = fx.exp_id, name = "water 01A (B0101)", setup_id = auto)
    image = HimalayaUI.create_exposure!(db; experiment_id = fx.exp_id, sample_id = water, setup_id = auto, filename = "a_tot")
    with_inproc_routes(db) do call
        resp = call("PUT", "/api/samples/$water/selection";
            body = JSON3.write(Dict(:exposure_ids => Int[])), headers = ["Content-Type" => "application/json"])
        @test resp.status == 422
        @test JSON3.read(String(resp.body)).error == "An autosampler sample has one image, so there is nothing to select."
        @test HimalayaUI.selected_image_ids(db, water) == [image]
    end
end

@testset "a rescan keeps a custom selection and grows an all-selected one" begin
    fx = _combine_fixture()
    db = fx.db
    other = HimalayaUI.create_sample!(db; experiment_id = fx.exp_id, name = "T")
    o1 = HimalayaUI.create_exposure!(db; experiment_id = fx.exp_id, sample_id = other, filename = "o1")
    _select!(db, vcat(fx.ids[1:2], [o1]))
    n1 = HimalayaUI.create_exposure!(db; experiment_id = fx.exp_id, sample_id = fx.s_id, filename = "n1")
    n2 = HimalayaUI.create_exposure!(db; experiment_id = fx.exp_id, sample_id = other, filename = "n2")
    HimalayaUI._keep_custom_selections!(db, [n1, n2])
    @test HimalayaUI.selected_image_ids(db, fx.s_id) == fx.ids[1:2]
    @test HimalayaUI.selected_image_ids(db, other) == [o1, n2]
end

@testset "an average whose images disagree on q records the error once" begin
    fx = _combine_fixture()
    db = fx.db
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(fx.analysis_dir, "c_tot.dat"); force = true)
    HimalayaUI.analyze_exposure!(db, fx.ids[3])
    cid = HimalayaUI.refresh_selection!(db, fx.s_id)
    err = first(HimalayaUI.sqlrows(DBInterface.execute(db,
        "SELECT analysis_error FROM exposures WHERE id = ?", [cid]))).analysis_error
    @test occursin("q grids", err)
    @test !HimalayaUI.combined_due(db, cid)
end

@testset "migration adds deselected and keeps reviewed samples on the image they were worked on" begin
    fx = _combine_fixture()
    db = fx.db
    plain = HimalayaUI.create_sample!(db; experiment_id = fx.exp_id, name = "P")
    p = [HimalayaUI.create_exposure!(db; experiment_id = fx.exp_id, sample_id = plain, filename = "p$i") for i in 1:2]
    ix = Int(DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO indices (exposure_id, phase, basis, score) VALUES (?, 'Pn3m', 0.1, 0.9)", [fx.ids[2]])))
    DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [fx.ids[2], ix])
    DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN deselected")
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_EXPOSURES_DESELECTED])
    HimalayaUI.migrate_schema!(db)
    @test "deselected" in HimalayaUI.cols_of(db, "exposures")
    @test HimalayaUI.selected_image_ids(db, fx.s_id) == [fx.ids[2]]
    @test HimalayaUI.selected_image_ids(db, plain) == p
    HimalayaUI.migrate_schema!(db)
    @test HimalayaUI.selected_image_ids(db, fx.s_id) == [fx.ids[2]]
end

@testset "migration keeps the image holding the call selected over a newer curated image" begin
    fx = _combine_fixture()
    db = fx.db
    ix = Int(DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO indices (exposure_id, phase, basis, score) VALUES (?, 'Pn3m', 0.1, 0.9)", [fx.ids[1]])))
    DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [fx.ids[1], ix])
    DBInterface.execute(db,
        "INSERT INTO user_actions (entity_type, entity_id, action, payload) VALUES ('exposure', ?, 'peak_added', '{}')", [fx.ids[3]])
    DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN deselected")
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_EXPOSURES_DESELECTED])
    HimalayaUI.migrate_schema!(db)
    @test HimalayaUI.selected_image_ids(db, fx.s_id) == [fx.ids[1]]
end

function _rerun_selection_migration!(db)
    DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN deselected")
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_EXPOSURES_DESELECTED])
    HimalayaUI.migrate_schema!(db)
end

_peak_edit!(db, id) = DBInterface.execute(db,
    "INSERT INTO user_actions (entity_type, entity_id, action, payload) VALUES ('exposure', ?, 'peak_added', '{}')", [id])

@testset "migration falls back to the newest image with peak edits when no kept image holds the call" begin
    fx = _combine_fixture()
    _peak_edit!(fx.db, fx.ids[1]); _peak_edit!(fx.db, fx.ids[2])
    _rerun_selection_migration!(fx.db)
    @test HimalayaUI.selected_image_ids(fx.db, fx.s_id) == [fx.ids[2]]

    fx = _combine_fixture()
    ix = Int(DBInterface.lastrowid(DBInterface.execute(fx.db,
        "INSERT INTO indices (exposure_id, phase, basis, score) VALUES (?, 'Pn3m', 0.1, 0.9)", [fx.ids[3]])))
    DBInterface.execute(fx.db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [fx.ids[3], ix])
    DBInterface.execute(fx.db, "UPDATE exposures SET status = 'rejected' WHERE id = ?", [fx.ids[3]])
    _peak_edit!(fx.db, fx.ids[1])
    _rerun_selection_migration!(fx.db)
    @test HimalayaUI.selected_image_ids(fx.db, fx.s_id) == [fx.ids[1]]
end

@testset "analyze skips averages of an earlier selection" begin
    fx = _combine_fixture()
    db = fx.db
    stale = HimalayaUI.refresh_selection!(db, fx.s_id)
    _select!(db, fx.ids[1:2])
    current = HimalayaUI.refresh_selection!(db, fx.s_id)
    @test current != stale
    DBInterface.execute(db, "UPDATE exposures SET status = 'rejected' WHERE id = ?", [fx.ids[3]])
    DBInterface.execute(db, "UPDATE exposures SET trace_hash = NULL, analysis_inputs_hash = NULL WHERE experiment_id = ?", [fx.exp_id])
    out = mktemp() do path, io
        skipped = redirect_stdout(io) do
            HimalayaUI._analyze_experiment!(db, fx.exp_id)
        end
        @test skipped == 0
        close(io)
        read(path, String)
    end
    @test occursin("average of 3 images (not the current selection)", out)
    @test occursin("Analyzing S / average of 2 images ... done", out)
    hashed(id) = !ismissing(first(HimalayaUI.sqlrows(DBInterface.execute(db,
        "SELECT analysis_inputs_hash AS h FROM exposures WHERE id = ?", [id]))).h)
    @test hashed(current)
    @test !hashed(stale)
end

@testset "backfill_averages! builds each selection's missing average once" begin
    fx = _combine_fixture()
    db = fx.db
    other = HimalayaUI.create_sample!(db; experiment_id = fx.exp_id, name = "T")
    HimalayaUI.create_exposure!(db; experiment_id = fx.exp_id, sample_id = other, filename = "t_tot")
    @test HimalayaUI.analyzed_exposure_id(db, fx.s_id) === nothing
    @test HimalayaUI.backfill_averages!(db) == 1
    cid = HimalayaUI.analyzed_exposure_id(db, fx.s_id)
    @test cid isa Int
    runs() = Int(first(HimalayaUI.sqlrows(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM user_actions WHERE action = 'analyze_run'"))).n)
    before = runs()
    HimalayaUI.backfill_averages!(db)
    @test runs() == before
    @test _kind_count(db, fx.s_id) == 1
end

@testset "export lists the current average, not averages of an earlier selection" begin
    fx = _combine_fixture()
    db = fx.db
    full = HimalayaUI.refresh_selection!(db, fx.s_id)
    _select!(db, fx.ids[1:2])
    pair = HimalayaUI.refresh_selection!(db, fx.s_id)
    @test pair != full
    ex = only(HimalayaUI.build_export(db, fx.exp_id))[:exposures]
    averages = filter(e -> e[:kind] == "combined", ex)
    @test length(averages) == 1
    @test only(averages)[:analyzed]
    @test length(only(averages)[:sources]) == 2
end

@testset "a bridged sample averages only its own images" begin
    fx = _combine_fixture()
    db = fx.db
    data = joinpath(@__DIR__, "..", "..", "..", "test", "data")
    mate = HimalayaUI.create_sample!(db; experiment_id = fx.exp_id, name = "S in another setup")
    mate_ids = map(("m1_tot", "m2_tot")) do f
        cp(joinpath(data, "cubic_tot.dat"), joinpath(fx.analysis_dir, f * ".dat"))
        id = HimalayaUI.create_exposure!(db; experiment_id = fx.exp_id, sample_id = mate, filename = f,
            exposure_time = 5.0, transmitted_counts = 500.0, counting_time = 5.0)
        HimalayaUI.analyze_exposure!(db, id)
        id
    end
    DBInterface.execute(db, "INSERT INTO sample_bridges (sample_id, group_id) VALUES (?, 1), (?, 1)", [fx.s_id, mate])
    @test HimalayaUI.bridge_mates(db, fx.s_id) == [mate]
    @test _sources(db, HimalayaUI.refresh_selection!(db, fx.s_id)) == fx.ids
    @test _sources(db, HimalayaUI.refresh_selection!(db, mate)) == sort(collect(mate_ids))
end

using Test, Dates, SQLite, DBInterface, Tables, JSON3, HimalayaUI
isdefined(@__MODULE__, :write_prp) || include("test_ingestion_core.jl")

_as_count(db, sql, args...) = Int(only(Tables.rowtable(DBInterface.execute(db, sql, collect(args)))).c)

_as_stem(label, well, kind, scan) = "$(label)_$(well)_$(kind)$(scan)"

function _as_write_dat(path)
    q = collect(0.01:0.01:0.6)
    open(path, "w") do io
        for x in q
            println(io, x, " ", 1000 / (1 + 100 * x^2) + 10, " ", 3.0)
        end
    end
end

_as_add_acquisition!(root, (label, well, kind, scan); kw...) =
    _as_add_stem!(root, _as_stem(label, well, kind, scan); kw...)

function _as_add_stem!(root, stem; frames = 2, minute = 0)
    data = mkpath(joinpath(root, "data"))
    tot = mkpath(joinpath(root, "analysis", "automatic_analysis", "tot_files"))
    t = DateTime(2026, 4, 26, 12, 0, 0) + Minute(minute)
    for f in 1:frames
        frame = lpad(f, 3, '0')
        write_prp(joinpath(data, "$(stem)_0_$(frame).prp");
            timestamp = Dates.format(t + Second(f), "dd u yyyy HH:MM:SS"), horizontal_position_mm = 58.9)
        write(joinpath(data, "$(stem)_0_$(frame).tif"), "fake tif")
    end
    _as_write_dat(joinpath(tot, "$(stem)_tot.dat"))
    stem
end

function _as_add_at!(root, label, scan)
    write(joinpath(mkpath(joinpath(root, "data")), "AT_SolSAXS_$(label)_$(scan).txt"), "load 30 5 1 3\n")
end

function _as_write_folder!(root, acquisitions; at = Tuple{String,Int}[])
    for (i, a) in enumerate(acquisitions)
        _as_add_acquisition!(root, a; minute = i)
    end
    for (label, scan) in at
        _as_add_at!(root, label, scan)
    end
    root
end

const _AS_RUN_A = [("wat", "01A", "B", 100), ("wat", "01B", "S", 101), ("wat", "01C", "S", 102)]
const _AS_RUN_B = [("pep", "01A", "B", 103), ("pep", "01B", "S", 104)]

@testset "setup mode" begin
    @testset "a setup stores its mode, and an undecided one stores NULL" begin
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "x", path = "/x")
        a = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = "/x/a", analysis_dir = "/x/a")
        b = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = "/x/b", analysis_dir = "/x/b",
                                     position = 1, mode = "autosampler")
        @test ismissing(HimalayaUI.setup_row(db, a).mode)
        @test HimalayaUI.setup_row(db, b).mode == "autosampler"
        @test_throws ArgumentError HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = "/x/c",
            analysis_dir = "/x/c", position = 2, mode = "powder")
    end

    @testset "detect_setup_mode goes by the autosampler's run files" begin
        auto = _as_write_folder!(mktempdir(), _AS_RUN_A; at = [("wat", 100)])
        @test HimalayaUI.detect_setup_mode(joinpath(auto, "data")) == "autosampler"

        well_named = _as_write_folder!(mktempdir(), _AS_RUN_A)
        @test HimalayaUI.detect_setup_mode(joinpath(well_named, "data")) == "capillary"

        cap = mktempdir()
        data = mkpath(joinpath(cap, "data"))
        for i in 1:6
            write(joinpath(data, "DOPE_pH7_37C_S$(lpad(i, 4, '0'))_0_001.tif"), "fake tif")
        end
        @test HimalayaUI.detect_setup_mode(data) == "capillary"

        @test HimalayaUI.detect_setup_mode(mktempdir()) === nothing
        @test HimalayaUI.detect_setup_mode(joinpath(cap, "missing")) === nothing
        onlyat = mktempdir()
        _as_add_at!(onlyat, "wat", 100)
        @test HimalayaUI.detect_setup_mode(joinpath(onlyat, "data")) == "autosampler"
    end

    @testset "the layout of an autosampler folder points at the averaged traces" begin
        root = _as_write_folder!(mktempdir(), _AS_RUN_A; at = [("wat", 100)])
        r = HimalayaUI.resolve_experiment_layout(root)
        @test r.mode == "autosampler"
        @test r.analysis_dir == joinpath(root, "analysis", "automatic_analysis", "tot_files")
        @test r.image_pattern == "{name}_0_001.tif"
        @test r.metadata_pattern == "{name}_0_001.prp"
        @test r.integration_pattern == "{name}_tot.dat"

        cap = mktempdir(); mkpath(joinpath(cap, "data"))
        write(joinpath(cap, "data", "MO_S0001_0_001.tif"), "fake tif")
        @test HimalayaUI.resolve_experiment_layout(cap).mode == "capillary"
        @test HimalayaUI.resolve_experiment_layout(mktempdir()).mode === nothing
    end

    @testset "migration marks setups that already have exposures as capillary" begin
        path = joinpath(mktempdir(), "m.db")
        db = HimalayaUI.open_db(path)
        eid = HimalayaUI.create_experiment!(db; name = "x", path = "/x")
        used = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = "/x/a", analysis_dir = "/x/a")
        empty_setup = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = "/x/b", analysis_dir = "/x/b",
                                               position = 1, included = false, image_pattern = "{name}_0_009.tif")
        sid = HimalayaUI.create_sample!(db; experiment_id = eid, name = "S", setup_id = used)
        HimalayaUI.create_exposure!(db; experiment_id = eid, setup_id = used, sample_id = sid, filename = "A_S1_0_001")
        DBInterface.execute(db, "ALTER TABLE setups DROP COLUMN mode")
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_SETUP_MODE])
        SQLite.close(db)

        db2 = HimalayaUI.open_db(path)
        @test HimalayaUI.setup_row(db2, used).mode == "capillary"
        @test ismissing(HimalayaUI.setup_row(db2, empty_setup).mode)
        DBInterface.execute(db2, "UPDATE setups SET mode = 'autosampler' WHERE id = ?", [used])
        SQLite.close(db2)

        db3 = HimalayaUI.open_db(path)
        @test HimalayaUI.setup_row(db3, used).mode == "autosampler"
        SQLite.close(db3)
    end
end

function _as_grouping(root; at = Tuple{String,Int}[])
    data = joinpath(root, "data")
    layout = HimalayaUI.autosampler_layout(root)
    metas = HimalayaUI.scan_directory(data, layout.dir;
        tif_pattern = layout.image_pattern, prp_pattern = layout.metadata_pattern,
        dat_pattern = layout.integration_pattern)
    HimalayaUI.group_autosampler(metas, data)
end

_as_names(res) = [s.name for l in res.loads for s in l.samples]

function _as_live_row_db(root)
    db = fresh_db()
    eid = HimalayaUI.create_experiment!(db; name = "as", path = root)
    sid = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = joinpath(root, "data"),
        analysis_dir = joinpath(root, "analysis"), root_dir = root,
        image_pattern = "{name}_0_009.tif", metadata_pattern = "{name}_0_009.prp",
        integration_pattern = "{name}_tot.dat")
    (db, eid, sid)
end

@testset "autosampler grouping" begin
    @testset "one sample per acquisition, loads split by label and by run files" begin
        runs = vcat(_AS_RUN_A, _AS_RUN_B, [("pep", "01A", "B", 105), ("pep", "01B", "S", 106)])
        root = _as_write_folder!(mktempdir(), runs; at = [("wat", 100), ("pep", 103), ("pep", 105)])
        write(joinpath(root, "data", "agbe_S0099_0_001.tif"), "fake tif")
        res = _as_grouping(root)
        @test [length(l.samples) for l in res.loads] == [1, 3, 2, 2]
        @test [l.load_index for l in res.loads] == [99, 100, 103, 105]
        @test [s.slot_index for l in res.loads for s in l.samples] == collect(99:106)
        @test _as_names(res) == ["agbe_S0099", "wat 01A (B100)", "wat 01B (S101)", "wat 01C (S102)",
            "pep 01A (B103)", "pep 01B (S104)", "pep 01A (B105)", "pep 01B (S106)"]
        @test all(s -> length(s.exposures) == 1, (s for l in res.loads for s in l.samples))
        first_exposure = res.loads[2].samples[1].exposures[1]
        @test first_exposure.stem == "wat_01A_B100"
        @test first_exposure.scan_id == 100
        @test res.loads[2].samples[2].exposures[1].scan_id == 101
        @test all(s -> s.grouping_source == "auto_acquisition", (s for l in res.loads for s in l.samples))
    end

    @testset "without run files a change of label starts a load" begin
        root = _as_write_folder!(mktempdir(), vcat(_AS_RUN_A, _AS_RUN_B))
        res = _as_grouping(root)
        @test [length(l.samples) for l in res.loads] == [3, 2]
        @test [l.load_index for l in res.loads] == [100, 103]
    end

    @testset "a same-label run needs its run file to be a separate load" begin
        root = _as_write_folder!(mktempdir(), vcat(_AS_RUN_A, [("wat", "01A", "B", 103)]))
        @test length(_as_grouping(root).loads) == 1
    end
end

@testset "autosampler scan" begin
    runs = vcat(_AS_RUN_A, _AS_RUN_B)
    function scanned()
        root = _as_write_folder!(mktempdir(), runs; at = [("wat", 100), ("pep", 103)])
        db, eid, sid = _as_live_row_db(root)
        r = HimalayaUI.scan_and_group!(db, eid)
        (; root, db, eid, sid, r)
    end

    @testset "an undecided setup with the wrong patterns is put right on its first scan" begin
        (; root, db, eid, sid, r) = scanned()
        st = HimalayaUI.setup_row(db, sid)
        @test st.mode == "autosampler"
        @test st.image_pattern == "{name}_0_001.tif"
        @test st.metadata_pattern == "{name}_0_001.prp"
        @test st.integration_pattern == "{name}_tot.dat"
        @test st.analysis_dir == joinpath(root, "analysis", "automatic_analysis", "tot_files")
        @test r.status == :ok
        @test r.added_exposures == 5 && r.added_samples == 5 && r.added_loads == 2
        @test _as_count(db, "SELECT COUNT(*) AS c FROM samples WHERE setup_id = ?", sid) == 5
        @test _as_count(db, "SELECT COUNT(DISTINCT sample_id) AS c FROM exposures WHERE setup_id = ?", sid) == 5
        filenames = [x.filename for x in Tables.rowtable(DBInterface.execute(db,
            "SELECT filename FROM exposures WHERE setup_id = ? ORDER BY id", [sid]))]
        @test filenames == ["wat_01A_B100", "wat_01B_S101", "wat_01C_S102", "pep_01A_B103", "pep_01B_S104"]
    end

    @testset "a capillary setup whose exposures carry well tokens stays capillary" begin
        (; root, db, eid, sid) = scanned()
        DBInterface.execute(db, "UPDATE setups SET mode = 'capillary' WHERE id = ?", [sid])
        analysis = joinpath(root, "analysis")
        DBInterface.execute(db, "UPDATE setups SET analysis_dir = ? WHERE id = ?", [analysis, sid])
        HimalayaUI.scan_and_group!(db, eid)
        st = HimalayaUI.setup_row(db, sid)
        @test st.mode == "capillary"
        @test st.analysis_dir == analysis
    end

    @testset "nothing is analyzed and nothing is reported as a problem" begin
        (; db, eid, sid) = scanned()
        @test _as_count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ? AND (trace_hash IS NOT NULL OR analysis_error IS NOT NULL OR analysis_version IS NOT NULL)", sid) == 0
        @test _as_count(db, "SELECT COUNT(*) AS c FROM auto_peaks") == 0
        @test isempty(HimalayaUI.ingest_problems(db, eid))
        @test isempty(HimalayaUI.pending_analysis_ids(db, eid))
        @test isempty(HimalayaUI.pending_analysis_ids(db, eid; retry_failed = true))
        @test !HimalayaUI.ingest_has_work(db, eid)
    end

    @testset "new images count as work, analysis never does" begin
        (; root, db, eid) = scanned()
        _as_add_acquisition!(root, ("pep", "01C", "S", 105))
        @test HimalayaUI.ingest_has_work(db, eid)
        HimalayaUI.scan_and_group!(db, eid)
        @test !HimalayaUI.ingest_has_work(db, eid)
        @test isempty(HimalayaUI.pending_analysis_ids(db, eid))
    end

    @testset "each sample still resolves an image and its trace" begin
        (; db, eid, sid) = scanned()
        samples = Tables.rowtable(DBInterface.execute(db, "SELECT id FROM samples WHERE setup_id = ? ORDER BY id", [sid]))
        shown = [HimalayaUI.analyzed_exposure_id(db, Int(s.id)) for s in samples]
        @test all(!isnothing, shown)
        @test length(unique(shown)) == 5
        with_inproc_routes(db) do call
            r = call("GET", "/api/exposures/$(first(shown))/trace")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test length(body.q) == length(body.I) > 10
        end
    end

    @testset "a rescan keeps the existing samples and loads where they were" begin
        root = _as_write_folder!(mktempdir(), _AS_RUN_A; at = [("wat", 100)])
        db, eid, sid = _as_live_row_db(root)
        HimalayaUI.scan_and_group!(db, eid)
        snapshot() = [(x.id, x.name, x.load_id, x.slot_index) for x in Tables.rowtable(DBInterface.execute(db,
            "SELECT id, name, load_id, slot_index FROM samples WHERE setup_id = ? ORDER BY id", [sid]))]
        before = snapshot()
        @test length(before) == 3

        for a in (("wat", "01D", "S", 103), ("pep", "01A", "B", 104), ("pep", "01B", "S", 105))
            _as_add_acquisition!(root, a; minute = a[4])
        end
        _as_add_at!(root, "pep", 104)
        r = HimalayaUI.scan_and_group!(db, eid)
        after = snapshot()
        @test after[1:3] == before
        @test [x[2] for x in after[4:end]] == ["wat 01D (S103)", "pep 01A (B104)", "pep 01B (S105)"]
        @test r.added_exposures == 3
        @test _as_count(db, "SELECT COUNT(*) AS c FROM loads WHERE setup_id = ?", sid) == 2
        @test HimalayaUI.scan_and_group!(db, eid).added_exposures == 0
    end
end

@testset "automatic analysis passes skip autosampler setups" begin
    root = _as_write_folder!(mktempdir(), _AS_RUN_A; at = [("wat", 100)])
    db, eid, sid = _as_live_row_db(root)
    HimalayaUI.scan_and_group!(db, eid)
    sample_ids = [Int(x.id) for x in Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM samples WHERE setup_id = ? ORDER BY id", [sid]))]
    untouched() = _as_count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ? AND (trace_hash IS NOT NULL OR analysis_version IS NOT NULL)", sid) == 0

    @testset "the command line analyze" begin
        redirect_stdout(devnull) do
            HimalayaUI._analyze_experiment!(db, eid)
        end
        @test untouched()
    end

    @testset "selecting several images makes no average" begin
        DBInterface.execute(db, "UPDATE exposures SET sample_id = ? WHERE setup_id = ? AND sample_id != ?",
            [sample_ids[1], sid, sample_ids[1]])
        @test HimalayaUI.refresh_selection!(db, sample_ids[1]) === nothing
        @test _as_count(db, "SELECT COUNT(*) AS c FROM exposures WHERE kind = 'combined'") == 0
        @test HimalayaUI.backfill_averages!(db) >= 0
        @test _as_count(db, "SELECT COUNT(*) AS c FROM exposures WHERE kind = 'combined'") == 0
        @test untouched()
    end

    @testset "the peak model pass leaves them alone" begin
        DBInterface.execute(db, """UPDATE exposures SET trace_hash = 'h', analysis_version = ?
                                    WHERE setup_id = ?""", [HimalayaUI.ANALYSIS_VERSION, sid])
        @test isempty(HimalayaUI.peak_model_pass_ids(db, 1))
        DBInterface.execute(db, "UPDATE setups SET mode = 'capillary' WHERE id = ?", [sid])
        @test length(HimalayaUI.peak_model_pass_ids(db, 1)) == 3
    end
end

function _as_wait_ingest(db, eid; timeout = 60)
    t0 = time()
    while time() - t0 < timeout
        st = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT ingest_status FROM experiments WHERE id = ?", [eid]))).ingest_status
        isequal(st, "scanning") || HimalayaUI.ingest_running(eid) || return st
        sleep(0.05)
    end
    error("ingest did not finish")
end

_as_json(call, method, path, body) = call(method, path;
    headers = ["Content-Type" => "application/json"], body = Vector{UInt8}(JSON3.write(body)))

@testset "setup mode routes" begin
    @testset "the folder preview reports each setup's mode and layout" begin
        auto = _as_write_folder!(mktempdir(), _AS_RUN_A; at = [("wat", 100)])
        cap = mktempdir(); mkpath(joinpath(cap, "data"))
        write(joinpath(cap, "data", "MO_S0001_0_001.tif"), "fake tif")
        db = fresh_db()
        with_inproc_routes(db) do call
            r = JSON3.read(String(call("GET", "/api/fs/preview?path=$(HimalayaUI.HTTP.escapeuri(auto))").body))
            @test only(r.setups).mode == "autosampler"
            @test only(r.setups).patterns.image == "{name}_0_001.tif"
            @test only(r.setups).patterns.integration == "{name}_tot.dat"
            @test only(r.setups).analysis_dir == joinpath(auto, "analysis", "automatic_analysis", "tot_files")
            @test only(r.setups).exposures == 3
            @test only(r.setups).without_trace == 0
            c = JSON3.read(String(call("GET", "/api/fs/preview?path=$(HimalayaUI.HTTP.escapeuri(cap))").body))
            @test only(c.setups).mode == "capillary"
        end
    end

    @testset "a new experiment takes the detected mode unless the request names one" begin
        auto = _as_write_folder!(mktempdir(), _AS_RUN_A; at = [("wat", 100)])
        db = fresh_db()
        with_inproc_routes(db) do call
            r = _as_json(call, "POST", "/api/experiments", Dict(:data_dir => joinpath(auto, "data")))
            @test r.status == 202
            eid = Int(JSON3.read(String(r.body)).id)
            _as_wait_ingest(db, eid)
            st = only(HimalayaUI.setups_of(db, eid))
            @test st.mode == "autosampler"
            @test st.integration_pattern == "{name}_tot.dat"
            @test _as_count(db, "SELECT COUNT(*) AS c FROM samples WHERE experiment_id = ?", eid) == 3
            shown = JSON3.read(String(call("GET", "/api/experiments/$eid").body))
            @test only(shown.setups).mode == "autosampler"

            other = _as_write_folder!(mktempdir(), _AS_RUN_A)
            r2 = _as_json(call, "POST", "/api/experiments", Dict(:data_dir => joinpath(other, "data"), :mode => "capillary"))
            eid2 = Int(JSON3.read(String(r2.body)).id)
            _as_wait_ingest(db, eid2)
            @test only(HimalayaUI.setups_of(db, eid2)).mode == "capillary"

            bad = _as_json(call, "POST", "/api/experiments",
                Dict(:data_dir => joinpath(_as_write_folder!(mktempdir(), _AS_RUN_A), "data"), :mode => "powder"))
            @test bad.status == 400
            @test occursin("mode", JSON3.read(String(bad.body)).error)
        end
    end

    @testset "a beamtime's setups and added setups carry a mode" begin
        bt = mktempdir()
        a = _as_write_folder!(joinpath(bt, "setupA"), _AS_RUN_A; at = [("wat", 100)])
        db = fresh_db()
        with_inproc_routes(db) do call
            spec(root; kw...) = Dict(:root => root, :data_dir => joinpath(root, "data"), kw...)
            r = _as_json(call, "POST", "/api/experiments", Dict(:path => bt, :setups => [spec(a)]))
            eid = Int(JSON3.read(String(r.body)).id)
            _as_wait_ingest(db, eid)
            @test only(HimalayaUI.setups_of(db, eid)).mode == "autosampler"

            b = _as_write_folder!(joinpath(bt, "setupB"), _AS_RUN_B; at = [("pep", 103)])
            @test _as_json(call, "POST", "/api/experiments/$eid/setups", Dict(:setups => [spec(b; mode = "powder")])).status == 400
            added = _as_json(call, "POST", "/api/experiments/$eid/setups", Dict(:setups => [spec(b; mode = "autosampler")]))
            @test added.status == 202
            sid = only(JSON3.read(String(added.body)).setup_ids)
            _as_wait_ingest(db, eid)
            @test HimalayaUI.setup_row(db, sid).mode == "autosampler"
            @test _as_count(db, "SELECT COUNT(*) AS c FROM samples WHERE setup_id = ?", sid) == 2
        end
    end

    @testset "a setup found by a rescan takes its mode from the folder" begin
        bt = mktempdir()
        a = _as_write_folder!(joinpath(bt, "setupA"), _AS_RUN_A; at = [("wat", 100)])
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = realpath(bt))
        HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = realpath(joinpath(a, "data")),
            analysis_dir = realpath(joinpath(a, "analysis")), root_dir = realpath(a), mode = "autosampler")
        _as_write_folder!(joinpath(bt, "setupB"), _AS_RUN_B; at = [("pep", 103)])
        added = HimalayaUI.add_new_setups!(db, eid)
        @test HimalayaUI.setup_row(db, only(added)).mode == "autosampler"
    end
end

function _as_snapshot(db, sid)
    [(x.id, x.name, x.load_id, x.slot_index) for x in Tables.rowtable(DBInterface.execute(db,
        "SELECT id, name, load_id, slot_index FROM samples WHERE setup_id = ? ORDER BY id", [sid]))]
end

_as_most_per_sample(db, sid) = _as_count(db, """
    SELECT COALESCE(MAX(n), 0) AS c FROM (SELECT COUNT(*) AS n FROM exposures
     WHERE setup_id = ? AND kind = 'file' GROUP BY sample_id)""", sid)

@testset "autosampler samples keep their place as files arrive" begin
    @testset "an acquisition with no scan number never joins another sample" begin
        root = _as_write_folder!(mktempdir(), vcat(_AS_RUN_A, _AS_RUN_B); at = [("wat", 100), ("pep", 103)])
        _as_add_stem!(root, "calibration")
        db, eid, sid = _as_live_row_db(root)
        HimalayaUI.scan_and_group!(db, eid)
        before = _as_snapshot(db, sid)
        @test length(before) == 6
        @test "calibration" in [x[2] for x in before]

        for a in (("pep", "01C", "S", 105), ("wat", "01A", "B", 106))
            _as_add_acquisition!(root, a)
        end
        _as_add_at!(root, "wat", 106)
        HimalayaUI.scan_and_group!(db, eid)
        after = _as_snapshot(db, sid)
        @test length(after) == 8
        @test issubset(before, after)
        @test _as_most_per_sample(db, sid) == 1
    end

    @testset "a metadata file that arrives late does not shift the others" begin
        root = _as_write_folder!(mktempdir(), _AS_RUN_A; at = [("wat", 100)])
        late = joinpath(root, "data", "wat_01B_S101_0_001.prp")
        held = read(late, String)
        rm(late)
        db, eid, sid = _as_live_row_db(root)
        HimalayaUI.scan_and_group!(db, eid)
        before = _as_snapshot(db, sid)
        @test length(before) == 2

        write(late, held)
        _as_add_acquisition!(root, ("wat", "01D", "S", 103))
        HimalayaUI.scan_and_group!(db, eid)
        after = _as_snapshot(db, sid)
        @test length(after) == 4
        @test issubset(before, after)
        @test _as_most_per_sample(db, sid) == 1
    end

    @testset "a run file that appears after its first acquisition was scanned corrupts nothing" begin
        root = _as_write_folder!(mktempdir(), _AS_RUN_A)
        db, eid, sid = _as_live_row_db(root)
        _as_add_at!(root, "wat", 100)
        HimalayaUI.scan_and_group!(db, eid)
        before = _as_snapshot(db, sid)
        @test _as_count(db, "SELECT COUNT(*) AS c FROM loads WHERE setup_id = ?", sid) == 1

        _as_add_at!(root, "wat", 103)
        _as_add_acquisition!(root, ("wat", "01A", "B", 103))
        _as_add_acquisition!(root, ("wat", "01B", "S", 104))
        HimalayaUI.scan_and_group!(db, eid)
        after = _as_snapshot(db, sid)
        @test length(after) == 5
        @test issubset(before, after)
        @test _as_most_per_sample(db, sid) == 1
        @test _as_count(db, "SELECT COUNT(*) AS c FROM loads WHERE setup_id = ?", sid) == 2
    end
end

@testset "an undecided or calibration-only setup is settled when the evidence arrives" begin
    @testset "no files means no mode" begin
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "x", path = "/x")
        empty_dir = mkpath(joinpath(mktempdir(), "data"))
        sid = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = empty_dir, analysis_dir = empty_dir)
        HimalayaUI.scan_and_group!(db, eid)
        @test ismissing(HimalayaUI.setup_row(db, sid).mode)
        write(joinpath(empty_dir, "MO_S0001_0_001.tif"), "fake tif")
        HimalayaUI.scan_and_group!(db, eid)
        @test HimalayaUI.setup_row(db, sid).mode == "capillary"
    end

    @testset "a calibration shot taken first does not lock the setup to capillary" begin
        root = mktempdir()
        _as_add_stem!(root, "agbe_S0090")
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = root)
        sid = HimalayaUI.add_setup_from_folder!(db, eid, root)
        HimalayaUI.scan_and_group!(db, eid)
        @test HimalayaUI.setup_row(db, sid).mode == "capillary"
        @test _as_count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", sid) == 1

        for a in _AS_RUN_A
            _as_add_acquisition!(root, a)
        end
        _as_add_at!(root, "wat", 100)
        HimalayaUI.scan_and_group!(db, eid)
        st = HimalayaUI.setup_row(db, sid)
        @test st.mode == "autosampler"
        @test st.analysis_dir == joinpath(root, "analysis", "automatic_analysis", "tot_files")
        @test st.integration_pattern == "{name}_tot.dat"
        @test _as_count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", sid) == 4
        @test _as_count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ? AND filename LIKE 'agbe%'", sid) == 1
        @test _as_most_per_sample(db, sid) == 1
        @test _as_count(db, """SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?
                                AND filename LIKE 'wat%' AND (trace_hash IS NOT NULL OR analysis_version IS NOT NULL)""", sid) == 0
    end
end

@testset "the setup mode migration also runs on a database from before setups existed" begin
    path = joinpath(mktempdir(), "legacy.db")
    legacy = SQLite.DB(path)
    HimalayaUI.create_schema!(legacy)
    add_legacy_experiment_columns!(legacy)
    DBInterface.execute(legacy, """INSERT INTO experiments (id, name, path, data_dir, analysis_dir)
        VALUES (7, 'used', '/bt/a/data', '/bt/a/data', '/bt/a/analysis')""")
    DBInterface.execute(legacy, """INSERT INTO experiments (id, name, path, data_dir, analysis_dir)
        VALUES (8, 'empty', '/bt/b/data', '/bt/b/data', '/bt/b/analysis')""")
    DBInterface.execute(legacy, "INSERT INTO samples (id, experiment_id, name) VALUES (700, 7, 'MO (S01P01)')")
    DBInterface.execute(legacy, """INSERT INTO exposures (id, experiment_id, sample_id, filename)
        VALUES (7000, 7, 700, 'MO_S0001_0_001')""")
    SQLite.close(legacy)

    db = HimalayaUI.open_db(path)
    @test only(HimalayaUI.setups_of(db, 7)).mode == "capillary"
    @test ismissing(only(HimalayaUI.setups_of(db, 8)).mode)
    SQLite.close(db)
end

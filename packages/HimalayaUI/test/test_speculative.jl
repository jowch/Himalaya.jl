using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using Himalaya

@testset "speculative.jl unit tests" begin
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="x")

    # Insert two manual peaks via peak_curations(kind='add')
    # at q=1.0 (ratio 1) and q=2.0 (ratio 2 of Lamellar)
    res = DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 1.0)", [e_id])
    p1 = Int(DBInterface.lastrowid(res))
    res = DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 2.0)", [e_id])
    p2 = Int(DBInterface.lastrowid(res))

    # 2-peak Lamellar
    new_id = HimalayaUI.insert_speculative_index!(db, e_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))

    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM indices WHERE id = ?", [new_id]))
    @test length(rows) == 1
    @test String(rows[1].kind)  == "speculative"
    @test String(rows[1].phase) == "Lamellar"
    @test rows[1].basis ≈ 1.0 atol=1e-9
    # 2-point fit through (1,1) and (2,2) is exact ⇒ R² = 1
    @test rows[1].r_squared ≈ 1.0 atol=1e-9
    @test rows[1].score > 0  # nonzero coverage × consistency

    ip_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM index_peaks WHERE index_id = ? ORDER BY ratio_position", [new_id]))
    @test length(ip_rows) == 2
    @test ip_rows[1].ratio_position == 1
    @test ip_rows[2].ratio_position == 2
    # Both are curation kind
    @test all(r -> String(r.peak_kind) == "curation", ip_rows)

    # Snap helper — uses peak_curations UNION ALL auto_peaks (via the route query)
    peak_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, q FROM peak_curations WHERE exposure_id = ? AND kind = 'add'
           UNION ALL SELECT id, q FROM auto_peaks WHERE exposure_id = ?""",
        [e_id, e_id]))
    snaps = HimalayaUI.compute_snap(peak_rows, Himalaya.Lamellar, 1.0, 1)
    # Ratio 1 anchored at q=1.0 means basis = 1.0; ratio 2 predicted at q=2.0 should snap to p2.
    snap_r2 = first(filter(s -> s.ratio_position == 2, snaps))
    @test snap_r2.suggested_peak_id == p2
    @test snap_r2.predicted_q ≈ 2.0 atol=1e-9
    # Higher ratios have no candidates in this synthetic peak set
    snap_r3 = first(filter(s -> s.ratio_position == 3, snaps))
    @test snap_r3.suggested_peak_id === nothing

    # resolve_phase
    @test HimalayaUI.resolve_phase("Pn3m")          === Himalaya.Pn3m
    @test HimalayaUI.resolve_phase("Himalaya.Pn3m") === Himalaya.Pn3m
    @test HimalayaUI.resolve_phase("NotAPhase")     === nothing
end

@testset "pipeline preserves speculative across re-analyze" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    # Pick first two auto peaks for a synthetic Lamellar speculative
    peaks = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 2", [e_id]))
    @test length(peaks) >= 2
    p1, p2 = Int(peaks[1].id), Int(peaks[2].id)

    new_id = HimalayaUI.insert_speculative_index!(db, e_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))

    # D-10: speculative indices survive re-analyze by virtue of kind='speculative'
    # (only kind='auto' indices are wiped/rebuilt), independent of any active-set
    # membership — so no group/assignment setup is needed to exercise this.

    pre_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM indices WHERE id = ?", [new_id]))
    pre_basis = Float64(pre_rows[1].basis)

    # Re-run analyze
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    # Index still exists with the same id (kind='speculative' preserved across wipe)
    post_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM indices WHERE id = ?", [new_id]))
    @test length(post_rows) == 1
    @test String(post_rows[1].kind) == "speculative"
    # Basis should be very close — same q values resolved to new auto peak ids
    @test post_rows[1].basis ≈ pre_basis rtol=0.01

    # index_peaks rows re-resolved with current peak ids
    ip_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT ip.*, COALESCE(ap.q, pc.q) AS q
           FROM index_peaks ip
           LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
           LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
           WHERE ip.index_id = ? ORDER BY ip.ratio_position""", [new_id]))
    @test length(ip_rows) == 2
    # (D-10: the speculative index is no longer auto-added to a custom group, so
    # there's no membership to assert here. The pipeline's custom-group re-attach
    # machinery itself is covered in test_pipeline.jl.)
end

# ── Helpers for the next two testsets ───────────────────────────────────────
# Build a synthetic exposure with curation-add peaks at known q-values so we can
# verify auto-discovery behaviour without depending on real-data peak shapes.
function _spec_synthetic_exposure(qs::Vector{Float64})
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="x")
    pids = Int[]
    for (i, q) in enumerate(qs)
        res = DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)",
            [e_id, q])
        push!(pids, Int(DBInterface.lastrowid(res)))
    end
    (; db, exposure_id = e_id, peak_ids = pids)
end

# Run the speculative-preservation half of `_persist_analysis_inner!` against
# a synthetic exposure (no .dat file). Mirrors what `analyze_exposure!` does
# minus auto peak detection. The `eff` NamedTuple is built from the curation
# rows so auto-discovery iterates over the right peaks.
function _spec_run_reanalyze!(db::SQLite.DB, exposure_id::Int)
    SQLite.transaction(db) do
        peaks_result = (q = Float64[], indices = Int[],
                        prominence = Float64[], sharpness = Float64[])
        I_full = Float64[]
        q_full = Float64[]
        # Build eff from current curation-add peaks (no auto peaks in synthetic exposure)
        add_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, q FROM peak_curations WHERE exposure_id = ? AND kind = 'add'",
            [exposure_id]))
        eff = (
            q         = Float64[Float64(r.q)  for r in add_rows],
            sharpness = Float64[0.0            for _ in add_rows],
            peak_id   = Int[Int(r.id)          for r in add_rows],
            peak_kind = Symbol[:curation       for _ in add_rows],
        )
        HimalayaUI._persist_analysis_inner!(db, exposure_id, q_full, I_full,
                                            peaks_result, Himalaya.Index[],
                                            Himalaya.Index[], eff)
    end
end

@testset "speculative auto-discovers new peaks during reanalysis" begin
    # Lamellar at basis = 0.05: predicted q's are 0.05, 0.10, 0.15, 0.20, ...
    # Place peaks at 0.05 (rp 1), 0.10 (rp 2), 0.15 (rp 3 — initially unseeded).
    fx = _spec_synthetic_exposure([0.05, 0.10, 0.15])
    p1, p2, p3 = fx.peak_ids[1], fx.peak_ids[2], fx.peak_ids[3]

    # Build the speculative with only rp 1 and rp 2 assigned.
    spec_id = HimalayaUI.insert_speculative_index!(fx.db, fx.exposure_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))
    @test length(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT * FROM index_peaks WHERE index_id = ?", [spec_id]))) == 2

    # Reanalyze — auto-discovery should snap p3 into rp 3.
    _spec_run_reanalyze!(fx.db, fx.exposure_id)

    ip = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT ratio_position, peak_id FROM index_peaks WHERE index_id = ? ORDER BY ratio_position",
        [spec_id]))
    @test length(ip) == 3
    @test [Int(r.ratio_position) for r in ip] == [1, 2, 3]
    @test Int(ip[3].peak_id) == p3
    rows = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT status, kind FROM indices WHERE id = ?", [spec_id]))
    @test String(rows[1].status) == "candidate"
    @test String(rows[1].kind)   == "speculative"
end

@testset "speculative revives from stale when matching peaks return" begin
    # Same Lamellar setup: peaks at 0.05, 0.10, 0.15.
    fx = _spec_synthetic_exposure([0.05, 0.10, 0.15])
    p1, p2, p3 = fx.peak_ids[1], fx.peak_ids[2], fx.peak_ids[3]
    spec_id = HimalayaUI.insert_speculative_index!(fx.db, fx.exposure_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))

    # Force stale: remove anchor peaks (delete curation-add rows for p1, p2)
    # so re-resolve drops them. With only p3 (which doesn't have an existing
    # assignment) the index goes stale.
    DBInterface.execute(fx.db,
        "DELETE FROM peak_curations WHERE id IN (?, ?)", [p1, p2])
    _spec_run_reanalyze!(fx.db, fx.exposure_id)
    @test String(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT status FROM indices WHERE id = ?", [spec_id]))[1].status) == "stale"
    @test length(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT * FROM index_peaks WHERE index_id = ?", [spec_id]))) == 0

    # Re-add the peaks. Snapshot is empty (index_peaks was wiped during the stale
    # cycle), so basis_for_snap falls back to the persisted `basis` on the
    # indices row. Auto-discovery then snaps all three peaks.
    res = DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.05)", [fx.exposure_id])
    new_p1 = Int(DBInterface.lastrowid(res))
    res = DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.10)", [fx.exposure_id])
    new_p2 = Int(DBInterface.lastrowid(res))
    _spec_run_reanalyze!(fx.db, fx.exposure_id)
    @test String(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT status FROM indices WHERE id = ?", [spec_id]))[1].status) == "candidate"
    ip = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT ratio_position FROM index_peaks WHERE index_id = ? ORDER BY ratio_position",
        [spec_id]))
    @test length(ip) >= 2
end

@testset "speculative create is atomic: both indices row and user_actions row exist" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    peaks = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 2", [e_id]))
    p1 = Int(peaks[1].id)
    p2 = Int(peaks[2].id)

    with_test_server(db) do port, base
        body = Dict(:phase => "Lamellar",
                    :anchor_peak_id => p1, :anchor_ratio => 1,
                    :additional => [Dict(:ratio_position => 2, :peak_id => p2)])
        r = HTTP.post("$base/api/exposures/$e_id/speculative";
            body = JSON3.write(body),
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"])
        @test r.status == 200
        new_ix = JSON3.read(String(r.body))
        new_id = Int(new_ix.id)

        # Both the indices row and the user_actions row must exist with matching index_id.
        idx_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM indices WHERE id = ? AND kind = 'speculative'", [new_id]))
        @test length(idx_rows) == 1

        evt_rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT id FROM user_actions
               WHERE action = 'speculative_created'
                 AND json_extract(payload, '\$.index_id') = ?""", [new_id]))
        @test length(evt_rows) == 1
    end
end

@testset "speculative delete is atomic: both indices row and user_actions row removed/created together" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    peaks = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 2", [e_id]))
    p1 = Int(peaks[1].id)
    p2 = Int(peaks[2].id)

    new_id = HimalayaUI.insert_speculative_index!(db, e_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))

    with_test_server(db) do port, base
        r = HTTP.delete("$base/api/indices/$new_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 200

        # indices row gone
        idx_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM indices WHERE id = ?", [new_id]))
        @test isempty(idx_rows)

        # speculative_deleted event was recorded with matching index_id
        evt_rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT id FROM user_actions
               WHERE action = 'speculative_deleted'
                 AND json_extract(payload, '\$.index_id') = ?""", [new_id]))
        @test length(evt_rows) == 1
    end
end

@testset "speculative HTTP routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    peaks = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 2", [e_id]))
    p1 = Int(peaks[1].id)
    p2 = Int(peaks[2].id)

    with_test_server(db) do port, base
        # Snap endpoint
        r = HTTP.get("$base/api/exposures/$e_id/speculative-snap?phase=Lamellar&anchor_peak_id=$p1&anchor_ratio=1")
        @test r.status == 200
        snaps = JSON3.read(String(r.body))
        @test length(snaps) >= 2
        @test snaps[1].ratio_position == 1
        @test snaps[1].is_anchor === true

        # Bad phase
        r = HTTP.get("$base/api/exposures/$e_id/speculative-snap?phase=Bogus&anchor_peak_id=$p1";
                     status_exception = false)
        @test r.status == 400

        # Missing anchor
        r = HTTP.get("$base/api/exposures/$e_id/speculative-snap?phase=Lamellar";
                     status_exception = false)
        @test r.status == 400

        # Create speculative
        body = Dict(:phase => "Lamellar",
                    :anchor_peak_id => p1, :anchor_ratio => 1,
                    :additional => [Dict(:ratio_position => 2, :peak_id => p2)],
                    :active => true)
        r = HTTP.post("$base/api/exposures/$e_id/speculative";
            body = JSON3.write(body),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 200
        new_ix = JSON3.read(String(r.body))
        @test new_ix.kind == "speculative"
        @test new_ix.phase == "Lamellar"
        @test length(new_ix.peaks) == 2
        new_id = Int(new_ix.id)

        # D-10: `active:true` is accepted but no longer auto-adds to any set (the
        # legacy custom group + GET /groups are retired); "make active" is now an
        # explicit POST /assignment/members, covered in test_assignments.jl.

        # DELETE rejects auto indices
        auto_ix_id = Int(first(filter(ix -> ix.kind != "speculative",
                                       JSON3.read(String(HTTP.get("$base/api/exposures/$e_id/indices").body)))).id)
        r = HTTP.delete("$base/api/indices/$auto_ix_id"; status_exception = false)
        @test r.status == 403

        # DELETE allows speculative + cleans up
        r = HTTP.delete("$base/api/indices/$new_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        r2 = HTTP.get("$base/api/indices/$new_id"; status_exception = false)
        @test r2.status == 404
        # Group membership and index_peaks gone
        rs = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM index_group_members WHERE index_id = ?", [new_id]))
        @test isempty(rs)
        rs = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM index_peaks WHERE index_id = ?", [new_id]))
        @test isempty(rs)
    end
end

@testset "insert_speculative_index! inherits exposure analysis_inputs_hash (issue #35 Bug 6)" begin
    # Pre-fix, the INSERT into `indices` set inputs_hash = NULL on speculative
    # rows. StaleIndicesBanner gates on (index.inputs_hash !== exposure
    # .analysis_inputs_hash), so a speculative create immediately registered as
    # stale and the banner spuriously fired after every speculative create. The
    # fix in src/speculative.jl reads the exposure's current analysis_inputs_hash
    # and writes it onto the new index — they share the effective peak set, so
    # the inherited hash is correct by construction.
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="x")

    expected_hash = "deadbeef" ^ 8  # 64 hex chars, matches SHA-256 fingerprint shape
    DBInterface.execute(db,
        "UPDATE exposures SET analysis_inputs_hash = ? WHERE id = ?",
        [expected_hash, e_id])

    res = DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 1.0)", [e_id])
    p1 = Int(DBInterface.lastrowid(res))
    res = DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 2.0)", [e_id])
    p2 = Int(DBInterface.lastrowid(res))

    new_id = HimalayaUI.insert_speculative_index!(db, e_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))

    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT inputs_hash FROM indices WHERE id = ?", [new_id]))
    @test length(rows) == 1
    @test String(rows[1].inputs_hash) == expected_hash
end

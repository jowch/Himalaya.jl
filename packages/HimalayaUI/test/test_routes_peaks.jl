using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "peaks routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        # GET peaks
        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) > 0
        @test all(p -> p.source == "auto", list)

        initial_indices = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM indices WHERE exposure_id = ?", [e_id]))
        @test !isempty(initial_indices)

        # POST manual peak
        r = HTTP.post("$base/api/exposures/$e_id/peaks";
            body = JSON3.write(Dict(:q => 0.5)),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 201
        body = JSON3.read(String(r.body))
        @test body.q == 0.5
        @test body.source == "manual"
        peak_id = body.id

        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        list = JSON3.read(String(r.body))
        manual = filter(p -> p.source == "manual", list)
        @test length(manual) == 1

        # DELETE peak
        r = HTTP.delete("$base/api/peaks/$peak_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        list = JSON3.read(String(r.body))
        # After deleting the only manual peak, no manual peaks should remain.
        # (Note: auto peak ids and curation ids are in separate sequences and
        #  may collide; check by source rather than by id.)
        @test !any(p -> p.source == "manual", list)

        # 404 delete
        r = HTTP.delete("$base/api/peaks/99999";
            headers = ["X-Username" => "alice"], status_exception = false)
        @test r.status == 404

        # ── PATCH excluded on auto peaks ─────────────────────────────────
        all_peaks = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 1",
            [e_id]))
        auto_id = Int(all_peaks[1].id)

        # Toggle excluded=true
        r = HTTP.patch("$base/api/peaks/$auto_id";
            body = JSON3.write(Dict(:excluded => true)),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.excluded == true
        @test body.source   == "auto"

        # Toggle back to false
        r = HTTP.patch("$base/api/peaks/$auto_id";
            body = JSON3.write(Dict(:excluded => false)),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.excluded == false

        # PATCHing a manual peak's excluded → 400 (manual peaks delete instead)
        r = HTTP.post("$base/api/exposures/$e_id/peaks";
            body = JSON3.write(Dict(:q => 0.42)),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        manual_id = JSON3.read(String(r.body)).id
        r = HTTP.patch("$base/api/peaks/$manual_id";
            body = JSON3.write(Dict(:excluded => true)),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            status_exception = false)
        @test r.status == 400

        # PATCH unknown peak → 404
        r = HTTP.patch("$base/api/peaks/99999";
            body = JSON3.write(Dict(:excluded => true)),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            status_exception = false)
        @test r.status == 404

        # ── Reanalysis preserves excluded state ──────────────────────────
        # Re-mark the auto peak as excluded, then re-run analysis.
        r = HTTP.patch("$base/api/peaks/$auto_id";
            body = JSON3.write(Dict(:excluded => true)),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 200
        excluded_q = JSON3.read(String(r.body)).q

        HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

        # After reanalysis, the auto peak at the same q still exists in auto_peaks,
        # and an exclude curation row should still be present for it.
        # Use the joined view (get_peaks_for_exposure) via GET /api/exposures/:id/peaks.
        r_after = HTTP.get("$base/api/exposures/$e_id/peaks")
        post_list = JSON3.read(String(r_after.body))
        match = filter(p -> abs(Float64(p.q) - excluded_q) < 1e-6, post_list)
        @test !isempty(match)
        @test match[1].excluded == true
    end
end

# ── New R2.2 curation-lifecycle testsets ─────────────────────────────────────

@testset "POST /peaks then GET returns manual peak with source='manual'" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        r = HTTP.post("$base/api/exposures/$e_id/peaks";
            body = JSON3.write(Dict(:q => 0.77)),
            headers = ["Content-Type" => "application/json", "X-Username" => "bob"])
        @test r.status == 201
        new_id = JSON3.read(String(r.body)).id

        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        list = JSON3.read(String(r.body))
        manual = filter(p -> p.source == "manual", list)
        @test length(manual) == 1
        @test manual[1].id == new_id
        @test manual[1].q ≈ 0.77
        @test manual[1].excluded == false
    end
end

@testset "PATCH excluded=true inserts curation; excluded=false removes it" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        # Pick an auto peak
        auto_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 1", [e_id]))
        auto_id = Int(auto_rows[1].id)
        auto_q  = Float64(auto_rows[1].q)

        # Exclude it
        r = HTTP.patch("$base/api/peaks/$auto_id";
            body = JSON3.write(Dict(:excluded => true)),
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).excluded == true

        # Verify a curation row was inserted
        excl = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude' AND ABS(q - ?) < 1e-6",
            [e_id, auto_q]))
        @test length(excl) == 1

        # Un-exclude it
        r = HTTP.patch("$base/api/peaks/$auto_id";
            body = JSON3.write(Dict(:excluded => false)),
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).excluded == false

        # Curation row should be gone
        excl2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude' AND ABS(q - ?) < 1e-6",
            [e_id, auto_q]))
        @test isempty(excl2)
    end
end

@testset "PATCH excluded=true is idempotent (no duplicate curation rows)" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        auto_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 1", [e_id]))
        auto_id = Int(auto_rows[1].id)
        auto_q  = Float64(auto_rows[1].q)

        # Two consecutive PATCH { excluded: true }
        for _ in 1:2
            r = HTTP.patch("$base/api/peaks/$auto_id";
                body = JSON3.write(Dict(:excluded => true)),
                headers = ["Content-Type" => "application/json", "X-Username" => "alice"])
            @test r.status == 200
        end

        # Must have exactly one exclude curation row
        excl = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude' AND ABS(q - ?) < 1e-6",
            [e_id, auto_q]))
        @test length(excl) == 1
    end
end

@testset "DELETE manual peak removes curation row; subsequent GET doesn't return it" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        # Add a manual peak
        r = HTTP.post("$base/api/exposures/$e_id/peaks";
            body = JSON3.write(Dict(:q => 0.88)),
            headers = ["Content-Type" => "application/json", "X-Username" => "carol"])
        @test r.status == 201
        manual_id = JSON3.read(String(r.body)).id

        # Delete it
        r = HTTP.delete("$base/api/peaks/$manual_id";
            headers = ["X-Username" => "carol"])
        @test r.status == 204

        # Should no longer appear in GET (no manual peaks)
        r = HTTP.get("$base/api/exposures/$e_id/peaks")
        list = JSON3.read(String(r.body))
        @test !any(p -> p.source == "manual", list)

        # curation row gone
        cr = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peak_curations WHERE id = ?", [manual_id]))
        @test isempty(cr)
    end
end

@testset "POST /peaks response id round-trips to the correct curation row" begin
    # Regression for the read-back race: the response id must point to a
    # peak_curations row whose q matches the request, not a concurrent add.
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        # Two sequential POSTs at different q values.
        r1 = HTTP.post("$base/api/exposures/$e_id/peaks";
            body = JSON3.write(Dict(:q => 0.111)),
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"])
        @test r1.status == 201
        b1 = JSON3.read(String(r1.body))

        r2 = HTTP.post("$base/api/exposures/$e_id/peaks";
            body = JSON3.write(Dict(:q => 0.222)),
            headers = ["Content-Type" => "application/json", "X-Username" => "bob"])
        @test r2.status == 201
        b2 = JSON3.read(String(r2.body))

        # Each response id must point to the row whose q matches the request.
        row1 = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, q FROM peak_curations WHERE id = ?", [Int(b1.id)]))
        @test length(row1) == 1
        @test Float64(row1[1].q) ≈ 0.111 atol=1e-9

        row2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, q FROM peak_curations WHERE id = ?", [Int(b2.id)]))
        @test length(row2) == 1
        @test Float64(row2[1].q) ≈ 0.222 atol=1e-9

        # The two ids must be distinct.
        @test Int(b1.id) != Int(b2.id)
    end
end

@testset "DELETE auto peak returns 400" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_test_server(db) do port, base
        auto_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 1", [e_id]))
        auto_id = Int(auto_rows[1].id)

        r = HTTP.delete("$base/api/peaks/$auto_id";
            headers = ["X-Username" => "alice"],
            status_exception = false)
        @test r.status == 400
    end
end

# ── M2.2 BACKEND: with_idempotency + sync reanalyze + post_state broadcast ──

function _setup_peaks_test(tmp)
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    return db, e_id, analysis_dir
end

function _attach_sse_subscriber!()
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    return pending, sub
end

function _detach_sse_subscriber!(sub, pending)
    lock(HimalayaUI.SSE_LOCK) do
        filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
    end
    close(pending)
end

function _drain_curation_frames(pending; timeout=0.5)
    frames = String[]
    deadline = time() + timeout
    while time() < deadline
        if isready(pending)
            push!(frames, take!(pending))
        else
            sleep(0.01)
        end
    end
    # Filter to curation-only frames (subscriber init may emit other types).
    filter(f -> occursin("event: curation", f), frames)
end

function _frame_kind(frame::String)
    m = match(r"\"kind\":\"([^\"]+)\"", frame)
    m === nothing ? nothing : m.captures[1]
end

@testset "M2.2: POST /peaks emits one enriched SSE frame with post_state" begin
    tmp = mktempdir()
    db, e_id, _ = _setup_peaks_test(tmp)

    with_test_server(db) do port, base
        pending, sub = _attach_sse_subscriber!()
        try
            # Drain any existing frames (e.g. analyze_run from setup).
            _drain_curation_frames(pending; timeout=0.05)

            r = HTTP.post("$base/api/exposures/$e_id/peaks";
                body = JSON3.write(Dict(:q => 0.55)),
                headers = ["Content-Type"  => "application/json",
                           "X-Username"    => "alice",
                           "X-Client-Op-Id" => "uuid-m22-post-1"])
            @test r.status == 201
            body = JSON3.read(String(r.body))
            @test haskey(body, :event_id)
            @test haskey(body, :view_row_id)
            @test haskey(body, :analysis_inputs_hash)

            db_hash = HimalayaUI.read_inputs_hash(db, e_id)
            @test String(body.analysis_inputs_hash) == db_hash

            # Exactly one enriched curation frame for peak_added.
            frames = _drain_curation_frames(pending)
            peak_added = filter(f -> _frame_kind(f) == "peak_added", frames)
            @test length(peak_added) == 1
            frame = peak_added[1]
            @test occursin("\"client_op_id\":\"uuid-m22-post-1\"", frame)
            @test occursin("\"post_state\"", frame)
            @test occursin("\"analysis_inputs_hash\"", frame)
            @test occursin("\"indices\"", frame)
            # No analyze_run frame (defer_broadcast suppressed it).
            @test !any(f -> _frame_kind(f) == "analyze_run", frames)
        finally
            _detach_sse_subscriber!(sub, pending)
        end
    end
end

@testset "M2.2: POST /peaks retry with same op_id replays cache, no dup events/broadcasts" begin
    tmp = mktempdir()
    db, e_id, _ = _setup_peaks_test(tmp)

    with_test_server(db) do port, base
        pending, sub = _attach_sse_subscriber!()
        try
            _drain_curation_frames(pending; timeout=0.05)

            r1 = HTTP.post("$base/api/exposures/$e_id/peaks";
                body = JSON3.write(Dict(:q => 0.66)),
                headers = ["Content-Type"   => "application/json",
                           "X-Username"     => "alice",
                           "X-Client-Op-Id" => "uuid-m22-retry"])
            @test r1.status == 201
            r1_body_str = String(copy(r1.body))

            frames1 = _drain_curation_frames(pending)
            @test length(filter(f -> _frame_kind(f) == "peak_added", frames1)) == 1

            # Same op_id → cached body returned, no new event row, no broadcast.
            r2 = HTTP.post("$base/api/exposures/$e_id/peaks";
                body = JSON3.write(Dict(:q => 0.66)),
                headers = ["Content-Type"   => "application/json",
                           "X-Username"     => "alice",
                           "X-Client-Op-Id" => "uuid-m22-retry"])
            @test r2.status == 201
            @test r1_body_str == String(copy(r2.body))

            # Drain — should be empty.
            frames2 = _drain_curation_frames(pending; timeout=0.2)
            @test isempty(filter(f -> _frame_kind(f) == "peak_added", frames2))

            # Exactly one peak_added user_actions row.
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM user_actions WHERE action = 'peak_added' AND entity_id = ?",
                [e_id]))
            @test length(rows) == 1
        finally
            _detach_sse_subscriber!(sub, pending)
        end
    end
end

@testset "M2.2: PATCH excluded uses with_idempotency and runs analyze synchronously" begin
    tmp = mktempdir()
    db, e_id, _ = _setup_peaks_test(tmp)

    with_test_server(db) do port, base
        pending, sub = _attach_sse_subscriber!()
        try
            auto_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 1", [e_id]))
            auto_id = Int(auto_rows[1].id)

            hash_before = HimalayaUI.read_inputs_hash(db, e_id)
            _drain_curation_frames(pending; timeout=0.05)

            r = HTTP.patch("$base/api/peaks/$auto_id";
                body = JSON3.write(Dict(:excluded => true)),
                headers = ["Content-Type"   => "application/json",
                           "X-Username"     => "alice",
                           "X-Client-Op-Id" => "uuid-m22-patch"])
            @test r.status == 200
            r_body_str = String(copy(r.body))
            body = JSON3.read(r_body_str)
            @test body.excluded == true
            @test haskey(body, :event_id)
            @test haskey(body, :analysis_inputs_hash)

            hash_after = HimalayaUI.read_inputs_hash(db, e_id)
            @test hash_after !== nothing
            @test hash_after != hash_before  # exclusion changed effective peak set

            frames = _drain_curation_frames(pending)
            excl = filter(f -> _frame_kind(f) == "peak_excluded", frames)
            @test length(excl) == 1
            @test occursin("\"post_state\"", excl[1])

            # Retry returns cached body, no second event.
            r2 = HTTP.patch("$base/api/peaks/$auto_id";
                body = JSON3.write(Dict(:excluded => true)),
                headers = ["Content-Type"   => "application/json",
                           "X-Username"     => "alice",
                           "X-Client-Op-Id" => "uuid-m22-patch"])
            @test r2.status == 200
            @test r_body_str == String(copy(r2.body))

            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM user_actions WHERE action = 'peak_excluded' AND entity_id = ?",
                [e_id]))
            @test length(rows) == 1
        finally
            _detach_sse_subscriber!(sub, pending)
        end
    end
end

@testset "M2.2: DELETE manual peak emits peak_removed with post_state" begin
    tmp = mktempdir()
    db, e_id, _ = _setup_peaks_test(tmp)

    with_test_server(db) do port, base
        pending, sub = _attach_sse_subscriber!()
        try
            # First add a manual peak (without subscriber drained).
            r = HTTP.post("$base/api/exposures/$e_id/peaks";
                body = JSON3.write(Dict(:q => 0.42)),
                headers = ["Content-Type" => "application/json", "X-Username" => "alice"])
            @test r.status == 201
            manual_id = JSON3.read(String(r.body)).id
            _drain_curation_frames(pending; timeout=0.1)

            r = HTTP.delete("$base/api/peaks/$manual_id";
                headers = ["X-Username"     => "alice",
                           "X-Client-Op-Id" => "uuid-m22-del"])
            @test r.status == 204

            frames = _drain_curation_frames(pending)
            removed = filter(f -> _frame_kind(f) == "peak_removed", frames)
            @test length(removed) == 1
            @test occursin("\"post_state\"", removed[1])
            @test occursin("\"analysis_inputs_hash\"", removed[1])
        finally
            _detach_sse_subscriber!(sub, pending)
        end
    end
end

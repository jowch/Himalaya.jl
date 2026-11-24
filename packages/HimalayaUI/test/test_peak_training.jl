using Test, Random, SQLite, DBInterface, Tables
import JSON
using HimalayaUI
import HimalayaPeakModel as HPM

isdefined(@__MODULE__, :_pt_db) || include("peak_training_fixture.jl")

@testset "cold start" begin
    mktempdir() do dir
        b = _pt_db(dir; analyze = false)
        c = HimalayaUI.cold_start_check(b.db)
        @test c.ok && c.group_by === :experiment
        @test c.counts.curations == 36 * 4 + 12 && c.counts.curated_exposures == 36
        @test c.counts.experiments == 3 && c.counts.samples == 12
        @test isempty(c.reasons)

        c = HimalayaUI.cold_start_check(b.db; min_experiments = 4)
        @test c.ok && c.group_by === :sample
        c = HimalayaUI.cold_start_check(b.db; min_experiments = 4, min_samples = 13)
        @test !c.ok && length(c.reasons) == 1 && occursin("experiments", c.reasons[1]) && occursin("samples", c.reasons[1])
        c = HimalayaUI.cold_start_check(b.db; min_curations = 200, min_exposures = 37)
        @test !c.ok && length(c.reasons) == 2
        @test any(r -> occursin("200", r) && occursin("156", r), c.reasons)
        @test any(r -> occursin("37", r) && occursin("36", r), c.reasons)

        DBInterface.execute(b.db, "DELETE FROM peak_curations WHERE kind = 'add'")
        c = HimalayaUI.cold_start_check(b.db)
        @test !c.ok && c.counts.curations == 12 && c.counts.curated_exposures == 12
        DBInterface.close!(b.db)
    end
end

@testset "gate" begin
    ok = (peaks = 5.0, ff_peaks = 1.0, phase = (10, 20), adds = (50, 80), excluded = (3, 10))
    hy = (peaks = 5.0, ff_peaks = 1.0, phase = (10, 20), adds = (40, 80), excluded = (5, 10))
    g = HimalayaUI.evaluate_gate(ok, hy)
    @test g.passed && isempty(g.reasons) && all(c -> c.passed, g.checks)
    @test Set(c.name for c in g.checks) == Set(["peaks_per_exposure", "peaks_per_form_factor", "phase_calls", "add_recall", "excludes_detected"])

    bad(; kw...) = HimalayaUI.evaluate_gate(merge(ok, NamedTuple(kw)), hy)
    @test HimalayaUI.evaluate_gate((; ok..., peaks = 5.5), hy).passed
    g = bad(peaks = 5.51); @test !g.passed && length(g.reasons) == 1 && occursin("peaks per exposure", g.reasons[1])
    g = bad(ff_peaks = 1.11); @test !g.passed && occursin("form-factor", g.reasons[1])
    g = bad(phase = (9, 20)); @test !g.passed && occursin("phase", g.reasons[1])
    @test bad(phase = (11, 20)).passed
    g = bad(adds = (39, 80)); @test !g.passed && occursin("add recall", g.reasons[1])
    g = bad(excluded = (6, 10)); @test !g.passed && occursin("exclude", g.reasons[1])
    g = bad(peaks = 9.0, adds = (1, 80)); @test length(g.reasons) == 2

    none = (; ok..., phase = (0, 0), adds = (0, 0), excluded = (0, 0))
    g = HimalayaUI.evaluate_gate(none, (; hy..., phase = (0, 0), adds = (0, 0), excluded = (0, 0)))
    @test g.passed
    @test all(c -> c.name in ("peaks_per_exposure", "peaks_per_form_factor") || c.skipped, g.checks)
    g = HimalayaUI.evaluate_gate((; ok..., ff_peaks = NaN), (; hy..., ff_peaks = NaN))
    @test g.passed
    g = HimalayaUI.evaluate_gate(ok, hy; peaks_slack = 0.0, ff_slack = 0.0)
    @test g.passed
    @test !HimalayaUI.evaluate_gate((; ok..., peaks = 5.1), hy; peaks_slack = 0.0).passed
end

@testset "top call matching" begin
    q = [0.07, 0.14, 0.21, 0.28]
    sh = fill(1.0, 4)
    @test HimalayaUI.top_call_matches(q, sh, [("Lamellar", 0.07)])
    @test HimalayaUI.top_call_matches(q, sh, [("Hexagonal", 0.07), ("Lamellar", 0.0705)])
    @test !HimalayaUI.top_call_matches(q, sh, [("Lamellar", 0.075)])
    @test !HimalayaUI.top_call_matches(q, sh, [("Pn3m", 0.07)])
    @test !HimalayaUI.top_call_matches(Float64[], Float64[], [("Lamellar", 0.07)])
    @test !HimalayaUI.top_call_matches([0.07], [1.0], [("Lamellar", 0.07)])
end

@testset "confirmed calls and shown-model filter" begin
    mktempdir() do dir
        b = _pt_db(dir; experiments = 1, samples = 1, per = 2, analyze = false)
        e1, e2 = b.eids
        calls = HimalayaUI.confirmed_calls(b.db, b.eids)
        @test Set(keys(calls)) == Set(b.eids) && calls[e1][1][1] == "Lamellar"
        DBInterface.execute(b.db, "UPDATE assignments SET state = 'form_factor' WHERE exposure_id = ?", [e2])
        @test collect(keys(HimalayaUI.confirmed_calls(b.db, b.eids))) == [e1]
        DBInterface.execute(b.db, "INSERT OR REPLACE INTO sample_calls (sample_id, state) VALUES (?, 'null')", [b.sids[1]])
        @test isempty(HimalayaUI.confirmed_calls(b.db, b.eids))
        DBInterface.execute(b.db, "UPDATE indices SET kind = 'custom'")
        DBInterface.execute(b.db, "DELETE FROM sample_calls")
        DBInterface.execute(b.db, "UPDATE assignments SET state = 'indexed'")
        @test isempty(HimalayaUI.confirmed_calls(b.db, b.eids))

        @test HimalayaUI.no_model_shown_ids(b.db) == Set(b.eids)
        DBInterface.execute(b.db, "UPDATE peak_reviews SET shown_model_id = 7 WHERE exposure_id = ?", [e1])
        @test HimalayaUI.no_model_shown_ids(b.db) == Set([e2])
        HimalayaUI.record_review!(b.db, e1; direct = true)
        @test HimalayaUI.no_model_shown_ids(b.db) == Set(b.eids)
        DBInterface.close!(b.db)
    end
end

@testset "confirmed calls: speculative indices always count, auto ones need a user assignment event" begin
    mktempdir() do dir
        b = _pt_db(dir; experiments = 1, samples = 1, per = 2, analyze = false)
        e1, e2 = b.eids
        DBInterface.execute(b.db, "DELETE FROM user_actions")
        @test isempty(HimalayaUI.confirmed_calls(b.db, b.eids))
        @test Set(keys(HimalayaUI.confirmed_calls(b.db, b.eids; auto = :all))) == Set(b.eids)

        DBInterface.execute(b.db,
            "INSERT INTO user_actions (action, entity_type, entity_id) VALUES ('assignment_set_state', 'exposure', ?)", [e1])
        DBInterface.execute(b.db,
            "INSERT INTO user_actions (action, entity_type, entity_id) VALUES ('assignment_add', 'sample', ?)", [e2])
        @test collect(keys(HimalayaUI.confirmed_calls(b.db, b.eids))) == [e1]

        DBInterface.execute(b.db, "UPDATE indices SET kind = 'speculative' WHERE exposure_id = ?", [e2])
        calls = HimalayaUI.confirmed_calls(b.db, b.eids)
        @test Set(keys(calls)) == Set(b.eids) && calls[e2][1][1] == "Lamellar"
        DBInterface.execute(b.db, "UPDATE assignments SET state = 'null' WHERE exposure_id = ?", [e2])
        @test collect(keys(HimalayaUI.confirmed_calls(b.db, b.eids))) == [e1]
        DBInterface.close!(b.db)
    end
end

@testset "run_peak_training on a read-only synthetic database" begin
    mktempdir() do dir
        b = _pt_db(dir)
        DBInterface.close!(b.db)
        ro = HimalayaUI.open_readonly_db(b.path)
        @test_throws SQLite.SQLiteException DBInterface.execute(ro, "DELETE FROM peak_curations")
        cache_dir = HimalayaUI.peak_feature_cache_dir(b.path)
        @test cache_dir == joinpath(dir, "cache", "peak-features", "h")
        @test HimalayaUI.peak_feature_cache_dir(joinpath(dir, "other.db")) != cache_dir
        recipe = HPM.default_recipe(n_trees = 25)

        r = HimalayaUI.run_peak_training(ro; recipe, cache_dir, io = devnull,
                                         peaks_slack = 100, ff_slack = 100)
        r = JSON.parse(JSON.json(r))
        @test r.status == "passed"
        @test r.counts.exposures == 36 && r.counts.labelled_exposures == 36 && r.counts.group_by == "experiment"
        @test r.counts.labels.positive > 0 && r.counts.labels.negative + r.counts.labels.weak_negative > 0
        @test r.watermark == HimalayaUI.label_watermark(ro)
        @test r.core_version == string(pkgversion(HimalayaUI.Himalaya))
        @test r.feature_version == HPM.feature_version(recipe.features)
        @test r.recipe_id == HPM.recipe_id(HPM.from_dict(HPM.Recipe, r.recipe))
        @test r.threshold in 0.3:0.05:0.9 && r.nms == 8
        @test length(readdir(cache_dir)) == 36
        cmp = r.metrics.comparison
        @test cmp.hybrid.exposures == 36 && cmp.model.exposures == 36
        @test cmp.phase_auto == "user_events"
        @test cmp.phase_targets == 36 && cmp.model.phase[2] == 36 && cmp.hybrid.phase[2] == 36
        @test cmp.phase_targets_mode == "reachable"
        @test cmp.phase_target_counts.all == 36 && cmp.phase_target_counts.reachable == 36
        @test cmp.phase_target_counts.per_phase.Lamellar.reachable == 36
        @test cmp.model.adds[2] == cmp.hybrid.adds[2] > 0
        @test cmp.model.adds_in_band == cmp.hybrid.adds_in_band
        @test r.metrics.gate.passed == (r.status == "passed")
        @test length(r.metrics.sweep) == length(0.3:0.05:0.9)
        @test r.metrics.cv.candidate_auc > 0.5
        @test length(r.metrics.cv.by_threshold) == 5
        @test sort(collect(String, keys(r.metrics.gate.caps))) == ["peaks_per_exposure", "peaks_per_form_factor"]
        if r.status == "passed"
            f = HPM.from_dict(HPM.Fitted, r.fitted)
            @test HPM.canonical(HPM.to_dict(f)) == HPM.canonical(r.fitted)
            @test isempty(r.metrics.gate.reasons)
        else
            @test r.fitted === nothing && !isempty(r.metrics.gate.reasons)
        end

        again = JSON.parse(JSON.json(HimalayaUI.run_peak_training(ro; recipe, cache_dir, io = devnull,
                                                                    peaks_slack = 100, ff_slack = 100)))
        @test again.metrics.comparison == r.metrics.comparison && again.threshold == r.threshold
        @test length(readdir(cache_dir)) == 36

        tight = JSON.parse(JSON.json(HimalayaUI.run_peak_training(ro; recipe, io = devnull,
                                                                    peaks_slack = -100)))
        @test tight.status == "rejected" && tight.fitted === nothing && tight.threshold == 0.9
        @test "no threshold keeps peaks within the caps" in tight.metrics.gate.reasons
        DBInterface.close!(ro)
    end
end

@testset "run_peak_training refuses to train below the cold-start minimum" begin
    mktempdir() do dir
        b = _pt_db(dir; experiments = 2, samples = 2, per = 3, analyze = false)
        DBInterface.close!(b.db)
        ro = HimalayaUI.open_readonly_db(b.path)
        cache_dir = HimalayaUI.peak_feature_cache_dir(b.path)
        r = JSON.parse(JSON.json(HimalayaUI.run_peak_training(ro; cache_dir, io = devnull)))
        @test r.status == "rejected" && r.fitted === nothing
        @test !r.metrics.gate.passed && length(r.metrics.gate.reasons) == 3
        @test r.counts.curated_exposures == 12 && r.counts.experiments == 2
        @test r.recipe_id == HPM.recipe_id(HPM.default_recipe())
        @test !isdir(cache_dir)
        DBInterface.close!(ro)
    end
end

@testset "peak_eval keeps per-exposure phase rows alongside the summary" begin
    mktempdir() do dir
        b = _pt_db(dir)
        DBInterface.close!(b.db)
        ro = HimalayaUI.open_readonly_db(b.path)
        exps = HimalayaUI.labeled_exposures(ro; include_unlabelled = true, io = devnull)
        ctx = HimalayaUI.eval_context(ro, exps)
        @test length(ctx.label_ids) == 36 && length(ctx.targets) == 36
        picks = [HimalayaUI.Himalaya.findpeaks(e.q, e.I, e.σ; method = :hybrid).indices for e in exps]
        ev = HimalayaUI.peak_eval(exps, picks; ctx.targets, ctx.label_ids)
        @test isequal(ev.summary, HimalayaUI.peak_summary(exps, picks; ctx.targets, ctx.label_ids))
        @test length(ev.phase_rows) == ev.summary.phase[2] == 36
        @test count(r -> r.ok, ev.phase_rows) == ev.summary.phase[1]
        @test all(r -> r.phase == "Lamellar", ev.phase_rows)
        @test Set(r.id for r in ev.phase_rows) == Set(b.eids)
        @test ev.peak_count == sum(length, picks) && ev.ff_peak_count == 0
        none = HimalayaUI.peak_eval(exps, [Int[] for _ in exps]; ctx.targets, ctx.label_ids)
        @test none.summary.phase == (0, 36) && all(r -> !r.ok, none.phase_rows)
        DBInterface.close!(ro)
    end
end

@testset "phase targets keep only calls reachable from the user's own peaks" begin
    mktempdir() do dir
        b = _pt_db(dir; experiments = 1, samples = 1, per = 3)
        e1, e2, e3 = b.eids
        DBInterface.execute(b.db, "UPDATE indices SET basis = basis * 2 WHERE exposure_id = ?", [e2])
        DBInterface.execute(b.db, "UPDATE indices SET phase = 'Hexagonal' WHERE exposure_id = ?", [e3])
        DBInterface.close!(b.db)
        ro = HimalayaUI.open_readonly_db(b.path)
        exps = HimalayaUI.labeled_exposures(ro; include_unlabelled = true, io = devnull)

        ctx = HimalayaUI.eval_context(ro, exps)
        @test Set(keys(ctx.targets)) == Set([e1, e2])
        @test ctx.counts.all == 3 && ctx.counts.reachable == 2
        @test ctx.counts.per_phase["Lamellar"] == (reachable = 2, all = 2)
        @test ctx.counts.per_phase["Hexagonal"] == (reachable = 0, all = 1)
        @test ctx.label_ids == Set(b.eids)

        every = HimalayaUI.eval_context(ro, exps; phase_targets = :all)
        @test Set(keys(every.targets)) == Set(b.eids)
        @test every.counts == ctx.counts
        @test_throws ArgumentError HimalayaUI.eval_context(ro, exps; phase_targets = :typo)

        picks = HimalayaUI.oracle_picks(exps)
        ev = HimalayaUI.peak_eval(exps, picks; ctx.targets, ctx.label_ids)
        @test Dict(r.id => r.ok for r in ev.phase_rows) == Dict(e1 => true, e2 => false)
        ev = HimalayaUI.peak_eval(exps, picks; every.targets, every.label_ids)
        @test ev.summary.phase == (1, 3)
        DBInterface.close!(ro)
    end
end

@testset "phase targets keep the full call list of a partly reachable exposure" begin
    mktempdir() do dir
        b = _pt_db(dir; experiments = 1, samples = 1, per = 2)
        e1, e2 = b.eids
        member = "(SELECT index_id FROM assignment_members WHERE exposure_id = ?)"
        q0 = only(r.basis for r in HimalayaUI.sqlrows(DBInterface.execute(b.db,
            "SELECT basis FROM indices WHERE id = $member", [e1])))
        DBInterface.execute(b.db, "UPDATE indices SET phase = 'Hexagonal' WHERE id = $member", [e1])
        ix = Int(DBInterface.lastrowid(DBInterface.execute(b.db,
            "INSERT INTO indices (exposure_id, phase, basis, kind) VALUES (?, 'Lamellar', ?, 'auto')", [e1, q0])))
        DBInterface.execute(b.db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [e1, ix])
        DBInterface.close!(b.db)
        ro = HimalayaUI.open_readonly_db(b.path)
        exps = HimalayaUI.labeled_exposures(ro; include_unlabelled = true, io = devnull)

        (; targets, counts) = HimalayaUI.eval_context(ro, exps)
        @test [p for (p, _) in targets[e1]] == ["Hexagonal", "Lamellar"] && length(targets[e2]) == 1
        @test counts.all == 2 && counts.reachable == 2
        @test counts.per_phase["Hexagonal"] == (reachable = 1, all = 1)
        @test counts.per_phase["Lamellar"] == (reachable = 1, all = 1)
        ev = HimalayaUI.peak_eval(exps, HimalayaUI.oracle_picks(exps); targets, label_ids = Set(b.eids))
        rows = Dict(r.id => r for r in ev.phase_rows)
        @test rows[e1].phase == "Hexagonal" && rows[e1].ok && rows[e2].ok
        every = HimalayaUI.eval_context(ro, exps; phase_targets = :all)
        @test every.targets == targets
        DBInterface.close!(ro)
    end
end

@testset "a detector's peaks can make a call reachable that the user's peaks cannot" begin
    mktempdir() do dir
        b = _pt_db(dir; experiments = 1, samples = 1, per = 2)
        e1, e2 = b.eids
        exps = HimalayaUI.labeled_exposures(b.db; include_unlabelled = true, io = devnull)
        e = exps[1]
        q0 = 0.06
        pk = unique(argmin(abs.(e.q .- q0 * r)) for r in (1, sqrt(3), 2, sqrt(7), 3))
        sh = HimalayaUI._sharpness_at(HimalayaUI.Himalaya.sharpness(e.I), pk)
        hex = [c for c in HimalayaUI.Himalaya.indexpeaks(e.q[pk], sh) if nameof(HimalayaUI.Himalaya.phase(c)) == :Hexagonal]
        @test !isempty(hex)
        basis = HimalayaUI.Himalaya.basis(first(hex))
        DBInterface.execute(b.db, "UPDATE indices SET phase = 'Hexagonal', basis = ?
            WHERE id = (SELECT index_id FROM assignment_members WHERE exposure_id = ?)", [basis, e1])
        DBInterface.close!(b.db)
        ro = HimalayaUI.open_readonly_db(b.path)
        exps = HimalayaUI.labeled_exposures(ro; include_unlabelled = true, io = devnull)

        @test Set(keys(HimalayaUI.eval_context(ro, exps).targets)) == Set([e2])
        picks = [i == 1 ? pk : Int[] for i in eachindex(exps)]
        ctx = HimalayaUI.eval_context(ro, exps; detector_picks = [picks])
        @test Set(keys(ctx.targets)) == Set(b.eids) && ctx.counts.reachable == 2
        ev = HimalayaUI.peak_eval(exps, picks; ctx.targets, ctx.label_ids)
        @test Dict(r.id => r.ok for r in ev.phase_rows)[e1]
        DBInterface.close!(ro)
    end
end
